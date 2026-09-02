#!/usr/bin/env python3
"""Migrate the TIGRESS full-data share from symlinks to real data in place.

Transformation per model:

* VTK files (``id*/*.vtk``) are MOVED from the native run directory into the
  share tree with ``os.rename`` (a server-side rename on the same filesystem --
  no bytes are copied, no extra space is used). A symbolic link is then left at
  the native path so local tooling (pyathena) still resolves the original path.
* ``starpar/*.starpar.vtk``, ``hst/<run>.hst``, ``hst/<run>.sn`` and
  ``<run>.par`` are COPIED into the share tree. The native originals are kept.

After this runs, the share tree holds the only physical copy of each VTK file,
and the native ``id*`` directories hold symlinks pointing into the share tree.

The default action is a DRY RUN: it classifies every file and prints what would
happen without touching anything. Pass ``--execute`` to act. The operation is
idempotent and resumable -- re-running skips files already in the target state
and repairs any that were interrupted mid-step. ``--validate`` audits the final
state; ``--finalize`` regenerates ``MANIFEST.tsv`` and ``RELEASE.json``.

Nothing is ever copied over, or renamed onto, an existing *real* file: the only
thing at a share destination that this script will replace is a symlink (the old
export) or a size-mismatched partial copy. Native VTK data is only ever removed
by an atomic rename that simultaneously creates it at the share path, so there
is no moment at which the bytes do not exist somewhere.
"""

import argparse
import json
import os
import re
import shutil
import stat
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


SHARE_ROOT = Path("/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share")
REPOSITORY = Path(__file__).resolve().parents[1]
README_TEMPLATE = REPOSITORY / "FULL_CADENCE_SHARE_README.md"
CONTENT_VERSION = "2026-09-01-realdata-v1"
LAYOUT = "native-vtk-moved-real-files"

VTK_NAME = re.compile(r"^(?P<run>.+?)(?:-id(?P<rank>\d+))?\.(?P<snap>\d{4})\.vtk$")
ID_DIR = re.compile(r"^id\d+$")


class Model:
    def __init__(self, shared_name, run_name, native, ranks):
        self.shared_name = shared_name
        self.run_name = run_name
        self.native = Path(native)
        self.ranks = ranks

    @property
    def share(self):
        return SHARE_ROOT / self.shared_name


MODELS = [
    Model(
        "R8_2pc",
        "R8_2pc_rst",
        "/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-R8/R8_2pc_rst",
        56,
    ),
    Model(
        "R8_4pc",
        "R8_4pc_newacc",
        "/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-local/R8_4pc_newacc",
        28,
    ),
]


# ---------------------------------------------------------------------------
# discovery
# ---------------------------------------------------------------------------

class Entry:
    __slots__ = ("op", "native", "share", "relative", "model")

    def __init__(self, op, native, share, relative, model):
        self.op = op            # "move" or "copy"
        self.native = native    # Path in the native run directory
        self.share = share      # Path in the share tree
        self.relative = relative
        self.model = model


def discover(model, snap_floor=None, snap_ceil=None, warnings=None):
    """Return the ordered list of Entry objects for one model."""
    entries = []
    warnings = warnings if warnings is not None else []

    # VTK -> move
    id_dirs = sorted(
        (p for p in model.native.iterdir() if ID_DIR.match(p.name) and p.is_dir()),
        key=lambda p: int(p.name[2:]),
    )
    for id_dir in id_dirs:
        for path in sorted(id_dir.iterdir()):
            if not path.name.endswith(".vtk"):
                continue
            match = VTK_NAME.match(path.name)
            if not match:
                warnings.append(f"skip non-standard vtk: {path}")
                continue
            snap = int(match.group("snap"))
            if snap_floor is not None and snap < snap_floor:
                continue
            if snap_ceil is not None and snap > snap_ceil:
                continue
            relative = Path(id_dir.name) / path.name
            entries.append(Entry("move", path, model.share / relative, relative, model))

    # starpar -> copy
    starpar = model.native / "starpar"
    if starpar.is_dir():
        for path in sorted(starpar.iterdir()):
            if not path.name.endswith(".starpar.vtk"):
                continue
            match = VTK_NAME.match(path.name.replace(".starpar.vtk", ".vtk"))
            if match:
                snap = int(match.group("snap"))
                if snap_floor is not None and snap < snap_floor:
                    continue
                if snap_ceil is not None and snap > snap_ceil:
                    continue
            relative = Path("starpar") / path.name
            entries.append(Entry("copy", path, model.share / relative, relative, model))

    # hst files -> copy
    for suffix in ("hst", "sn"):
        native = model.native / "hst" / f"{model.run_name}.{suffix}"
        relative = Path("hst") / native.name
        entries.append(Entry("copy", native, model.share / relative, relative, model))

    # par -> copy
    native = model.native / f"{model.run_name}.par"
    relative = Path(native.name)
    entries.append(Entry("copy", native, model.share / relative, relative, model))

    return entries


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------

def lstat_or_none(path):
    try:
        return os.lstat(path)
    except FileNotFoundError:
        return None


def _kind(st):
    if st is None:
        return "missing"
    if stat.S_ISLNK(st.st_mode):
        return "symlink"
    if stat.S_ISREG(st.st_mode):
        return "file"
    if stat.S_ISDIR(st.st_mode):
        return "dir"
    return "special"


def classify(entry):
    """Return (state, detail). state in {done, move, copy, relink, error}."""
    n = lstat_or_none(entry.native)
    s = lstat_or_none(entry.share)
    nk, sk = _kind(n), _kind(s)

    if entry.op == "move":
        # target: native symlink -> share real file
        if nk == "symlink" and sk == "file":
            if os.readlink(entry.native) == str(entry.share):
                return "done", None
            return "error", f"native symlink -> {os.readlink(entry.native)} (expected {entry.share})"
        if nk == "missing" and sk == "file":
            return "relink", "native missing, share present -> recreate native symlink"
        if nk == "file" and sk in ("missing", "symlink"):
            return "move", None
        if nk == "file" and sk == "file":
            return "error", "both native and share are real files (conflict; refusing)"
        if nk == "missing" and sk == "missing":
            return "error", "both native and share are missing (data absent!)"
        return "error", f"unexpected states native={nk} share={sk}"

    # copy
    if nk != "file":
        return "error", f"native is {nk}, expected a real file to copy"
    if sk == "file" and n.st_size == s.st_size:
        return "done", None
    return "copy", None


# ---------------------------------------------------------------------------
# actions
# ---------------------------------------------------------------------------

def do_move(entry, log):
    share = entry.share
    share.parent.mkdir(parents=True, exist_ok=True)
    s = lstat_or_none(share)
    if s is not None and not stat.S_ISLNK(s.st_mode):
        raise RuntimeError(f"refusing move: share exists and is not a symlink: {share}")
    # atomic rename: replaces an existing symlink (or nothing) with the real file.
    os.rename(entry.native, share)          # raises OSError(EXDEV) rather than copying
    if os.path.lexists(entry.native):
        raise RuntimeError(f"native still present after rename: {entry.native}")
    os.symlink(str(share), entry.native)
    _verify_move(entry)
    log(f"move  {entry.model.shared_name}/{entry.relative}")


def do_relink(entry, log):
    if os.path.lexists(entry.native):
        raise RuntimeError(f"native unexpectedly present: {entry.native}")
    if not stat.S_ISREG(os.lstat(entry.share).st_mode):
        raise RuntimeError(f"share is not a real file, cannot relink: {entry.share}")
    os.symlink(str(entry.share), entry.native)
    _verify_move(entry)
    log(f"relink {entry.model.shared_name}/{entry.relative}")


def _verify_move(entry):
    ns = os.lstat(entry.native)
    ss = os.lstat(entry.share)
    if not stat.S_ISLNK(ns.st_mode):
        raise RuntimeError(f"native is not a symlink after move: {entry.native}")
    if not stat.S_ISREG(ss.st_mode):
        raise RuntimeError(f"share is not a real file after move: {entry.share}")
    if os.readlink(entry.native) != str(entry.share):
        raise RuntimeError(f"native symlink target wrong: {entry.native}")
    if os.stat(entry.native).st_size != ss.st_size:
        raise RuntimeError(f"native symlink resolves to wrong size: {entry.native}")


def do_copy(entry, log):
    share = entry.share
    share.parent.mkdir(parents=True, exist_ok=True)
    tmp = share.with_name(f".{share.name}.part-{os.getpid()}")
    with open(entry.native, "rb") as src, open(tmp, "wb") as dst:
        shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
        dst.flush()
        os.fsync(dst.fileno())
    shutil.copystat(entry.native, tmp)
    os.replace(tmp, share)                  # atomic; replaces symlink or partial file
    src_size = os.lstat(entry.native).st_size
    dst_size = os.lstat(share).st_size
    if src_size != dst_size:
        raise RuntimeError(f"copy size mismatch for {share}: {dst_size} != {src_size}")
    log(f"copy  {entry.model.shared_name}/{entry.relative}")


# ---------------------------------------------------------------------------
# preflight self-tests (execute only)
# ---------------------------------------------------------------------------

def preflight(model):
    model.share.mkdir(parents=True, exist_ok=True)
    tag = f".migrate_selftest_{os.getpid()}"
    # rename native -> share must be a true rename (same filesystem)
    src = model.native / tag
    dst = model.share / tag
    src.write_bytes(b"selftest")
    try:
        os.rename(src, dst)
    except OSError as error:
        src.unlink(missing_ok=True)
        raise RuntimeError(
            f"[{model.shared_name}] rename native->share failed ({error}); "
            f"cannot migrate without a real copy. Aborting."
        )
    dst.unlink(missing_ok=True)
    # symlink creation in a native id directory
    link = model.native / "id0" / tag
    os.symlink("/dev/null", link)
    link.unlink()


# ---------------------------------------------------------------------------
# manifest / release
# ---------------------------------------------------------------------------

def utc(ts):
    return datetime.fromtimestamp(ts, tz=timezone.utc).isoformat()


def write_manifest(model, entries):
    manifest = model.share / "MANIFEST.tsv"
    with manifest.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("path\top\tsize_bytes\tmtime_utc\n")
        for entry in entries:
            st = os.stat(entry.share)
            stream.write(
                f"{entry.relative.as_posix()}\t{entry.op}\t{st.st_size}\t{utc(st.st_mtime)}\n"
            )


def write_release(all_entries):
    models = {}
    for model in MODELS:
        rows = all_entries[model.shared_name]
        moved = [e for e in rows if e.op == "move"]
        copied = [e for e in rows if e.op == "copy"]
        snaps = sorted({int(VTK_NAME.match(e.native.name).group("snap")) for e in moved})
        models[model.shared_name] = {
            "run_name": model.run_name,
            "vtk_files_moved": len(moved),
            "copied_files": len(copied),
            "first_snapshot": snaps[0] if snaps else None,
            "last_snapshot": snaps[-1] if snaps else None,
            "snapshot_count": len(snaps),
            "moved_bytes": sum(os.stat(e.share).st_size for e in moved),
        }
    release = {
        "content_version": CONTENT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "layout": LAYOUT,
        "collection": {
            "uuid": "db422b43-55bc-41d6-876f-ffe7e6582b47",
            "root_filesystem_path": "/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share",
            "shared_path": "/",
            "access": "read on /",
            "web_app_url": (
                "https://app.globus.org/file-manager?origin_id="
                "db422b43-55bc-41d6-876f-ffe7e6582b47&origin_path=%2F"
            ),
        },
        "target": str(SHARE_ROOT),
        "models": models,
    }
    (SHARE_ROOT / "RELEASE.json").write_text(
        json.dumps(release, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------

def validate(all_entries):
    problems = 0
    for model in MODELS:
        for entry in all_entries[model.shared_name]:
            state, detail = classify(entry)
            if entry.op == "move" and state != "done":
                print(f"  FAIL move {entry.relative}: {state} {detail or ''}")
                problems += 1
            if entry.op == "copy" and state != "done":
                print(f"  FAIL copy {entry.relative}: {state} {detail or ''}")
                problems += 1
    if problems:
        print(f"validation FAILED: {problems} problem(s)")
        return False
    print("validation OK: every VTK is a native symlink to a real share file; "
          "every copy present with matching size.")
    return True


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def gather(args, warnings):
    return {
        m.shared_name: discover(m, args.snap_floor, args.snap_ceil, warnings)
        for m in MODELS
    }


def summarize(all_entries, warnings):
    print(f"{'model':8s} {'move':>7s} {'copy':>6s} {'done':>7s} {'relink':>7s} "
          f"{'error':>6s} {'move_TiB':>9s} {'copy_MiB':>9s}")
    grand = dict(move=0, copy=0, done=0, relink=0, error=0, mb=0, tb=0.0)
    errors = []
    for model in MODELS:
        c = dict(move=0, copy=0, done=0, relink=0, error=0)
        move_bytes = 0
        copy_bytes = 0
        for entry in all_entries[model.shared_name]:
            state, detail = classify(entry)
            if state == "error":
                c["error"] += 1
                errors.append((entry, detail))
                continue
            c[state] += 1
            if state == "move":
                move_bytes += os.lstat(entry.native).st_size
            elif state == "copy":
                copy_bytes += os.lstat(entry.native).st_size
        print(f"{model.shared_name:8s} {c['move']:7d} {c['copy']:6d} {c['done']:7d} "
              f"{c['relink']:7d} {c['error']:6d} {move_bytes/1024**4:9.2f} "
              f"{copy_bytes/1024**2:9.2f}")
        for k in ("move", "copy", "done", "relink", "error"):
            grand[k] += c[k]
        grand["tb"] += move_bytes / 1024**4
        grand["mb"] += copy_bytes / 1024**2
    print(f"{'TOTAL':8s} {grand['move']:7d} {grand['copy']:6d} {grand['done']:7d} "
          f"{grand['relink']:7d} {grand['error']:6d} {grand['tb']:9.2f} {grand['mb']:9.2f}")
    if warnings:
        print(f"\n{len(warnings)} warning(s):")
        for w in warnings[:20]:
            print(f"  ! {w}")
    if errors:
        print(f"\n{len(errors)} ERROR state(s) -- execute is blocked until resolved:")
        for entry, detail in errors[:40]:
            print(f"  x {entry.model.shared_name}/{entry.relative}: {detail}")
    return grand["error"] == 0


def run_execute(all_entries):
    log_path = SHARE_ROOT.parent / f".migrate-{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.log"
    stream = open(log_path, "a", encoding="utf-8")

    def log(message):
        stream.write(f"{datetime.now(timezone.utc).isoformat()}\t{message}\n")
        stream.flush()

    print(f"logging to {log_path}")
    for model in MODELS:
        print(f"preflight {model.shared_name} ...")
        preflight(model)

    done = moved = copied = relinked = 0
    t0 = time.time()
    for model in MODELS:
        entries = all_entries[model.shared_name]
        for i, entry in enumerate(entries):
            state, detail = classify(entry)
            if state == "error":
                raise RuntimeError(f"{entry.model.shared_name}/{entry.relative}: {detail}")
            if state == "done":
                done += 1
            elif state == "move":
                do_move(entry, log); moved += 1
            elif state == "relink":
                do_relink(entry, log); relinked += 1
            elif state == "copy":
                do_copy(entry, log); copied += 1
            if (moved + copied + relinked) and (moved + copied + relinked) % 500 == 0:
                rate = (moved + copied + relinked) / max(time.time() - t0, 1e-9)
                print(f"  {model.shared_name}: {i+1}/{len(entries)}  "
                      f"moved={moved} copied={copied} relinked={relinked}  "
                      f"({rate:.0f}/s)")
        write_manifest(model, entries)
    write_release(all_entries)
    stream.close()
    print(f"done: moved={moved} copied={copied} relinked={relinked} already_done={done}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--execute", action="store_true", help="perform the migration")
    group.add_argument("--validate", action="store_true", help="audit the final state")
    group.add_argument("--finalize", action="store_true",
                       help="regenerate MANIFEST.tsv and RELEASE.json only")
    parser.add_argument("--snap-floor", type=int, default=None,
                        help="only include snapshots >= this number")
    parser.add_argument("--snap-ceil", type=int, default=None,
                        help="only include snapshots <= this number")
    args = parser.parse_args()

    warnings = []
    all_entries = gather(args, warnings)

    if args.validate:
        return 0 if validate(all_entries) else 1
    if args.finalize:
        for model in MODELS:
            write_manifest(model, all_entries[model.shared_name])
        write_release(all_entries)
        print("regenerated MANIFEST.tsv and RELEASE.json")
        return 0

    ok = summarize(all_entries, warnings)
    if args.execute:
        if not ok:
            print("\nrefusing to execute while ERROR states exist.")
            return 1
        run_execute(all_entries)
        print("\nvalidating ...")
        return 0 if validate(all_entries) else 1
    print("\nDRY RUN only. Re-run with --execute to perform the migration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
