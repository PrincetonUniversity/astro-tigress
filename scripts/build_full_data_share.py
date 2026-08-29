#!/usr/bin/env python3
"""Build and validate the raw-only TIGRESS symlink share tree."""

import argparse
import json
import os
import stat
import sys
from collections import namedtuple
from datetime import datetime, timezone
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
TARGET = Path(
    "/tigerdata/EOSTRIKE/TIGRESS-classic/TIGRESS-full-data-share"
)
README_TEMPLATE = REPOSITORY / "FULL_CADENCE_SHARE_README.md"
CONTENT_VERSION = "2026-08-29-v1"


Model = namedtuple(
    "Model", ("shared_name", "run_name", "source", "first", "last", "ranks")
)
Entry = namedtuple("Entry", ("model", "source", "relative"))


MODELS = (
    Model(
        "R8_2pc",
        "R8_2pc_rst",
        Path(
            "/tigerdata/EOSTRIKE/TIGRESS-classic/"
            "TIGRESS-R8/R8_2pc_rst"
        ),
        285,
        448,
        56,
    ),
    Model(
        "R8_4pc",
        "R8_4pc_newacc",
        Path(
            "/tigerdata/EOSTRIKE/TIGRESS-classic/"
            "TIGRESS-local/R8_4pc_newacc"
        ),
        200,
        650,
        28,
    ),
)


def iter_entries(model: Model):
    """Yield the exact allowlisted source and export path pairs."""
    yield Entry(
        model,
        model.source / f"{model.run_name}.par",
        Path(f"{model.run_name}.par"),
    )
    for suffix in ("hst", "sn"):
        relative = Path("hst") / f"{model.run_name}.{suffix}"
        yield Entry(model, model.source / relative, relative)

    for number in range(model.first, model.last + 1):
        snapshot = f"{number:04d}"
        relative = Path("starpar") / (
            f"{model.run_name}.{snapshot}.starpar.vtk"
        )
        yield Entry(model, model.source / relative, relative)

        for rank in range(model.ranks):
            directory = f"id{rank}"
            if rank == 0:
                filename = f"{model.run_name}.{snapshot}.vtk"
            else:
                filename = f"{model.run_name}-id{rank}.{snapshot}.vtk"
            relative = Path(directory) / filename
            yield Entry(model, model.source / relative, relative)


def require_regular_file(path: Path) -> os.stat_result:
    """Return lstat data after rejecting missing, linked, or special files."""
    try:
        metadata = path.lstat()
    except FileNotFoundError as error:
        raise RuntimeError(f"Required source file is missing: {path}") from error
    if not stat.S_ISREG(metadata.st_mode):
        raise RuntimeError(f"Required path is not a regular file: {path}")
    return metadata


def inventory():
    """Validate sources and return entries plus per-model summary data."""
    parent_device = TARGET.parent.stat().st_dev
    entries_by_model = {}
    summaries = {}

    for model in MODELS:
        if model.source.stat().st_dev != parent_device:
            raise RuntimeError(
                f"Source and target parent are on different filesystems: "
                f"{model.source}"
            )

        entries = list(iter_entries(model))
        seen = set()
        total_bytes = 0
        for entry in entries:
            if entry.relative in seen:
                raise RuntimeError(
                    f"Duplicate export path for {model.shared_name}: "
                    f"{entry.relative}"
                )
            seen.add(entry.relative)
            total_bytes += require_regular_file(entry.source).st_size

        expected = 3 + (model.last - model.first + 1) * (model.ranks + 1)
        if len(entries) != expected:
            raise RuntimeError(
                f"Internal count error for {model.shared_name}: "
                f"{len(entries)} != {expected}"
            )

        entries_by_model[model.shared_name] = entries
        summaries[model.shared_name] = {
            "run_name": model.run_name,
            "first_snapshot": model.first,
            "last_snapshot": model.last,
            "snapshot_count": model.last - model.first + 1,
            "processor_count": model.ranks,
            "symlink_count": len(entries),
            "apparent_bytes": total_bytes,
        }

    return entries_by_model, summaries


def utc_mtime(metadata: os.stat_result) -> str:
    return datetime.fromtimestamp(
        metadata.st_mtime, tz=timezone.utc
    ).isoformat()


def write_manifest(root: Path, entries) -> None:
    manifest = root / "MANIFEST.tsv"
    with manifest.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("path\tsize_bytes\tmtime_utc\tdevice\tinode\n")
        for entry in entries:
            metadata = entry.source.stat()
            stream.write(
                f"{entry.relative.as_posix()}\t{metadata.st_size}\t"
                f"{utc_mtime(metadata)}\t{metadata.st_dev}\t"
                f"{metadata.st_ino}\n"
            )


def validate(root: Path, entries_by_model) -> None:
    """Validate exact content and symlink targets of a built tree."""
    if not root.is_dir() or root.is_symlink():
        raise RuntimeError(f"Share root is not a real directory: {root}")

    expected_root = {"README.md", "RELEASE.json"}
    expected_root.update(model.shared_name for model in MODELS)
    actual_root = {path.name for path in root.iterdir()}
    if actual_root != expected_root:
        raise RuntimeError(
            f"Unexpected root content: {sorted(actual_root ^ expected_root)}"
        )
    for filename in ("README.md", "RELEASE.json"):
        require_regular_file(root / filename)

    for model in MODELS:
        model_root = root / model.shared_name
        if not model_root.is_dir() or model_root.is_symlink():
            raise RuntimeError(f"Model root is not a real directory: {model_root}")
        entries = entries_by_model[model.shared_name]
        expected_links = {entry.relative for entry in entries}
        expected_files = {Path("MANIFEST.tsv")}
        expected_directories = set()
        for relative in expected_links | expected_files:
            expected_directories.update(relative.parents)
        expected_directories.discard(Path("."))

        for path in model_root.rglob("*"):
            relative = path.relative_to(model_root)
            metadata = path.lstat()
            if stat.S_ISLNK(metadata.st_mode):
                if relative not in expected_links:
                    raise RuntimeError(f"Unexpected export symlink: {path}")
            if stat.S_ISDIR(metadata.st_mode):
                if relative not in expected_directories:
                    raise RuntimeError(f"Unexpected export directory: {path}")
            elif stat.S_ISREG(metadata.st_mode):
                if relative not in expected_files:
                    raise RuntimeError(f"Unexpected export file: {path}")
            elif stat.S_ISLNK(metadata.st_mode):
                pass
            else:
                raise RuntimeError(f"Unexpected special export path: {path}")

        for entry in entries:
            destination = model_root / entry.relative
            source_stat = require_regular_file(entry.source)
            if not destination.is_symlink():
                raise RuntimeError(
                    f"Export is not a symlink: {destination}"
                )
            if Path(os.readlink(str(destination))) != entry.source:
                raise RuntimeError(f"Incorrect symlink target: {destination}")
            destination_stat = destination.stat()
            if (
                source_stat.st_dev,
                source_stat.st_ino,
                source_stat.st_size,
            ) != (
                destination_stat.st_dev,
                destination_stat.st_ino,
                destination_stat.st_size,
            ):
                raise RuntimeError(f"Symlink resolves incorrectly: {destination}")

        manifest = model_root / "MANIFEST.tsv"
        line_count = sum(1 for _ in manifest.open(encoding="utf-8"))
        if line_count != len(entries) + 1:
            raise RuntimeError(
                f"Manifest count mismatch for {model.shared_name}: "
                f"{line_count - 1} != {len(entries)}"
            )


def build(entries_by_model, summaries) -> None:
    """Build in a unique staging directory, validate, then rename."""
    if TARGET.exists() or TARGET.is_symlink():
        raise RuntimeError(f"Refusing to overwrite existing target: {TARGET}")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    staging = TARGET.with_name(f"{TARGET.name}.staging-{stamp}-{os.getpid()}")
    if staging.exists() or staging.is_symlink():
        raise RuntimeError(f"Staging path already exists: {staging}")

    staging.mkdir(mode=0o750)
    try:
        for model in MODELS:
            model_root = staging / model.shared_name
            model_root.mkdir()
            entries = entries_by_model[model.shared_name]
            for entry in entries:
                destination = model_root / entry.relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                os.symlink(str(entry.source), str(destination))
            write_manifest(model_root, entries)

        readme = README_TEMPLATE.read_text(encoding="utf-8")
        (staging / "README.md").write_text(readme, encoding="utf-8")
        release = {
            "content_version": CONTENT_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "layout": "native-processor-file-symlinks",
            "target": str(TARGET),
            "models": summaries,
        }
        (staging / "RELEASE.json").write_text(
            json.dumps(release, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        validate(staging, entries_by_model)
        staging.rename(TARGET)
    except Exception:
        print(
            f"Build failed; staging was retained for inspection: {staging}",
            file=sys.stderr,
        )
        raise

    print(f"Created and validated {TARGET}")


def print_summary(summaries) -> None:
    for name, summary in summaries.items():
        gib = summary["apparent_bytes"] / 1024**3
        print(
            f"{name}: {summary['first_snapshot']:04d}-"
            f"{summary['last_snapshot']:04d}, "
            f"{summary['snapshot_count']} snapshots, "
            f"{summary['symlink_count']} symlinks, {gib:.2f} GiB apparent"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--execute", action="store_true", help="build and publish the tree"
    )
    action.add_argument(
        "--validate-only", action="store_true", help="validate the target"
    )
    args = parser.parse_args()

    entries_by_model, summaries = inventory()
    print_summary(summaries)
    if args.execute:
        build(entries_by_model, summaries)
    elif args.validate_only:
        validate(TARGET, entries_by_model)
        print(f"Validated {TARGET}")
    else:
        print("Dry run only; pass --execute to create the share tree.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
