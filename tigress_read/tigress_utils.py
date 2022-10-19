import os
import shutil

def copy_file(fn_old, fn_new):
    """Copy file from old to new location. Make new directory if needed.
    Raise warning if old file does not exist.
    inputs:
        fn_old: string, old filename
        fn_new: string, new filename"""
    if fn_old == fn_new:
        raise ValueError('New filename cannot be the same as the original one.')
    if os.path.isfile(fn_old):
        if os.path.isfile(fn_new):
            print("New file already exists: "+fn_new)
            while True:
                answer = input("Overwrite? (y/n):")
                if answer == "y":
                    break
                elif answer == "n":
                    return          
        dir_new = os.path.dirname(fn_new)
        if not os.path.exists(dir_new):
            os.makedirs(dir_new)
        shutil.copyfile(fn_old, fn_new)
    else:
        raise RuntimeError("file does not exist: "+fn_old)
    return
