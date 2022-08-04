"""
Remove fields from a vtk file
"""

import read_vtk
import numpy as np
import os

def delete_fields_vtk(fname, fname_new, field_del, verbose=False):
    """
    Function to re-write a vtk file with removal of specified fields.
    Works only for joined vtk files.

    Parameters
    ----------
    fname : str
        Name of the original vtk file
    fname_new : str
        Name of the new vtk file
    field_del : (list of) str
        Field names to be removed
    verbose : bool
        Produce verbose messages if True.

    Returns nothing
    """
    
    if fname == fname_new:
        raise ValueError('New filename cannot be the same as the original one.')
    if os.path.isfile(fname):
        if os.path.isfile(fname_new):
            print("New file already exists: "+fname_new)
            while True:
                answer = input("Overwrite? (y/n):")
                if answer == "y":
                    break
                elif answer == "n":
                    return          
        dir_new = os.path.dirname(fname_new)
        if not os.path.exists(dir_new):
            os.makedirs(dir_new)
    else:
        raise RuntimeError("file does not exist: "+fname)

    field_del = np.atleast_1d(field_del)
    
    ds = read_vtk.read_vtk(fname)
    fmap = ds.grid[0]['field_map']
    # Remove fields to be deleted from field map
    for f_del in field_del:
        fmap.pop(f_del)

    if verbose:
        print('Original:', fname)
        print('New:', fname_new)
        print('Writing:', end=' ')
    with open(fname, 'rb') as fp, open(fname_new, 'wb') as fp2:
        # Write header
        fp2.write(fp.read(ds.grid[0]['data_offset']))
        for f in fmap:
            if verbose:
                print(f, end=' ')
            fp2.write(fmap[f]['title'])
            #if fmap[f]['title'].decode('ascii').startswith('SCALARS'):
            if fmap[f]['read_table']:
                fp2.write(b'LOOKUP_TABLE default\n')

            fp.seek(fmap[f]['offset'])
            fp.readline()
            if fmap[f]['read_table']:
                fp.readline()
            fp2.write(fp.read(fmap[f]['dsize']))

    return None
