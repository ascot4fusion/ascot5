if __name__ == "__main__":
    if len(sys.argv) == 4:
        field, subfield = sys.argv[2].split('/')
        copyinput(sys.argv[1], sys.argv[2], field, subfield)
    else:
        print 'This function should be called as:'
        print 'python copyinput.py source.h5 target.h5 copied/field'

def copyinput(fns,fnt,field,subfield):
    """Copy input field from one hdf5 file into another

    Keyword arguments:
    fns      -- path to source hdf5 file
    fnt      -- path to target hdf5 file
    field    -- main field where copying is done e.g. bfield
    subfield -- subfield to be copied e.g. B2D
    """


    # Get the target field and type from source
    fs = h5py.File(fns, "r")
    o = fs[field]
    types = o.attrs["type"]
    
 if not "/bfield" in f:
        o = f.create_group('bfield')

    # Check if the field in target exists (otherwise it is created)
    # Delete possible existing types and sub fields of same type as copied
    ft = h5py.File(fnt, "a")
    if not field in ft:
        ot = ft.create_group(field)
    else:
        ot = ft[field]
        del ot.attrs["type"]
        if subfield in ot:
            del ot[subfield]

    # Do the copying and set the type
    fs.copy(field + '/' + subfield, ot, name=subfield)
    ot.attrs["type"] = np.string_(types)

    fs.close()
    ft.close()
