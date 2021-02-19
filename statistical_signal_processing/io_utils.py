import tables

def read_whole_h5file(fname):
    
    h5file=tables.open_file(fname)
    output={}
    for group in h5file.walk_groups("/"):
        output[group._v_pathname]={}
        for array in h5file.list_nodes(group, classname = 'Array'):
            output[group._v_pathname][array.name]=array.read()
    h5file.close()
     
    return output
    
def read_partial_h5file(fname,groupstodo):
    
    h5file=tables.open_file(fname)
    output={}
    for group in h5file.walk_groups("/"):
        if (group._v_pathname) in groupstodo:
            output[group._v_pathname]={}
            for array in h5file.list_nodes(group, classname = 'Array'):
                output[group._v_pathname][array.name]=array.read()
    h5file.close()
     
    return output    
    
