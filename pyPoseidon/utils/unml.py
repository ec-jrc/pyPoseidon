def unml(nml,patch):

    for g in nml.keys():
        try:
            for (v,r) in patch.items():
                if v in nml[g].keys(): 
                    nml[g][v] = r
        except:
            pass

    return nml
