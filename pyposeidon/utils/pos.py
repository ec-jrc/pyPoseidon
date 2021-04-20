import numpy as np
import pandas as pd

def MakeFacesVectorized(Nr,Nc):

    out = np.empty((Nr-1,Nc-1,4),dtype=int)

    r = np.arange(Nr*Nc).reshape(Nr,Nc)

    out[:,:, 0] = r[:-1,:-1]
    out[:,:, 1] = r[:-1,1:]
    out[:,:, 2] = r[1:,1:]
    out[:,:, 3] = r[1:,:-1]

    out.shape =(-1,4)
    return out
    
    
def MakeFacesVectorized_periodic(Nr,Nc):

    out = np.empty((Nr-1,Nc,4),dtype=int)

    r = np.arange(Nr*Nc).reshape(Nr,Nc)

    out[:,:-1, 0] = r[:-1,:-1]
    out[:,:-1, 1] = r[:-1,1:]
    out[:,:-1, 2] = r[1:,1:]
    out[:,:-1, 3] = r[1:,:-1]
    
    out[:,-1, 0] = r[:-1,-1]
    out[:,-1, 1] = r[:-1,0]
    out[:,-1, 2] = r[1:,0]
    out[:,-1, 3] = r[1:,-1]

    

    out.shape =(-1,4)
    return out
    
    
def to_df(elems,nodes):
    # cells to polygons
    ap = nodes.loc[elems.a,['longitude','latitude']]
    bp = nodes.loc[elems.b,['longitude','latitude']]
    cp = nodes.loc[elems.c,['longitude','latitude']]
    dp = nodes.loc[elems.d,['longitude','latitude']]
    
    ap['z']=0
    bp['z']=0
    cp['z']=0
    dp['z']=0
    
    elems['ap'] = ap.values.tolist()
    elems['bp'] = bp.values.tolist()
    elems['cp'] = cp.values.tolist()
    elems['dp'] = dp.values.tolist()
    
    elems['va']=nodes.loc[elems.a,'d2'].values.tolist()
    elems['vb']=nodes.loc[elems.b,'d2'].values.tolist()
    elems['vc']=nodes.loc[elems.c,'d2'].values.tolist()
    elems['vd']=nodes.loc[elems.d,'d2'].values.tolist()
    
    return elems.drop(['a','b','c','d'],axis=1)

    
def to_df_uv(elems,nodes):
    # cells to polygons
    ap = nodes.loc[elems.a,['u','v']]
    bp = nodes.loc[elems.b,['u','v']]
    cp = nodes.loc[elems.c,['u','v']]
    dp = nodes.loc[elems.d,['u','v']]
    
    ap['z']=0
    bp['z']=0
    cp['z']=0
    dp['z']=0
    
    elems['ap'] = ap.values.tolist()
    elems['bp'] = bp.values.tolist()
    elems['cp'] = cp.values.tolist()
    elems['dp'] = dp.values.tolist()
    
    elems['va']=nodes.loc[elems.a,'d2'].values.tolist()
    elems['vb']=nodes.loc[elems.b,'d2'].values.tolist()
    elems['vc']=nodes.loc[elems.c,'d2'].values.tolist()
    elems['vd']=nodes.loc[elems.d,'d2'].values.tolist()
    
    return elems.drop(['a','b','c','d'],axis=1)
    

def to_df_3d(elems,nodes):
    # cells to polygons
    ap = nodes.loc[elems.a,['x','y','z']]
    bp = nodes.loc[elems.b,['x','y','z']]
    cp = nodes.loc[elems.c,['x','y','z']]
    dp = nodes.loc[elems.d,['x','y','z']]
        
    elems['ap'] = ap.values.tolist()
    elems['bp'] = bp.values.tolist()
    elems['cp'] = cp.values.tolist()
    elems['dp'] = dp.values.tolist()
    
    elems['va']=nodes.loc[elems.a,'d2'].values.tolist()
    elems['vb']=nodes.loc[elems.b,'d2'].values.tolist()
    elems['vc']=nodes.loc[elems.c,'d2'].values.tolist()
    elems['vd']=nodes.loc[elems.d,'d2'].values.tolist()
    
    return elems.drop(['a','b','c','d'],axis=1)
    
    
    
def to_sq(df,fpos):
    
    with open(fpos,'w') as f:
        f.write('//*********************************************************************\n')
        f.write('// *\n')
        f.write('// *  pyposeidon\n')
        f.write('// *\n') 
        f.write('// *  Scalar 2D post-processing view\n')
        f.write('// *\n')
        f.write('// *********************************************************************/\n\n')

        f.write('// This view contains a scalar field defined on quads.\n')
        f.write('\n')
        f.write('View "{}" {}\n'.format('bgmesh','{'))
        for idx, vals in df.iterrows():
            f.write('SQ({},{},{},{},{},{},{},{},{},{},{},{}){{{},{},{},{}}};\n'
                    .format(
                            vals.ap[0],vals.ap[1],vals.ap[2],
                            vals.bp[0],vals.bp[1],vals.bp[2],
                            vals.cp[0],vals.cp[1],vals.cp[2],
                            vals.dp[0],vals.dp[1],vals.dp[2],
                            vals.va,vals.vb,vals.vc,vals.vd
                    )
                   )

        f.write('{};\n'.format('}'))