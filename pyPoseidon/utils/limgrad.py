#LIMGRAD impose gradient limits on a discrete mesh-size function defined over a 2-simplex triangulation (from JIGSAW tools)
# TODO revisit and check
import numpy as np

    
def limgrad2(edge,elen,ffun,dfdx,imax):
    
    rfun = ffun.T.flatten()
    
    rfun = np.array([[hf] for hf in list(rfun)])
    
    eps = np.finfo(float).eps
    
    nnod = rfun.size
        
    #-- IVEC(NPTR(II,1):NPTR(II,2)) are edges adj. to II-TH node
    nvec = np.hstack([edge[:,0], edge[:,1]])

    ivec = np.hstack([np.arange(edge.shape[0]),np.arange(edge.shape[0])])


    nvec_ = np.sort(nvec, kind='mergesort')
    pidx = np.argsort(nvec, kind='mergesort') # to match with matlab/octave -> https://stackoverflow.com/questions/39484073/matlab-sort-vs-numpy-argsort-how-to-match-results
    ivec = ivec[pidx]
    nvec = nvec_
    
    
    mark = np.full(rfun.size, False, dtype=bool)
    mark[edge[:,0]]=True
    mark[edge[:,1]]=True

    dif = [nvec[i+1]-nvec[i] for i in range(nvec.shape[0]-1)]
    idxx = np.argwhere(np.array(dif) > 0).flatten()
    
    nptr=np.zeros((mark.size,2))
    nptr[mark,0] = np.append(np.array([0]),idxx+1)
    nptr[mark,1] = np.append(idxx, nnod-1)
    
    nptr = nptr.astype(int)

#----------------------------- ASET=ITER if node is "active"
    aset = np.zeros(nnod)
    
    ftol = min(rfun.flatten()) * np.sqrt(eps)
    

#----------------------------- exhaustive 'til all satisfied 
    
    for i in range(1,imax):
    
    #------------------------- find "active" nodes this pass
        aidx = np.argwhere(aset == i - 1 ) 
        aidx = aidx.flatten()
        
        if not aidx.any(): break
      
    #------------------------- reorder => better convergence

        aval = np.sort(rfun.reshape(ffun.shape).flatten()[aidx],kind='mergesort')
        idxx = np.argsort(rfun.reshape(ffun.shape).flatten()[aidx], kind='mergesort')
        
        aidx = aidx[idxx]
       
    #%------------------------- visit adj. edges and set DFDX
        for ipos in range(len(aidx)):
            npos = aidx[ipos]
        
            for jpos in range(nptr[npos,0], nptr[npos,1]+1):
                
                epos = ivec[jpos]
                
                nod1 = edge[epos,0]
                nod2 = edge[epos,1]
                
               # print ipos, jpos, epos, nod1, nod2

            #----------------- calc. limits about min.-value
                if rfun[nod1] > rfun[nod2]:
                
                    
                    fun1 = rfun[nod2] + elen[epos] * dfdx 
                
                    if rfun[nod1] > fun1+ftol :
                        rfun[nod1] = fun1
                        aset[nod1] = i
                else:
                
                    fun2 = rfun[nod1] + elen[epos] * dfdx 
                    
                    if   rfun[nod2] > fun2+ftol :
                        rfun[nod2] = fun2
                        aset[nod2] = i
     
    flag = i < imax
    
    return rfun,flag
    
