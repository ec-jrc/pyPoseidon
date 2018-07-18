import numpy as np
import pickle
import os
#from pyPoseidon.model import *
from pyPoseidon.utils.vis import *
from pyPoseidon.utils.obs import obs
from pyPoseidon.grid import *
import pyresample
import pandas as pd
import datetime
from pyPoseidon.utils.get_value import get_value
import holoviews as hv
import geoviews as gv
from cartopy import crs
import xarray as xr
import pyresample
import glob
import sys
import subprocess
import scipy.interpolate 

FFWriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264','-pix_fmt','yuv420p'])
   
class data:   
    
    def __init__(self,**kwargs):
    
       solver = kwargs.get('solver',None)
       impl=None
       if solver == 'd3d':
           self.impl = d3d(**kwargs)
       elif solver == 'schism':
           self.impl = schism(**kwargs)  
       else:
           sys.stdout.flush()
           sys.stdout.write('\n')
           sys.stdout.write('solver is not defined, exiting \n')
           sys.stdout.flush()
           sys.exit(1)             
       

class d3d(data):
        
    def __init__(self,**kwargs):
        
        rpath = kwargs.get('rpath','./')
        
        folders = kwargs.get('folders',None)  #[os.path.join(os.path.abspath(loc),name) for name in os.listdir(loc) if os.path.isdir(os.path.join(loc,name))]
        
        if folders :
            self.folders = folders
        else:
            self.folders = [rpath]
        
        #check if many tags present
        ifiles = glob.glob(self.folders[0]+'/*_info.pkl')
                
        if len(ifiles) > 1:
            #--------------------------------------------------------------------- 
              sys.stdout.flush()
              sys.stdout.write('\n')
              sys.stdout.write('more than one configuration, specify tag argument \n')
              sys.stdout.flush()
            #--------------------------------------------------------------------- 
                    
        tag = kwargs.get('tag', None)
        
        if tag :
            ifile = self.folders[0]+'/'+tag+'_info.pkl'
        else:
            ifile = ifiles[0]
        
        #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('reading data based on {} \n'.format(ifile))
        sys.stdout.flush()
        #--------------------------------------------------------------------- 
                       
        with open(ifile, 'r') as f:
                          self.info=pickle.load(f)  
                        
        grid=r2d.read_file(self.folders[0]+'/'+self.info['tag']+'.grd')
            
        deb=np.loadtxt(self.folders[0]+'/'+self.info['tag']+'.dep')
        
        #create mask
            
        d=deb[1:-1,1:-1]
        self.w=d==-999.
        
        b=deb[:-1,:-1]
        b[b==-999.]=np.nan
        
        self.dem = xr.Dataset({'bathymetry': (['latitude', 'longitude'], b)},
                    coords={'longitude': ('longitude', grid.lons[0,:]),   
                            'latitude': ('latitude', grid.lats[:,0])})
        
        self.grid = grid
            
            
            # read array
        
        nfiles = [folder +'/'+'trim-'+self.info['tag']+'.nc' for folder in self.folders]
        
        ds = xr.open_mfdataset(nfiles, autoclose=True)
                
        #grid points
        xg = ds.XCOR.values.T
        yg = ds.YCOR.values.T   
        self.x = xg
        self.y = yg
        
        self.dx = self.x[0,1]-self.x[0,0]
        self.dy = self.y[1,0]-self.y[0,0]        
        
        #pressure points       
        xz = self.x - self.dx/2.
        yz = self.y - self.dy/2.
        
        xz = xz[1:-1,1:-1]  
        yz = yz[1:-1,1:-1]  
        
        self.xh = np.ma.masked_array(xz, self.w) #mask land
        self.yh = np.ma.masked_array(yz, self.w)
        
        #velocity points
        xv = xz + self.dx/2.
        yu = yz + self.dy/2.
                
        self.times=ds.time.values
        
        ovars = kwargs.get('vars', ['S1','U1','V1','WINDU','WINDV','PATM'])
        
        keys = [d for d in ds.variables.keys() if d in ovars]
        
        zkeys = [d for d in keys if d not in ['U1','V1']]
        
        dic = {}  
                
        for key in zkeys:
            h =  ds[key][:,1:-1,1:-1]
            ha = h.transpose(h.dims[0],h.dims[2],h.dims[1]) #transpose lat/lon
            ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
            dic.update({key :(['time','YZ','XZ'], np.ma.masked_where(ww==True,ha))}) #mask land
        
        if 'U1' in keys:        
            # Get water velocities - U1
            h =  ds['U1'][:,0,1:-1,1:-1] 
            ha = h.transpose(h.dims[0],h.dims[2],h.dims[1]) #transpose lat/lon
            ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
            wu = np.ma.masked_where(ww==True,ha) #mask land
            dic.update({'U1' :(['time','YU','XU'], wu)}) 
        
            #interpolate on XZ, YZ
            orig = pyresample.geometry.SwathDefinition(lons=xz,lats=yu) # original points
            targ = pyresample.geometry.SwathDefinition(lons=xz,lats=yz) # target grid
            uz = []
            for k in range(wu.shape[0]):
                uz.append(pyresample.kd_tree.resample_nearest(orig,wu[k,:,:],targ,radius_of_influence=50000))
         
            dic.update({'U1Z' :(['time','YZ','XZ'], uz)})
        
        if 'V1' in keys:        
            # Get water velocities - V1
            h =  ds['V1'][:,0,1:-1,1:-1] 
            ha = h.transpose(h.dims[0],h.dims[2],h.dims[1]) #transpose lat/lon
            ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
            wv = np.ma.masked_where(ww==True,ha) #mask land
            dic.update({'V1' :(['time','YV','XV'], wv)}) 
        
            #interpolate on XZ, YZ
            orig = pyresample.geometry.SwathDefinition(lons=xv,lats=yz) # original points
            targ = pyresample.geometry.SwathDefinition(lons=xz,lats=yz) # target grid
            vz = []
            for k in range(wv.shape[0]):
                vz.append(pyresample.kd_tree.resample_nearest(orig,wv[k,:,:],targ,radius_of_influence=50000))
      
            dic.update({'V1Z' : (['time','YZ','XZ'], vz)})
       
                 
        self.Dataset=xr.Dataset(dic,                                              
                     coords={'XZ':(xz[0,:]),'YZ':(yz[:,0]),
                             'XU':(xz[0,:]),'YU':(yu[:,0]),
                             'XV':(xv[0,:]),'YV':(yz[:,0]),                   
                     'time':self.times})
                     
                     
        #clean duplicates
        self.Dataset = self.Dataset.sel(time=~self.Dataset.indexes['time'].duplicated())
                     
        dic = self.info.copy()   # start with x's keys and values
        dic.update(kwargs)    # modifies z with y's keys and values & returns None
        
        if 'sa_date' not in dic.keys():
            dic.update({'sa_date':self.Dataset.time.values[0]})
            
        if 'se_date' not in dic.keys():
            dic.update({'se_date':self.Dataset.time.values[-1]})
               
        self.obs = obs(**dic)
                
    def hview(self,var,**kwargs):
        
        return hv.Dataset(self.Dataset[var])
   
   
    def gview(self,var,**kwargs):
        
        return gv.Dataset(self.Dataset[var])
    
                                   
    def frames(self,var,**kwargs):
        
        if len(var) == 1 :  
            return contour(self.xh,self.yh,self.Dataset[var[0]],self.times,**kwargs)
        elif len(var) == 2:
            return quiver(self.xh,self.yh,self.Dataset[var[0]],self.Dataset[var[1]],self.times,**kwargs)


class schism(data):
    
            
    def __init__(self,**kwargs):     
        
        rpath = kwargs.get('rpath','./')
        
        folders = kwargs.get('folders',None)  #[os.path.join(os.path.abspath(loc),name) for name in os.listdir(loc) if os.path.isdir(os.path.join(loc,name))]
        
        if folders :
            self.folders = folders
        else:
            self.folders = [rpath]

    
        ncores = get_value(self,kwargs,'ncores',1)
    
        sys.stdout.write('Combining output\n')
        sys.stdout.flush()  
        
        
        inu = len(glob.glob(self.folders[0] + 'outputs/schout_000*_1.nc'))
        
        datai=[]
        
        for folder in self.folders:
                        
            #netcdf files combined by autocombine_MPI_elfe.pl
            cfiles = glob.glob(folder+'outputs/schout_0.nc') + glob.glob(folder+'outputs/schout_[!00]*.nc')
            cfiles.sort()
            if cfiles :
                #read grid data
                gr = xr.open_mfdataset(folder + 'outputs/schout_0.nc', autoclose=True, drop_variables=['time'])
                #read the rest..
                var = xr.open_mfdataset(cfiles[1:], autoclose=True, drop_variables=gr.variables.keys())
  
                datai.append(xr.merge([gr,var])) #merge all
                continue
            
            #already saved complete netcdf
            cfile = glob.glob(folder+'outputs/schout.nc')
            if cfile :
                datai.append(xr.open_mfdataset(cfile, autoclose=True))
                continue
            
            #COMBINE ON THE FLY   
            gtol = glob.glob(folder+'outputs/global*') #get the global_to_local file

            gindx = pd.read_csv(gtol[0],header=None,delim_whitespace=True) # read the file
            
            gindx = gindx.set_index(gindx.columns[0]) # set the global index as index
            
            gindx.columns=['dist'] # rename the column to dist[ribution]
            
            
            gfiles = glob.glob(folder+'outputs/local*') # Read the global node index distribution to the cores
            gfiles.sort() # sort them
            
            #create a dict from filenames to identify parts in the dataframes below
            keys=[]
            for name in gfiles:
                keys.append('node{}'.format(name.split('/')[-1].split('_')[-1]))
            
            
            header = pd.read_csv(gfiles[0],header=None,nrows=1,delim_whitespace=True,names=['ns_global','ne_global','np_global','nvrt','nproc','ntracers','T','S','GEN','AGE','SED3D','EcoSim','ICM','CoSINE','Feco','TIMOR','FABM']) # read the header
            
            [nedges, nelems, npoints] = header.loc[:,['ns_global','ne_global','np_global']].values.flatten() # get relevant values
            
            # get the number of elems from all files
            nels = []
            for i in range(len(gfiles)):
                ne = pd.read_csv(gfiles[i],skiprows=2, header=None, nrows = 1)
                nels.append(ne.values.flatten()[0].astype(int))
            
            #read and add them to pandas DataFrame 
            frames=np.empty(len(gfiles),dtype=object)
            for i in range(len(gfiles)):
                frames[i] = pd.read_csv(gfiles[i],skiprows=3,header=None, nrows=nels[i], \
                            names=['local','global_n'], delim_whitespace=True)

            elems = pd.concat(frames,keys=keys)

            #get the number of nodes from all files
            nq = []
            for i in range(len(gfiles)):
                nn = pd.read_csv(gfiles[i],skiprows=nels[i] + 3, header=None, nrows = 1)
                nq.append(nn.values.flatten()[0].astype(int))
            
            #read and add them to pandas DataFrame
            nframes=np.empty(len(gfiles),dtype=object)
            for i in range(len(gfiles)):
                nframes[i] = pd.read_csv(gfiles[i],skiprows=nels[i] + 4,header=None, nrows=nq[i], \
                            names=['local','global_n'], delim_whitespace=True)
    
            nodes = pd.concat(nframes,keys=keys)

            #get the number of edges from all files
            nw = []
            for i in range(len(gfiles)):
                nb = pd.read_csv(gfiles[i],skiprows=nels[i] + nq[i] + 4, header=None, nrows = 1)
                nw.append(nb.values.flatten()[0].astype(int))
            
            #read and add them to pandas DataFrame
            wframes=np.empty(len(gfiles),dtype=object)
            for i in range(len(gfiles)):
                wframes[i] = pd.read_csv(gfiles[i],skiprows=nels[i] + nq[i] + 5,header=None, nrows=nw[i],
                            names=['local','global_n'], delim_whitespace=True)
    
            re = pd.concat(wframes,keys=keys)
                
            #read headers
            h0 = pd.read_csv(gfiles[0],skiprows=nels[0] + nq[0] + nw[0] + 6, header=None, nrows = 1, delim_whitespace=True,
                            names=['start_year','start_month','start_day','start_hour','utc_start'])
        
            h1 = pd.read_csv(gfiles[0],skiprows=nels[0] + nq[0] + nw[0] + 7, header=None, nrows = 1, delim_whitespace=True, 
                            names = ['nrec','dtout','nspool','nvrt','kz','h0','h_s','h_c','theta_b','theta_f','ics'])
            
            ztots = ['ztot_'+str(i) for i in range(1,h1.loc[:,'kz']-1)]
            
            sigmas = ['sigma_'+str(i) for i in range(h1.loc[:,'nvrt'] - h1.loc[:,'kz'] + 1) ]
            
            h2 = pd.read_csv(gfiles[0],skiprows=nels[0] + nq[0] + nw[0] + 8, header=None, nrows = 1, delim_whitespace=True, names=ztots + sigmas)
            
            
            #combine headers
            header2 = pd.concat([h0, h1, h2], axis=1)
            
            #read lat/lon from all files
            gframes=np.empty(len(gfiles),dtype=object)
            for i in range(len(gfiles)):
                gframes[i] = pd.read_csv(gfiles[i],skiprows=nels[i] + nq[i] + nw[i] + 10, header=None, nrows = nq[i], delim_whitespace=True, names=['lon','lat','depth','kbp00'])
    
            grid = pd.concat(gframes,keys=keys)

            #read tesselation from all files
            eframes=np.empty(len(gfiles),dtype=object)
            for i in range(len(gfiles)):
                eframes[i] = pd.read_csv(gfiles[i],skiprows=nels[i] + nq[i] + nw[i] + nq[i] + 10, header=None, nrows = nels[i], delim_whitespace=True, names=['type','a','b','c'])
    
            tri = pd.concat(eframes,keys=keys)

            # read output netcdf files
            files = glob.glob(folder+'outputs/schout_0*_*.nc')
            files.sort()

            #store them in a list
            out=[]
            for i in range(len(keys)):
                ifiles = [f for f in files if '{}_'.format(i) in f]
                out.append(xr.open_mfdataset(ifiles))
                
            #convert times to Timestamp
            date = header2.loc[:,['start_year','start_month','start_day','start_hour','utc_start']]
            date.columns=['year','month','day','hour','utc'] # rename the columns
            #set the start timestamp
            sdate = pd.Timestamp(year=date.year.values, month=date.month.values, day=date.day.values, hour=date.hour.values, tz=date.utc.values)
            #get times as timestamps
            times = pd.to_datetime(out[0].time.values, unit='s', origin=sdate)

            # combine grid
            cnodes = nodes.global_n.drop_duplicates() # drop duplicate global nodes
            
            dgrid = grid.drop_duplicates() # keep only one of the duplicates
            dgrid.index = dgrid.index.droplevel() # drop multi-index
            
            dgrid = dgrid.reset_index(drop=True) # reset index
            
            dgrid.index = cnodes.values - 1 # reindex based on the global index -1 for the python convention
            
            dgrid= dgrid.sort_index() #sort with the global index
            
            dgrid = dgrid.reset_index(drop=True) #reset index
            #dgrid = dgrid.dropna() #failcheck
            
            #Combine Tessellation
            for key in keys:
                nod = nodes.loc[key].copy() # make a copy of the first node
                nod = nod.set_index('local') # reset the index using the local values (in essense set the index to start from 1...)
                tri.loc[key,'ga'] = nod.reindex(tri.loc[key,'a'].values).values
                tri.loc[key,'gb'] = nod.reindex(tri.loc[key,'b'].values).values
                tri.loc[key,'gc'] = nod.reindex(tri.loc[key,'c'].values).values
                
            tri.loc[:,'ga'] = tri.loc[:,'ga'].apply(pd.to_numeric(int)) # make integer
            tri.loc[:,'gb'] = tri.loc[:,'gb'].apply(pd.to_numeric(int)) # make integer
            tri.loc[:,'gc'] = tri.loc[:,'gc'].apply(pd.to_numeric(int)) # make integer    
            
            #Sort elements
            gt3 = tri.loc[:,['ga','gb','gc']].copy() # make a copy
            gt3.index = gt3.index.droplevel() # drop multi-index
            gt3 = gt3.reset_index(drop=True)
            gt3.index = elems.global_n.values # we set the index equal to the global_n column
            gt3 = gt3.sort_index() #sort them
            #add nan column in place of the fourth node. This needs to be tested for quadrilaterals
            gt3['gd']=np.nan
            gt3 = gt3.reset_index() # reset to add more columns without problems
            ## Add mean x, y of the elememts. To be used in the output
            gt3['x1'] = dgrid.loc[gt3['ga'].values - 1, 'lon'].values #lon of the index, -1 for python convention
            gt3['y1'] = dgrid.loc[gt3['ga'].values - 1, 'lat'].values #lat of the index
            gt3['x2'] = dgrid.loc[gt3['gb'].values - 1, 'lon'].values
            gt3['y2'] = dgrid.loc[gt3['gb'].values - 1, 'lat'].values
            gt3['x3'] = dgrid.loc[gt3['gc'].values - 1, 'lon'].values
            gt3['y3'] = dgrid.loc[gt3['gc'].values - 1, 'lat'].values


            gt3['xc'] =  gt3[['x1', 'x2', 'x3']].mean(axis=1) #mean lon of the element
            gt3['yc'] =  gt3[['y1', 'y2', 'y3']].mean(axis=1)

            
            ## min kbe
            gt3['kbe1'] = dgrid.loc[gt3['ga'] - 1,'kbp00'].values
            gt3['kbe2'] = dgrid.loc[gt3['gb'] - 1,'kbp00'].values
            gt3['kbe3'] = dgrid.loc[gt3['gc'] - 1,'kbp00'].values

            gt3['kbe'] = gt3[['kbe1', 'kbe2', 'kbe3']].min(axis=1)
            
            gt3 = gt3.set_index('index') # set index back 
            gt34 = gt3.loc[:,['ga','gb','gc','gd']].values # SCHISM_hgrid_face_nodes
            
            #Grid's Edges
            edk=[]
            for key in keys:
                eds=[]
                for [ga,gb,gc] in tri.loc[key,['ga','gb','gc']].values:
                    eds.append([gb,gc])
                    eds.append([gc,ga])
                    eds.append([ga,gb])
        
                eds = np.array(eds)

                df = pd.DataFrame(eds)
                idsd = df[df.apply(sorted, axis=1).duplicated()].index.values # find the duplicates
                df_ = df.drop(idsd) #drop them 
                df_.index = re.loc[key,'global_n'].values
    
                edk.append(df_)#make a pandas DataFrame
            
            edgs = pd.concat(edk) # We concatenate, however there are dublicate indices ...
            
            #see https://stackoverflow.com/questions/13035764/remove-rows-with-duplicate-indices-pandas-dataframe-and-timeseries
            edgs1 = edgs.reset_index().drop_duplicates(subset='index', keep='first').set_index('index') #drop duplicates 
            edgs1 = edgs1.sort_index() #sort index 
            
            edgs1 = edgs1.reset_index() #reset index to add columns
            
            #mean x, y 
            edgs1['x1'] = dgrid.loc[edgs1[0].values - 1, 'lon'].values #lon of the index, -1 for python convention
            edgs1['y1'] = dgrid.loc[edgs1[0].values - 1, 'lat'].values #lat of the index
            edgs1['x2'] = dgrid.loc[edgs1[1].values - 1, 'lon'].values
            edgs1['y2'] = dgrid.loc[edgs1[1].values - 1, 'lat'].values
 
            edgs1['xc'] =  edgs1[['x1', 'x2']].mean(axis=1) #mean of the edge index
            edgs1['yc'] =  edgs1[['y1', 'y2']].mean(axis=1)

            ## min bottom index
            edgs1['kbs1'] = dgrid.loc[edgs1[0] - 1,'kbp00'].values
            edgs1['kbs2'] = dgrid.loc[edgs1[1] - 1,'kbp00'].values

            edgs1['kbs'] = edgs1[['kbs1', 'kbs2']].min(axis=1)
            
            edgs1.set_index('index') # set index again 
                                         
            #Combine element based variables
            s = pd.Series(gindx.values.flatten()) # create a series with the elements node reference
            
            dat = pd.concat([s] * times.shape[0], axis=1) # concatenate to the number of time steps
            
            dat.columns=times # set columns names as the timestamps
            
            edic={}
            for var in out[0].variables.keys():
             #   print out[0][var].name, out[0][var].dims, len(out[0][var].dims)
                if ('nSCHISM_hgrid_face' in out[0][var].dims) & (len(out[0][var].dims) == 2):
                    wd = dat.copy()
                    for time in times: #all times
                        for i in range(len(keys)): # all components
                            wd.loc[dat[time] == i] = out[i][var].values.T #
                    vname = out[0][var].name
                    edic[vname]=wd.T.values
        
                elif ('nSCHISM_hgrid_face' in out[0][var].dims) & ('two' in out[0][var].dims):
                    wd = dat.copy()
                    for time in times: #all times
                        for i in range(len(keys)): # all components
                            wd.loc[dat[time] == i] = out[i][var].values[:,:,0].T # wetdry_elem variable
                    vx = wd.T
                
                    wd = dat.copy()
                    for time in times: #all times
                        for i in range(len(keys)): # all components
                            wd.loc[dat[time] == i] = out[i][var].values[:,:,1].T # wetdry_elem variable

                    vy = wd.T

                    vname = out[0][var].name
                    edic[vname] = np.dstack([vx.values,vy.values])
        
                elif ('nSCHISM_hgrid_face' in out[0][var].dims) & ('nSCHISM_vgrid_layers' in out[0][var].dims):
                    s=out[0][var].shape[2]
                    ars=[]
                    for l in range(s):
                        wd = dat.copy()
                        for time in times: #all times
                            for i in range(len(keys)): # all components
                                wd.loc[dat[time] == i] = out[i][var].values[:,:,0].T # wetdry_elem variable

                        ars.append(wd.T)

                    vname = out[0][var].name
                    edic[vname] = np.dstack([v.values for v in ars])


            

            #Combine node based variables
            drop_mask = grid.duplicated().values # save the mask
            
            vdic={}
            for var in out[0].variables.keys():
                #print out[0][var].name, out[0][var].dims, len(out[0][var].dims)
                if ('nSCHISM_hgrid_node' in out[0][var].dims) & (len(out[0][var].dims) == 2):
                    ars = [v[out[0][var].name] for v in out]
       
                    v = combine(ars, drop_mask, cnodes, times)
        
                    vdic[out[0][var].name] = v.values
        
                elif ('nSCHISM_hgrid_node' in out[0][var].dims) & ('two' in out[0][var].dims) & (len(out[0][var].dims) == 3):
                    ars1 = [v[out[0][var].name][:,:,0] for v in out]
                    vx = combine(ars1, drop_mask, cnodes, times)
                    ars2 = [v[out[0][var].name][:,:,1] for v in out]
                    vy = combine(ars2, drop_mask, cnodes, times)
        
                    vname = out[0][var].name
                    vdic[vname] = np.dstack([vx.values,vy.values])
        
                elif ('nSCHISM_hgrid_node' in out[0][var].dims) & ('nSCHISM_vgrid_layers' in out[0][var].dims) & (len(out[0][var].dims) == 3):
                    s=out[0][var].shape[2]
                    ars=[]
                    for l in range(s):
                        arsi = [v[out[0][var].name][:,:,l] for v in out]
                        ars.append(combine(arsi, drop_mask, cnodes, times))

                    vname = out[0][var].name
                    vdic[vname] = np.dstack([v.values for v in ars])

                elif ('nSCHISM_hgrid_node' in out[0][var].dims) & ('nSCHISM_vgrid_layers' in out[0][var].dims) & ('two' in out[0][var].dims) & (len(out[0][var].dims) == 4):
                    s=out[0][var].shape[2]
                    ars=[]
                    for l in range(s):
                        arsx = [v[out[0][var].name][:,:,l,0] for v in out]
                        vx = combine(arsx, drop_mask, cnodes, times)
                        arsy = [v[out[0][var].name][:,:,l,1] for v in out]
                        vy = combine(arsy, drop_mask, cnodes, times)
                        ars.append(np.dstack([vx.values,vy.values]))
        
                    vname = out[0][var].name
                    vdic[vname] = np.stack([a for a in ars], axis=2) # stack correctly
            
            #Create output structure
            sigms = header2.loc[:,sigmas].values.flatten() # get sigmas
            iwet_dry = kwargs.get('iwet_dry',0)  # defined by the user
            ihgrid_id = kwargs.get('ihgrid_id',-2147483647) # defined by user - 0,dummy_dim,ihgrid_id
            
            #element based
            xrdic={}
            for key in edic.iterkeys():
                xrdic.update({key:([x for x in out[0][key].dims],edic[key])})
            
            xe = xr.Dataset(xrdic,coords={u'time':times, u'sigma': sigms })
            
                
            #node based
            xrdic={}
            for key in vdic.iterkeys():
                xrdic.update({key:([x for x in out[0][key].dims],vdic[key])})
                
            xn = xr.Dataset(xrdic,coords={u'time':times, u'sigma': sigms })
            
            #grid based
            
            #compute cs
            klev = np.arange(header2.kz,header2.nvrt+1)
            k = klev-header2.kz.values


            cs=np.zeros(k)

            cs=(1-header2.theta_b.values)*np.sinh(header2.theta_f.values*sigms[k])/np.sinh(header2.theta_f.values)+ \
                header2.theta_b.values*(np.tanh(header2.theta_f.values*(sigms[k]+0.5))-np.tanh(header2.theta_f.values*0.5))/2/np.tanh(header2.theta_f.values*0.5)

            
            xg = xr.Dataset({u'SCHISM_hgrid' : ([u'one'], [ihgrid_id]),
                             u'SCHISM_hgrid_face_nodes' : ([u'nSCHISM_hgrid_face', u'nMaxSCHISM_hgrid_face_nodes'], gt34),
                             u'SCHISM_hgrid_edge_nodes' : ([u'nSCHISM_hgrid_edge', u'two'], edgs1[[0,1]].values),
                             u'SCHISM_hgrid_node_x' : ([u'nSCHISM_hgrid_node'], dgrid.lon.values) ,
                             u'SCHISM_hgrid_node_y' : ([u'nSCHISM_hgrid_node'], dgrid.lat.values) ,
                             u'SCHISM_hgrid_edge_x' : ([u'nSCHISM_hgrid_edge'], edgs1['xc'].values),
                             u'SCHISM_hgrid_edge_y' : ([u'nSCHISM_hgrid_edge'], edgs1['yc'].values ),
                             u'SCHISM_hgrid_face_x' : ([u'nSCHISM_hgrid_face'], gt3.loc[:,'xc'].values),
                             u'SCHISM_hgrid_face_y' : ([u'nSCHISM_hgrid_face'], gt3.loc[:,'yc'].values),
                             u'depth': ([u'nSCHISM_hgrid_node'], dgrid.depth.values) ,
                             u'minimum_depth': ([u'one'], header2.loc[:,'h0'].values),
                             u'coordinate_system_flag' : ([u'one'], header2.loc[:,'ics'].values),
                             u'sigma_h_c' : ([u'one'], header2.loc[:,'h_c'].values),
                             u'sigma_theta_b': ([u'one'], header2.loc[:,'theta_b'].values),
                             u'sigma_theta_f' : ([u'one'], header2.loc[:,'theta_f'].values),
                             u'sigma_maxdepth' : ([u'one'], header2.loc[:,'h_s'].values),
                             u'Cs' : (['sigma'], cs),
                             u'dry_value_flag' : ([u'one'], [iwet_dry] ),
                             u'ele_bottom_index': ([u'nSCHISM_hgrid_face'], gt3.kbe.values ),
                             u'node_bottom_index' : ([u'nSCHISM_hgrid_node'], dgrid.kbp00.values),
                             u'edge_bottom_index' : ([u'nSCHISM_hgrid_edge'], edgs1.kbs.values),
                             },
                                 coords={u'time':times, u'sigma': sigms })
                                 
                                 
            #Choose attrs
            if header2.ics.values == 1:
                lat_coord_standard_name = 'projection_x_coordinate'
                lon_coord_standard_name = 'projection_y_coordinate'
                x_units = 'm'
                y_units = 'm'
                lat_str_len = 23
                lon_str_len = 23
            else:
                lat_coord_standard_name = 'latitude'
                lon_coord_standard_name = 'longitude'
                x_units = 'degrees_north'
                y_units = 'degrees_east'
                lat_str_len = 8
                lon_str_len = 9
                
                
            #Attrs
            xg.SCHISM_hgrid_node_x.attrs = {'long_name' : 'node x-coordinate', 'standard_name' : lon_coord_standard_name , 'units' : x_units, 'mesh' : 'SCHISM_hgrid'}

            xg.SCHISM_hgrid_node_y.attrs = {'long_name' : 'node y-coordinate', 'standard_name' : lat_coord_standard_name , 'units' : y_units, 'mesh' : 'SCHISM_hgrid'}

            xg.depth.attrs = {'long_name' : 'Bathymetry', 'units' : 'meters', 'positive' : 'down', 'mesh' : 'SCHISM_hgrid', 'location' : 'node'}

            xg.sigma_h_c.attrs = {'long_name' : 'ocean_s_coordinate h_c constant', 'units' : 'meters', 'positive' : 'down'}

            xg.sigma_theta_b.attrs = {'long_name' : 'ocean_s_coordinate theta_b constant'}

            xg.sigma_theta_f.attrs = {'long_name' : 'ocean_s_coordinate theta_f constant'}

            xg.sigma_maxdepth.attrs = {'long_name' : 'ocean_s_coordinate maximum depth cutoff (mixed s over z bound...', 'units' : 'meters', 'positive' : 'down'}

            xg.Cs.attrs = {'long_name' : 'Function C(s) at whole levels', 'positive' : 'up' }

            xg.dry_value_flag.attrs = {'values' : '0: use last-wet value; 1: use junk'}

            xg.SCHISM_hgrid_face_nodes.attrs = {'long_name' : 'Horizontal Element Table', 'cf_role' : 'face_node_connectivity' , 'start_index' : 1}

            xg.SCHISM_hgrid_edge_nodes.attrs = {'long_name' : 'Map every edge to the two nodes that it connects', 'cf_role' : 'edge_node_connectivity' , 'start_index' : 1}

            xg.SCHISM_hgrid_edge_x.attrs = {'long_name' : 'x_coordinate of 2D mesh edge' , 'standard_name' : lon_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xg.SCHISM_hgrid_edge_y.attrs = {'long_name' : 'y_coordinate of 2D mesh edge' , 'standard_name' : lat_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xg.SCHISM_hgrid_face_x.attrs = {'long_name' : 'x_coordinate of 2D mesh face' , 'standard_name' : lon_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xg.SCHISM_hgrid_face_y.attrs = {'long_name' : 'y_coordinate of 2D mesh face' , 'standard_name' : lat_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xg.node_bottom_index.attrs = {'long_name' : 'bottom level index at each node' , 'units' : 'non-dimensional', 'mesh' : 'SCHISM_hgrid', 'location' : 'node',
                'start_index' : 1}
            xg.ele_bottom_index.attrs = {'long_name' : 'bottom level index at each element' , 'units' : 'non-dimensional', 'mesh' : 'SCHISM_hgrid', 'location' : 'elem',
                'start_index' : 1}
            xg.edge_bottom_index.attrs = {'long_name' : 'bottom level index at each edge' , 'units' : 'non-dimensional', 'mesh' : 'SCHISM_hgrid', 'location' : 'edge',
                'start_index' : 1}

            xg.SCHISM_hgrid.attrs = {'long_name' : 'Topology data of 2d unstructured mesh',
                                       'topology_dimension' : 2,
                                       'cf_role' : 'mesh_topology',
                                       'node_coordinates' : 'SCHISM_hgrid_node_x SCHISM_hgrid_node_y',
                                       'face_node_connectivity' : 'SCHISM_hgrid_face_nodes',
                                       'edge_coordinates' : 'SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y',
                                       'face_coordinates' : 'SCHISM_hgrid_face_x SCHISM_hgrid_face_y',
                                       'edge_node_connectivity' : 'SCHISM_hgrid_edge_nodes'
                                      }
            

            dat = xr.merge([xe,xn,xg])#MERGE
             
            dat.attrs = {'Conventions': 'CF-1.0, UGRID-1.0', 'title': 'SCHISM Model output', 'source': 'SCHISM model output version v10', 'references': 'http://ccrm.vims.edu/schismweb/',
                         'history': 'created by pyPoseidon', 'comment': 'SCHISM Model output', 'type': 'SCHISM Model output', 'VisIT_plugin': 'https://schism.water.ca.gov/library/-/document_library/view/3476283' }
            
            
            datai.append(dat) #append to list 
            
            savenc = kwargs.get('savenc',False)
            
            if savenc :
                dat.to_netcdf(folder+'outputs/schout.nc')                             
    

        self.Dataset = xr.merge(datai) #save final xarray
                
        tag = kwargs.get('tag', None)
        
        dic={}
        
        try:
            
            with open(self.folders[0]+tag +'_info.pkl', 'r') as f:
                      self.info=pickle.load(f)  
    
            dic = self.info.copy()   # start with x's keys and values
            dic.update(kwargs)    # modifies z with y's keys and values & returns None
    
        except:
            pass
    
        if 'sa_date' not in dic.keys():
            dic.update({'sa_date':self.Dataset.time.values[0]})
        
        if 'se_date' not in dic.keys():
            dic.update({'se_date':self.Dataset.time.values[-1]})
            
        if 'minlon' not in dic.keys():
            dic.update({'minlon':self.Dataset.SCHISM_hgrid_node_x.values.min()})

        if 'maxlon' not in dic.keys():
            dic.update({'maxlon':self.Dataset.SCHISM_hgrid_node_x.values.max()})
            
        if 'minlat' not in dic.keys():
            dic.update({'minlat':self.Dataset.SCHISM_hgrid_node_y.values.min()})
        
        if 'maxlat' not in dic.keys():
            dic.update({'maxlat':self.Dataset.SCHISM_hgrid_node_y.values.max()})
        
        self.obs = obs(**dic)        
        
    
    #Function for combining variables
def combine(ars, drop_mask, cnodes, times):
            pout = ars[0].to_pandas().T
            for f in ars[1:]:
                pout = pd.concat([pout, f.to_pandas().T])
        
            pout = pout.reset_index(drop=True) # reset index
        
            pout = pout.drop(pout[drop_mask].index) # drop duplicate nodes
        
            pout.index = cnodes.values - 1 # reindex based on the global index -1 for the python convention
        
            pout = pout.sort_index() #sort with the global index
        
        
            pout = pout.reset_index(drop=True)#reindex for final version
        
            pout.columns = times # set time stamp 
        
            return pout.T # transpose to set time as index         
    
    
                           
class point:
    
    def __init__(self,**kwargs):
                        
        self.lon = kwargs.get('lon', None) 
        self.lat = kwargs.get('lat', None) 
        self.data = kwargs.get('data', None)
            
    def tseries(self,**kwargs):
        
        var = kwargs.get('var', None)
        method = kwargs.get('method', 'nearest')
        
        plat=float(self.lat)
        plon=float(self.lon)
 
        i=np.abs(self.data.xh[0,:].data-plon).argmin()
        j=np.abs(self.data.yh[:,0].data-plat).argmin()
                       
        xb, yb = self.data.xh[j-5:j+5,i-5:i+5],self.data.yh[j-5:j+5,i-5:i+5]
                       
        vals = self.data.Dataset[var][:,j-5:j+5,i-5:i+5].values
                        
        orig = pyresample.geometry.SwathDefinition(lons=xb,lats=yb) # create original swath grid
                                   
        targ = pyresample.geometry.SwathDefinition(lons=np.array([plon,plon]),lats=np.array([plat,plat])) #  point
        
        svals = []
                
        if method == 'nearest':
            for k in range(vals.shape[0]):
                s = pyresample.kd_tree.resample_nearest(orig,vals[k,:,:],targ,radius_of_influence=50000,fill_value=999)
                svals.append(s[0])
        elif method == 'gauss':
            for k in range(vals.shape[0]):
                s = pyresample.kd_tree.resample_gauss(orig,vals[k,:,:],targ,radius_of_influence=50000,fill_value=999,sigmas=25000)                
                svals.append(s[0])
                

        pdata = pd.DataFrame({'time':self.data.Dataset[var].time, self.data.Dataset[var].name : svals})
        setattr(self, self.data.Dataset[var].name, pdata.set_index(['time']))
        
        
    def stseries(self,**kwargs):
        
        var = kwargs.get('var', None)
        method = kwargs.get('method', 'nearest')
        
        plat=float(self.lat)
        plon=float(self.lon)
        
        points = pd.concat([self.data.SCHISM_hgrid_node_x[:].to_dataframe(), self.data.SCHISM_hgrid_node_y[:].to_dataframe()], axis=1)
        values = self.data[var].values
        
        svals = []
        for k in range(self.data[var].time.size):
            svals.append( scipy.interpolate.griddata(points.values, values[k,:], (plat, plon), method=method) )   
        
        pdata = pd.DataFrame({'time':self.data[var].time, self.data[var].name : svals})
        setattr(self, self.data[var].name, pdata.set_index(['time']))
        
        
