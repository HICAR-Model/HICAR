import xesmf as xe
import numpy as np
import pyproj
#import ESMF
import xarray as xr

#ESMF.Manager(debug=True)
#######################################################
##  Following is code from William "Ryan" Currier   ##
########  email:  william.r.currier@noaa.gov   #######
######################################################

def get_latlon_b(ds,lon_str,lat_str,lon_dim,lat_dim):
    
    lon_raw = ds[lon_str][:,:]
    lat_raw = ds[lat_str][:,:]
   # if (len(ds[lon_str].shape) == 3):
   #     lon_raw = ds[lon_str][0,:,:]
   #     lat_raw = ds[lat_str][0,:,:]
        
    #### Longitude East-West Stagger
    diffWestEast    = lon_raw.diff(lon_dim)                         # take difference in longitudes across west-east (columns) dimension
    diffWestEast    = np.array(diffWestEast)                            # assign this to numpy array
    padding         = diffWestEast[:,-1].reshape(len(diffWestEast[:,-1]),1)  # get the last column of the differences and reshape so it can be appended
    diffWestEastAll = np.append(diffWestEast, padding, axis=1)/2               # append last value to all_data
    # dimensions now the same as ds
    lon_b_orig      = lon_raw-diffWestEastAll
    # lon_b needs to be staggered - add to the final row to bound
    last_add     = lon_raw[:,-1]+diffWestEastAll[:,-1]
    last_add     = np.array(last_add)
    last_add     = last_add.reshape(len(last_add),1)
    lon_b_append = np.append(np.array(lon_b_orig),last_add,1)
    last_add     = lon_b_append[0,:].reshape(1,len(lon_b_append[0,:]))
    lon_b_append = np.append(last_add,lon_b_append,axis=0)
    
    
    #### Latitude Stagger
    diffSouthNorth=-lat_raw.diff(lat_dim)                         # take difference in longitudes across west-east (columns) dimension
    diffSouthNorth=np.array(diffSouthNorth)                            # assign this to numpy array
    padding=diffSouthNorth[0,:].reshape(1,len(diffSouthNorth[0,:]))  # get the last column of the differences and reshape so it can be appended
    diffSouthNorthAll = np.append(padding,diffSouthNorth,axis=0)/2               # append last value to all_data
    # dimensions now the same as ds
    lat_b_orig      = lat_raw-diffSouthNorthAll
    # lat_b needs to be staggered - add to the first row to bound
    last_add     = lat_raw[0,:]+diffSouthNorthAll[0,:]
    last_add     = np.array(last_add)
    last_add     = last_add.reshape(1,len(last_add))
    lat_b_append = np.append(last_add,np.array(lat_b_orig),axis=0)
    last_add     = lat_b_append[:,-1].reshape(len(lat_b_append[:,-1]),1)
    lat_b_append = np.append(lat_b_append,last_add,axis=1)
    
   # if (len(ds[lon_str].shape) == 3):
   #     lon_b_append = np.tile(lon_b_append,ds[lon_str].shape[0])
   #     lat_b_append = np.tile(lat_b_append,ds[lat_str].shape[0])
        
    grid_with_bounds = {'lon': ds[lon_str].values,
                               'lat': ds[lat_str].values,
                               'lon_b': lon_b_append,
                               'lat_b': lat_b_append,
                              }
    return grid_with_bounds

def get_latlon_b_rect(ds,lon_str,lat_str,lon_dim,lat_dim):
    #### Longitude Stagger
    diffWestEast    = ds[lon_str].diff(lon_dim)                         # take difference in longitudes across west-east (columns) dimension
    diffWestEastArr = np.array(diffWestEast).reshape(len(diffWestEast),1)                            # assign this to numpy array
    padding         = diffWestEastArr[-1].reshape(1,1)  # get the last column of the differences and reshape so it can be appended
    diffWestEastAll = np.append(diffWestEastArr, padding, axis=0)/2               # append last value to all_data
    # # dimensions now the same as ds
    lon_b_orig      = ds[lon_str].values.reshape(len(ds[lon_str]),1)-diffWestEastAll
    # lon_b needs to be staggered - add to the final row to bound
    last_add        = ds[lon_str][-1].values+diffWestEastAll[-1]
    last_add        = np.array(last_add).reshape(len(last_add),1)
    last_add        = last_add.reshape(1,1)
    lon_b_append    = np.append(np.array(lon_b_orig),last_add,0)
    #### Latitude Stagger
    diffSouthNorth    = ds[lat_str].diff(lat_dim)                         # take difference in latitudes across west-east (columns) dimension
    diffSouthNorthArr = np.array(diffSouthNorth).reshape(len(diffSouthNorth),1)                            # assign this to numpy array
    padding         = diffSouthNorthArr[-1].reshape(1,1)  # get the last column of the differences and reshape so it can be appended
    diffSouthNorthAll = np.append(diffSouthNorthArr, padding, axis=0)/2               # append last value to all_data
    # # dimensions now the same as ds
    lat_b_orig      = ds[lat_str].values.reshape(len(ds[lat_str]),1)-diffSouthNorthAll
    # lat_b needs to be staggered - add to the final row to bound
    last_add        = ds[lat_str][-1].values+diffSouthNorthAll[-1]
    last_add        = np.array(last_add).reshape(len(last_add),1)
    last_add        = last_add.reshape(1,1)
    lat_b_append    = np.append(np.array(lat_b_orig),last_add,0)
    grid_with_bounds = {'lon': ds[lon_str],
                               'lat': ds[lat_str],
                               'lon_b': lon_b_append.reshape(len(lon_b_append),),
                               'lat_b': lat_b_append.reshape(len(lat_b_append),),
                              }
    return grid_with_bounds


#Usage examples
# wrf_grid_with_bounds = get_latlon_b(dsWRF['PREC_ACC_NC'],lon_str='XLONG',lat_str='XLAT',lon_dim='west_east',lat_dim='south_north')
# # wrf_grid_with_bounds
# ICAR_grid_with_bounds = get_latlon_b(ds,lon_str='lon',lat_str='lat',lon_dim='lon_x',lat_dim='lat_y')

# regridder = xe.Regridder(wrf_grid_with_bounds, ICAR_grid_with_bounds, 'conservative',reuse_weights=True) # input grid, grid you want to resample to, method
# dsWRF_out = regridder(dsWRF)



##############################################
##  Following is code from Dylan Reynolds   ##
######  email:  dylan.reynolds@slf.ch   ######
##############################################

def prepDsForESMF(ds,lon_str,lat_str,lon_dim,lat_dim):

    lon_raw = ds[lon_str][:,:]
    lat_raw = ds[lat_str][:,:]

    lon_b = np.zeros((lon_raw.shape[0],lon_raw.shape[1],4))
    lat_b = np.zeros(lon_b.shape)
    
    #### Longitude East-West Stagger
    diffWestEast    = lon_raw.diff(lon_dim)        # take difference in longitudes across west-east (columns) dimension
    diffWestEast    = np.array(diffWestEast)                            # assign this to numpy array
    padding         = diffWestEast[:,-1].reshape(len(diffWestEast[:,-1]),1)  # get the last column of the differences and reshape so it can be appended
    diffWestEastAll = np.append(diffWestEast, padding, axis=1)/2               # append last value to all_data
    
    diffSouthNorth=-lon_raw.diff(lat_dim)          # take difference in longitudes across west-east (columns) dimension
    diffSouthNorth=np.array(diffSouthNorth)                            # assign this to numpy array
    padding=diffSouthNorth[0,:].reshape(1,len(diffSouthNorth[0,:]))  # get the last column of the differences and reshape so it can be appended
    diffSouthNorthAll = np.append(padding,diffSouthNorth,axis=0)/2               # append last value to all_data
        
    # dimensions now the same as ds
    lon_b[:,:,0]      = lon_raw+diffWestEastAll+diffSouthNorthAll
    lon_b[:,:,1]      = lon_raw+diffWestEastAll-diffSouthNorthAll
    lon_b[:,:,2]      = lon_raw-diffWestEastAll-diffSouthNorthAll
    lon_b[:,:,3]      = lon_raw-diffWestEastAll+diffSouthNorthAll

    
    #### Latitude Stagger
    diffSouthNorth=-lat_raw.diff(lat_dim)          # take difference in longitudes across west-east (columns) dimension
    diffSouthNorth=np.array(diffSouthNorth)                            # assign this to numpy array
    padding=diffSouthNorth[0,:].reshape(1,len(diffSouthNorth[0,:]))  # get the last column of the differences and reshape so it can be appended
    diffSouthNorthAll = np.append(padding,diffSouthNorth,axis=0)/2               # append last value to all_data
    
    diffWestEast    = lat_raw.diff(lon_dim)        # take difference in longitudes across west-east (columns) dimension
    diffWestEast    = np.array(diffWestEast)                            # assign this to numpy array
    padding         = diffWestEast[:,-1].reshape(len(diffWestEast[:,-1]),1)  # get the last column of the differences and reshape so it can be appended
    diffWestEastAll = np.append(diffWestEast, padding, axis=1)/2               # append last value to all_data
    
    # dimensions now the same as ds
    lat_b[:,:,0]      = lat_raw+diffWestEastAll+diffSouthNorthAll
    lat_b[:,:,1]      = lat_raw+diffWestEastAll-diffSouthNorthAll
    lat_b[:,:,2]      = lat_raw-diffWestEastAll-diffSouthNorthAll
    lat_b[:,:,3]      = lat_raw-diffWestEastAll+diffSouthNorthAll

    ds_template = ds.copy()
    ds_template['lon_b'] = xr.DataArray(lon_b, dims=[lat_dim, lon_dim, 'bounds'])
    ds_template['lat_b'] = xr.DataArray(lat_b, dims=[lat_dim, lon_dim, 'bounds'])
    ds_template['lat'].attrs['bounds'] = "lat_b"
    ds_template['lon'].attrs['bounds'] = "lon_b"
    
    if ("mask" in ds_template.data_vars): ds_template['mask'].attrs['missing_value'] = 0

    return ds_template

def coordsToLatLon(ds,x_var,y_var,x_dim,y_dim,in_proj,out_proj,inplace=False):
    
    #Do pyproj reprojection from given crs to 4326 (or, other) for xesmf reproj
    
    #if not(inplace): ds = ds.copy()
    
    #Setup reprojection transformer
    #in_proj = pyproj.CRS(('EPSG:' + crs))
    #out_proj = pyproj.CRS(('EPSG:' + out_crs))
    tg = pyproj.transformer.TransformerGroup(in_proj, out_proj)
    for trans in tg.transformers:
        print(repr(trans))
    
    tformer = pyproj.Transformer.from_crs(in_proj, out_proj,always_xy=True)

    print(tformer)
    
    #If we have 1D, rectilinear data, convert to 2D grid (even if redundant)
    if (len(ds[x_var].dims) == 1):
        x_len = len(ds[x_var].values)
        y_len = len(ds[y_var].values)
        x_grid, y_grid = np.meshgrid(ds[x_var].values,ds[y_var].values)
        
        #reprojected x y coords
        xs, ys = tformer.transform(x_grid.flat, y_grid.flat) 
        
    else:
        x_len = ds[x_var].shape[1]
        y_len = ds[y_var].shape[0]
        xs, ys = tformer.transform(ds[x_var].values.flat,ds[y_var].values.flat) 

    #out output will be as 1D arrays, so reshape back to the input x-y shape
    xs = np.reshape(xs,(y_len,x_len))
    ys = np.reshape(ys,(y_len,x_len))

    #stuff coordinates into 'lat' and 'lon' variables so that xesmf can do a regridding
    old_attrs = ds[x_var].attrs
    ds['lon'] = ((y_dim,x_dim),xs)
    ds['lon'].attrs = old_attrs

    old_attrs = ds[y_var].attrs
    ds['lat'] = ((y_dim,x_dim),ys)
    ds['lat'].attrs = old_attrs
    
    if (inplace): return
    return ds


def prepGridsForRegridder(ds_in,to_ds,ds_in_crs=None,ds_in_xvar='lon',ds_in_yvar='lat',ds_in_xdim='x',ds_in_ydim='y',\
                          to_ds_xvar='lon',to_ds_yvar='lat',to_ds_xdim='lon_x',to_ds_ydim='lat_y'):
    
    WGS64_crs = pyproj.CRS('EPSG:4326')
    if not(ds_in_crs): ds_in_crs = WGS64_crs
    
    #Setup to_ds grid object assuming standard ICAR output (for our purposes)
    to_ds_grid_with_bounds = get_latlon_b(to_ds,lon_str=to_ds_xvar,lat_str=to_ds_yvar,lon_dim=to_ds_xdim,lat_dim=to_ds_ydim)
    
    #If ds_in_crs is not 4326 (what we expect to have), create lat/lon coordinates in ds
    #Also, if we have rectilinear input, coordsToLatLon will put it on curvilinear grid
    if (not(ds_in_crs==WGS64_crs) or (len(ds_in[ds_in_xvar].dims) == 1)): 
        print('projecting coords to WGS')
        ds_in_b = coordsToLatLon(ds_in,ds_in_xvar,ds_in_yvar,ds_in_xdim,ds_in_ydim,in_proj=ds_in_crs,out_proj=WGS64_crs)
        ds_in_grid_with_bounds = get_latlon_b(ds_in_b,lon_str='lon',lat_str='lat',lon_dim=ds_in_xdim,lat_dim=ds_in_ydim)
    
    #If we have lat/lon, and it is curvilinear:
    else:
        ds_in_grid_with_bounds = \
                get_latlon_b(ds_in,lon_str=ds_in_xvar,lat_str=ds_in_yvar,lon_dim=ds_in_xdim,lat_dim=ds_in_ydim)
        
    if ("mask" in ds_in.data_vars): ds_in_grid_with_bounds["mask"] = ds_in["mask"]
    if ("mask" in to_ds.data_vars): to_ds_grid_with_bounds["mask"] = to_ds["mask"]
        
    return ds_in_grid_with_bounds, to_ds_grid_with_bounds

def regridToDs(ds_in,to_ds,ds_in_crs=None,ds_in_xvar='lon',ds_in_yvar='lat',ds_in_xdim='x',ds_in_ydim='y',to_ds_xvar='lon',\
               to_ds_yvar='lat',to_ds_xdim='lon_x',to_ds_ydim='lat_y',weights_fn=None,method='conservative_normed',extrap=None):
    #Function takes in a dataset (ds_in) and regrids/reprojects to the target ds, which is assumed to be in EPSG 4326
    
    # ds_in -- the dataset that will be regridded
    # to_ds -- the dataset(or simply grid) which ds_in will be regridded to match
    # ds_in_crs -- the CRS of ds_in, provided as a pyproj CRS object
    # ds_in_xvar -- the x-coordinate variable name of ds_in
    # ds_in_yvar -- the y-coordinate variable name of ds_in
    # ds_in_xdim -- the x-dimension name of ds_in
    # ds_in_ydim -- the y-dimension name of ds_in
    # to_ds_xvar -- the x-coordinate variable name of to_ds
    # to_ds_yvar -- the y-coordinate variable name of to_ds
    # to_ds_xdim -- the x-dimension name of to_ds
    # to_ds_ydim -- the y-dimension name of to_ds
    
    #If ds_in has lat-lon on 3D grid, convert to 2D for regridding
    if (len(ds_in[ds_in_xvar].shape) == 3):
        temp_lat = ds_in[ds_in_xvar][0,:,:]
        temp_lon = ds_in[ds_in_yvar][0,:,:]
        ds_in.drop([ds_in_xvar,ds_in_yvar])
        ds_in[ds_in_xvar] = temp_lat
        ds_in[ds_in_yvar] = temp_lon
        
    ds_in_grid_with_bounds, to_ds_grid_with_bounds = \
prepGridsForRegridder(ds_in,to_ds,ds_in_crs,ds_in_xvar,ds_in_yvar,ds_in_xdim,\
                      ds_in_ydim,to_ds_xvar,to_ds_yvar,to_ds_xdim,to_ds_ydim)
                
    reuse = False
    if (weights_fn): reuse=True
    # input grid, grid you want to resample to, method
    if reuse:
        regridder = xe.Regridder(ds_in_grid_with_bounds, to_ds_grid_with_bounds, \
                             method=method,extrap_method=extrap,filename=weights_fn,reuse_weights=reuse)
        disc_regridder = xe.Regridder(ds_in_grid_with_bounds, to_ds_grid_with_bounds, \
                             method=method,extrap_method=extrap,filename=weights_fn,reuse_weights=reuse)

    else:
        regridder = xe.Regridder(ds_in_grid_with_bounds, to_ds_grid_with_bounds, \
                             method=method,extrap_method=extrap,reuse_weights=reuse)
        disc_regridder = xe.Regridder(ds_in_grid_with_bounds, to_ds_grid_with_bounds, \
                             method='nearest_s2d',extrap_method=extrap,reuse_weights=reuse)
        
    ds_out = regridder(ds_in)
    ds_out_disc = disc_regridder(ds_in)
    
    ds_out.attrs = ds_in.attrs
    
    possible_discrete = ['landuse','landmask','ridge_mask','valley_mask']
    
    for field in possible_discrete:
        if field in ds_out.keys():
            ds_out[field].values = ds_out_disc[field].values
    
    return ds_out, regridder

def unstaggerICAR(ICAR):
        
    var_arr = []
        
    for v in ICAR.data_vars:
         if (v.lower()=='u' or v.lower()=='v'):
            print('unstaggering: '+v)
            temp_attrs = ICAR[v].attrs
            new_dims = list(ICAR[v].dims)
            new_dims[-1] = 'lon_x'
            new_dims[-2] = 'lat_y'
            new_dims = tuple(new_dims)
            
            temp_vals = ICAR[v]
            
            if (len(ICAR[v].dims) == 3):
                #Gather values, handling staggared U/V
                if (v.lower()=='u'):
                    temp_vals = (temp_vals[:,:,1:] + temp_vals[:,:,:-1])/2
                    temp_vals = temp_vals.swap_dims({'lon_u':'lon_x'})
                elif (v.lower()=='v'):
                    temp_vals = (temp_vals[:,1:,:] + temp_vals[:,:-1,:])/2
                    temp_vals = temp_vals.swap_dims({'lat_v':'lat_y'})
            elif (len(ICAR[v].dims) == 4):
                #Gather values, handling staggared U/V
                if (v.lower()=='u'):
                    temp_vals = (temp_vals[:,:,:,1:] + temp_vals[:,:,:,:-1])/2
                    temp_vals = temp_vals.swap_dims({'lon_u':'lon_x'})
                elif (v.lower()=='v'):
                    temp_vals = (temp_vals[:,:,1:,:] + temp_vals[:,:,:-1,:])/2
                    temp_vals = temp_vals.swap_dims({'lat_v':'lat_y'})

            ICAR[(v+'_m')] = temp_vals
            ICAR[(v+'_m')].attrs = temp_attrs

    ICAR = ICAR.drop(['u','v'])
    return ICAR

def unstaggerWRF(WRF):
        
    var_arr = []
        
    for v in WRF.data_vars:
        if (( len(WRF[v].dims) > 2 ) and (('stag' in WRF[v].dims[-1]) or ('stag' in WRF[v].dims[-2])) ):
            
            var_arr = var_arr + [v]
            
            print('unstaggering: '+v)
            temp_attrs = WRF[v].attrs
            new_dims = list(WRF[v].dims)
            new_dims[-1] = 'west_east'
            new_dims[-2] = 'south_north'
            new_dims = tuple(new_dims)
            
            x_stag = False
            if ('stag' in WRF[v].dims[-1]): x_stag = True
            
            temp_vals = WRF[v]
            
            if (len(WRF[v].dims) == 3):
                #Gather values, handling staggared U/V
                if (x_stag):
                    temp_vals = (temp_vals[:,:,1:] + temp_vals[:,:,:-1])/2
                    temp_vals = temp_vals.swap_dims({'west_east_stag':'west_east'})
                else:
                    temp_vals = (temp_vals[:,1:,:] + temp_vals[:,:-1,:])/2
                    temp_vals = temp_vals.swap_dims({'south_north_stag':'south_north'})
            elif (len(WRF[v].dims) == 4):
                #Gather values, handling staggared U/V
                if (x_stag):
                    temp_vals = (temp_vals[:,:,:,1:] + temp_vals[:,:,:,:-1])/2
                    temp_vals = temp_vals.swap_dims({'west_east_stag':'west_east'})
                else:
                    temp_vals = (temp_vals[:,:,1:,:] + temp_vals[:,:,:-1,:])/2
                    temp_vals = temp_vals.swap_dims({'south_north_stag':'south_north'})

            WRF[(v+'_m')] = temp_vals
            WRF[(v+'_m')].attrs = temp_attrs

    
    
    WRF = WRF.drop(var_arr)
    return WRF

def ds_bounds(ds,xvar,yvar):
    left = np.amin(ds[xvar].values[np.where(ds[xvar].values > 0)])
    bottom = np.amin(ds[yvar].values[np.where(ds[yvar].values > 0)])
    right = np.amax(ds[xvar].values[np.where(ds[xvar].values > 0)])
    top = np.amax(ds[yvar].values[np.where(ds[yvar].values > 0)])

    return [left, top, right, bottom]