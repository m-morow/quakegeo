import rioxarray
import rasterio
import numpy as np
import xarray as xr
import os
import tables as tb #pytables
import h5py
import datetime

"""
---- Optimized version, using Pytables (chunked arrays) to reduce memory and processing time ----

Convert Sentinel-1 radar data into Dx, Dy, Dz

Pre-processing steps required (in QGIS or other):
- resample .tif to same resolution and convert to UTM
- extend .tif to cover the same area (set non-data to some identifiable value)

functions {geotiff2xarray, geotiff2look} were written by H Yin
functions {make_Gi, make_di, invert, get_mask, make_mask, hdf2raster} were written by M Tan
"""

def geotiff2xarray(geotiff_file):
    '''
    Parameters
    ----------
    geotiff_file : string
        full path to .tif file that you want to convert
        This code was built to handle the geotiffs produced by ISCE for range or azimuth SAR offsets but might be able to handle other types of geotiffs
    Returns
    -------
    Returns an xarray dataset. This must be flattened to a dataarray to be compatible with pygmt

    NOTES: Might need to be updated to handle LOS datasets
    '''
    da = rioxarray.open_rasterio(geotiff_file)          # Open the data array with rioxarray
    observed_ds = da.to_dataset('band')          # Take the band object and set it as a variable
    observed_ds = observed_ds.rename({1: 'los_obs','x': 'lon','y': 'lat'})          # Rename the variable to a more useful name

    # NaNs are read in as -140000 by default so let's convert all of these to NaNs 
    observed_ds = observed_ds.where(observed_ds != -140000)          # replace all values equal to -140000 with NaNs (azimuth)
    observed_ds = observed_ds.where(observed_ds != -23000)           # replace all values equal to -24000 with NaNs (range)

    return observed_ds

def geotiff2look(geotiff_file):

    '''
    Parameters
    ----------
    geotiff_file : string
        full path to .tif Look file. 
        This should be a *los.tif file which contains look vectors as bands 1
        This code was built to handle the los (look) geotiffs produced by ISCE
    Returns
    -------
    look_ds : xarray dataset object
        With additional data variables:
            * The LOS look vector (los_e, los_n, los_u)
            * The azimuth offset look vector (az_e, az_n, az_u)
    Ex: look_ds = geotiff2look('s1_2023018-20230209_p14_los.tif')
    '''

    # ## Calculate Average range Angle
    # ## Calculate Average azimuth Angle
    da = rioxarray.open_rasterio(geotiff_file)          # Open the dataset with rioxarray
    look_ds = da.to_dataset('band')                  # Covert our xarray.DataArray into a xarray.Dataset
    look_ds = look_ds.rename({1: 'range',2:'azimuth'})          # Rename the variable to a more useful name
    look_ds = look_ds.where(look_ds != 0)          # replace zeros with NaNs 

    # Print mean range and azimuth angles for user
    print('Mean range angle (deg): ', np.nanmean(look_ds.range))
    print('Mean azimuth angle (deg): ', np.nanmean(look_ds.azimuth))

    # Calculate LOS Look vector
    #### MT: HYin E and N component calculations switched in Bill's doc
    #### f=lambda a: a * a ---> {function object stores result} = {keyword} {argument:} {expression}
    # Calculate East component and write it to our dataset as a new variable
    look_ds = look_ds.assign(los_e=lambda look_ds: np.sin(np.radians(look_ds.azimuth))*np.sin(np.radians(look_ds.range))) #Bill's doc
    #look_ds = look_ds.assign(los_e=lambda look_ds: np.cos(np.radians(look_ds.azimuth)+np.radians(450))*np.sin(np.radians(look_ds.range))) #HYin
    # Calculate North component and write it to our dataset as a new variable
    look_ds = look_ds.assign(los_n=lambda look_ds: np.cos(np.radians(look_ds.azimuth))*np.sin(np.radians(look_ds.range))) #Bill's doc
    #look_ds = look_ds.assign(los_n=lambda look_ds: np.sin(np.radians(look_ds.azimuth)+np.radians(450))*np.sin(np.radians(look_ds.range))) #HYin
    # Calculate Vertical component and write it to our dataset as a new variable
    look_ds = look_ds.assign(los_u=lambda look_ds: -1*np.cos(np.radians(look_ds.range))) #Bill's doc
    #look_ds = look_ds.assign(los_u=lambda look_ds: np.cos(np.radians(look_ds.range))) #HYin
    #print('Mean LOS look vector (e, n, u): ( ',np.nanmean(look_ds.los_e),',' ,np.nanmean(look_ds.los_n),',',np.nanmean(look_ds.los_u),')')
    #print('**Note: Positive values indicate motion in the direction of the satellite')
    #print('This is the same as the range offset look vector')

    # Calculate Azimuth offset Look vector
    # Calculate East component and write it to our dataset as a new variable
    look_ds = look_ds.assign(az_e=lambda look_ds: np.sin(np.radians(look_ds.azimuth))) #Bill's doc
    #look_ds = look_ds.assign(az_e=lambda look_ds: np.cos(np.radians(look_ds.azimuth))) #HYin
    look_ds = look_ds.assign(az_n=lambda look_ds: np.cos(np.radians(look_ds.azimuth))) #Bill's doc
    #look_ds = look_ds.assign(az_n=lambda look_ds: np.sin(np.radians(look_ds.azimuth))) #HYin
    look_ds = look_ds.assign(az_u=lambda look_ds: (look_ds.azimuth)*0) #HYin, Bill
    #print('Mean Azimuth offset look vector (e, n, u): ( ',np.nanmean(look_ds.az_e),',' ,np.nanmean(look_ds.az_n),',',np.nanmean(look_ds.az_u),')')              
    #print('**Note: Positive values indicate motion in the along-track direction of the satellite (i.e. the direction that the satellite is moving is positive)')

    # Reformat the x and y data coordinates to lat and lon
        # Renames the variable within variable ra_{}, az_{}, azimuth, range to a more useful name
        # ex: az_e: <xarray.DataArray 'az_e' (lat: 1064, lon: 543)>
    look_ds = look_ds.rename({'x': 'lon','y': 'lat'}) 
    look_ds = look_ds.drop(labels=['range', 'azimuth'])         
 
    return look_ds

def make_Gi(i, look_b1, look_b2, look_b3):
    """
    builds "G" matrix at pixel, i
    input:
        i, pixel
        look_b1, los path (in this case, p14)
        look_b2, los path (in this case, p21)
        look_b3, los path (in this case, p116)
    output:
        6x3 matrix of float and/or NaN values
    """
    Gmat = np.zeros(shape=(6,3))

    #path 1
    Gmat[0, 0] = float(look_b1.data_vars['los_e'][i])
    Gmat[0, 1] = float(look_b1.data_vars['los_n'][i]) #idx=379020 has first non NaN value
    Gmat[0, 2] = float(look_b1.data_vars['los_u'][i])
    Gmat[1, 0] = float(look_b1.data_vars['az_e'][i])
    Gmat[1, 1] = float(look_b1.data_vars['az_n'][i])
    Gmat[1, 2] = float(look_b1.data_vars['az_u'][i])

    #path 2
    Gmat[2, 0] = float(look_b2.data_vars['los_e'][i])
    Gmat[2, 1] = float(look_b2.data_vars['los_n'][i]) #idx=379020 has first non NaN value
    Gmat[2, 2] = float(look_b2.data_vars['los_u'][i])
    Gmat[3, 0] = float(look_b2.data_vars['az_e'][i])
    Gmat[3, 1] = float(look_b2.data_vars['az_n'][i])
    Gmat[3, 2] = float(look_b2.data_vars['az_u'][i])

    #path 3
    Gmat[4, 0] = float(look_b3.data_vars['los_e'][i])
    Gmat[4, 1] = float(look_b3.data_vars['los_n'][i]) #idx=379020 has first non NaN value
    Gmat[4, 2] = float(look_b3.data_vars['los_u'][i])
    Gmat[5, 0] = float(look_b3.data_vars['az_e'][i])
    Gmat[5, 1] = float(look_b3.data_vars['az_n'][i])
    Gmat[5, 2] = float(look_b3.data_vars['az_u'][i])
    
    Gmat[Gmat > 3.4e38] = np.nan #manually set 3.403e+38 to NaN
    
    return Gmat

def invert(Gmat, dmat):
    #on the cluster, code failed with error "linalg svd did not converge"
    #included exception
    # A * B where   
    # (G^T * G)^-1 (G^T * d)
    try:
        result = np.linalg.lstsq(Gmat, dmat, rcond=None)[0]     
    except:
        result = np.array([[-9999],[-9999],[-9999]])

    return result
    
def make_di(i, az_b1, ra_b1, az_b2, ra_b2, az_b3, ra_b3):
    """
    builds "d" matrix at pixel, i
    input:
        i, pixel
        az_b{}, azimuth raster for path
        ra_b{}, range raster for path
    output:
        6x1 matrix of float and/or NaN values
    """
    dmat = np.zeros(shape=(6,1))

    dmat[0, 0] = float(ra_b1.data_vars['rng_offset'][i]) #idx=315468 has first non NaN value
    dmat[1, 0] = float(az_b1.data_vars['az_offset'][i])

    dmat[2, 0] = float(ra_b2.data_vars['rng_offset'][i]) #idx=315468 has first non NaN value
    dmat[3, 0] = float(az_b2.data_vars['az_offset'][i])

    dmat[4, 0] = float(ra_b3.data_vars['rng_offset'][i]) #idx=315468 has first non NaN value
    dmat[5, 0] = float(az_b3.data_vars['az_offset'][i])
    
    dmat[dmat > 3.4e38] = np.nan #manually set 3.403e+38 to NaN

    return dmat

def get_mask(d):
    """
    grabs the indices of NaN values in matrix d
    input:
        d (array)
    output:
        column array of indices that correspond to the rows in d that have an NaN value
    """
    maskidx = np.argwhere(np.isnan(d))

    return np.unique(maskidx[:,0])

def make_mask(maskidx, Gmat, dmat):
    """
    deletes rows in matrix that correspond to the NaN rows in matrix d
    input:
        maskidx, array of row indices to remove
        matrix array
    output:
        matrix with deleted rows (this is the final G and d matrices)
    """
    maskidx = list(maskidx)
    G_res = np.delete(Gmat, maskidx, 0)
    d_res = np.delete(dmat, maskidx, 0)
    
    if len(np.argwhere(np.isnan(G_res))) == 0: #checks if G_res has its own NaN values
        return G_res, d_res #deletes row based on maskidx list
    else:
        mask = np.unique(np.argwhere(np.isnan(G_res))[:,0])
        return np.delete(G_res, mask, 0), np.delete(d_res, mask, 0)

def hdf2raster(metadata, h5f, setname):
    """
    converts H5 file to geotiff using reference raster metadata
    input:
        metadata, reference raster (I used azimuth)
        h5f, filepath as string
        setname, pytables identifier for dataset
                 ex: x = Dx_h5.create_carray(rootx, 'x'), where 'x' is the setname
    output:
        geotiff, saved to the same filepath as h5f
    """
    with h5py.File(h5f, 'r') as hf:
        ds_arr = hf[setname][()]

    outraster = h5f[:-3] + '.tif' #use h5f filepath and filename to create geotiff raster
    with rasterio.open(outraster, mode='w',**metadata) as dst:
        dst.write(ds_arr.astype(rasterio.float32), indexes=1) # the number one is the number of bands
    
    dst.close()

##############################
##    User defined inputs   ## 
##############################

#working directories
#dir = '/Users/mtan/Documents/Data/Turkey/invertDispFields/UTM/extended/' #input raster directory
dir = '/home/mtan/data/turkey/invertDisp'
#hdir = '/Users/mtan/Documents/Data/Turkey/invertDispFields/UTM/clusterResults' #where H5 files and converted geotiffs will be saved
hdir = '/home/mtan/data/turkey/invertDisp/results'
logfile = 'log_5363x4299.txt' #code writes a log file (code start&end time, output raster dims)

#line of sight rasters
los_1 = 'p14_los_test4cluster.tif'
los_2 = 'p21_los_test4cluster.tif'
los_3 = 'p116_los_test4cluster.tif'

#azimuth rasters
azi_1 = 'p14_azimuth_test4cluster.tif'
azi_2 = 'p21_azimuth_test4cluster.tif'
azi_3 = 'p116_azimuth_test4cluster.tif'

#range rasters
ran_1 = 'p14_range_test4cluster.tif'
ran_2 = 'p21_range_test4cluster.tif'
ran_3 = 'p116_range_test4cluster.tif'

#HD5F filenames to save to
Dxh5 = 'Dx-test4cluster-5363x4299v2.h5'
Dyh5 = 'Dy-test4cluster-5363x4299v2.h5'
Dzh5 = 'Dz-test4cluster-5363x4299v2.h5'

#pixel resolution of input rasters
ncol = 4299 #manually set for how many columns will be evaluated in code, AKA pixel width
nrow = 5363 #AKA pixel height
pix = nrow * ncol #total pixels evaluated

#ask user to check if inputs are correct
print("Are the following inputs correct?")
print("     raster dir   = {}".format(dir))
print("     hdf5 dir     = {}".format(hdir))
print("     out tif res. = {} x {} (width x height)".format(ncol, nrow))
print("\n")
print("     los raster 1 = {}".format(los_1))
print("     los raster 2 = {}".format(los_2))
print("     los raster 3 = {}".format(los_3))
print("\n")
print("     azi raster 1 = {}".format(azi_1))
print("     azi raster 2 = {}".format(azi_2))
print("     azi raster 3 = {}".format(azi_3))
print("\n")
print("     ran raster 1 = {}".format(ran_1))
print("     ran raster 2 = {}".format(ran_2))
print("     ran raster 3 = {}".format(ran_3))
print("\n")
print("     HD5F Dx file = {}".format(Dxh5))
print("     HD5F Dy file = {}".format(Dyh5))
print("     HD5F Dz file = {}".format(Dzh5))

if input("Press 'y' for YES or any key for NO") != "y":
    print("     Code terminated.")
    exit()
else:
    print("Reading rasters and initializing outfiles ...")

##############################
##      Read in rasters     ## 
##############################

#read in los rasters----
look_ds_Gi = geotiff2look(os.path.join(dir, los_1))
look_ds_Gi2 = geotiff2look(os.path.join(dir, los_2))
look_ds_Gi3 = geotiff2look(os.path.join(dir, los_3))
#-----------------------

#read in range and azimuth rasters----
azimuth_offsets_ds = geotiff2xarray(os.path.join(dir, azi_1))
azimuth_offsets_ds = azimuth_offsets_ds.rename({'los_obs':'az_offset'})
range_offsets_ds = geotiff2xarray(os.path.join(dir, ran_1))
range_offsets_ds = range_offsets_ds.rename({'los_obs':'rng_offset'})

azimuth_offsets_ds2 = geotiff2xarray(os.path.join(dir, azi_2))
azimuth_offsets_ds2 = azimuth_offsets_ds2.rename({'los_obs':'az_offset'})
range_offsets_ds2 = geotiff2xarray(os.path.join(dir, ran_2))
range_offsets_ds2 = range_offsets_ds2.rename({'los_obs':'rng_offset'})

azimuth_offsets_ds3 = geotiff2xarray(os.path.join(dir, azi_3))
azimuth_offsets_ds3 = azimuth_offsets_ds3.rename({'los_obs':'az_offset'})
range_offsets_ds3 = geotiff2xarray(os.path.join(dir, ran_3))
range_offsets_ds3 = range_offsets_ds3.rename({'los_obs':'rng_offset'})
#--------------------------------------

#stack Xarray datasets on the lon/lat dimension----
stackedGi = look_ds_Gi.stack(gridcells=('lon','lat'))
stackedGi2 = look_ds_Gi2.stack(gridcells=('lon','lat'))
stackedGi3 = look_ds_Gi3.stack(gridcells=('lon','lat'))

stacked_partial_d_ra = range_offsets_ds.stack(gridcells=('lon','lat'))
stacked_partial_d_az = azimuth_offsets_ds.stack(gridcells=('lon','lat'))

stacked_partial_d_ra2 = range_offsets_ds2.stack(gridcells=('lon','lat'))
stacked_partial_d_az2 = azimuth_offsets_ds2.stack(gridcells=('lon','lat'))

stacked_partial_d_ra3 = range_offsets_ds3.stack(gridcells=('lon','lat'))
stacked_partial_d_az3 = azimuth_offsets_ds3.stack(gridcells=('lon','lat'))
#---------------------------------------------------

#####################################
##     Set up for output files     ## 
#####################################

nchunk = 0 #begin counter for column

#initialize arrays to temporarily hold data before writing to h5f files
D_x = []
D_y = []
D_z = []

#initialize HDF5 files to write LSQ inversion to----
filters = tb.Filters(complevel=5) #compression level for H5 files where 0=none/default, 9=highest

Dx_h5 = tb.open_file(os.path.join(hdir, Dxh5), mode='w', title='x array')
rootx = Dx_h5.root
x = Dx_h5.create_carray(rootx, 'x', tb.Float32Atom(), shape=(nrow, ncol), filters=filters) #'x' is setname
#shape=(# elements in column, # elements in row)

Dy_h5 = tb.open_file(os.path.join(hdir, Dyh5), mode='w', title='y array')
rooty = Dy_h5.root
y = Dy_h5.create_carray(rooty, 'y', tb.Float32Atom(), shape=(nrow, ncol), filters=filters) #'y' is setname

Dz_h5 = tb.open_file(os.path.join(hdir, Dzh5), mode='w', title='z array')
rootz = Dz_h5.root
z = Dz_h5.create_carray(rootz, 'z', tb.Float32Atom(), shape=(nrow, ncol), filters=filters) #'z' is setname
#----------------------------------------------------

#get metadata from reference tif file----
inputpath_name = os.path.join(dir, azi_1)
raster = rasterio.open(inputpath_name)
out_meta = raster.profile

#explicitly set parameters
out_meta.update({'driver':'GTiff',
                 'width':ncol, #ncol or raster.shape[1]
                 'height':nrow, #nrow or raster.shape[0]
                 'dtype':'float32',
                 'count': 1, #bands
                 'crs':raster.crs, 
                 'transform':raster.transform,
                 'nodata':-9999})
#-----------------------------------------

#####################################
##      Begin of calculations      ## 
#####################################

print("Begin calculations ...") #this will take the longest time
begin = datetime.datetime.now() #get code starting time

for i in range(0, pix, 1):
    Gi = make_Gi(i, stackedGi, stackedGi2, stackedGi3)
    di = make_di(i, stacked_partial_d_az, stacked_partial_d_ra, stacked_partial_d_az2, stacked_partial_d_ra2, stacked_partial_d_az3, stacked_partial_d_ra3)
    di_mask = get_mask(di)
    Gi_masked, di_masked = make_mask(di_mask, Gi, di)
    result = invert(Gi_masked, di_masked)
    D_x = np.append(D_x, result[0]) 
    D_y = np.append(D_y, result[1])
    D_z = np.append(D_z, result[2])

    if len(D_x) == nrow: #if 'D_x' == column length
        x[:, nchunk] = D_x.reshape(nrow) #reshapes array to column of nrow elements long
        y[:, nchunk] = D_y.reshape(nrow)
        z[:, nchunk] = D_z.reshape(nrow)
        nchunk += 1 #identifies no. columns evaluated
        D_x, D_y, D_z = [], [], [] #reset arrays

#close h5f files
Dx_h5.close()
Dy_h5.close()
Dz_h5.close()

print("     ... code successfully wrote to h5 files.")

#####################################
##   Convert H5 files to Geotiffs  ## 
#####################################

print("Writing to rasters ...")

#Dx
hdf2raster(out_meta, h5f=os.path.join(hdir, Dxh5), setname='x')

#Dy
hdf2raster(out_meta, h5f=os.path.join(hdir, Dyh5), setname='y')

#Dz
hdf2raster(out_meta, h5f=os.path.join(hdir, Dzh5), setname='z')

end = datetime.datetime.now() #get code ending time

#write to logfile
with open(os.path.join(hdir, logfile), 'w') as f: 
    f.write('begin: {}'.format(str(begin)))
    f.write('\n')
    f.write('end: {}'.format(str(end)))
    f.write('\n')
    f.write('ncol x nrow: {} x {}'.format(ncol, nrow))
    f.write('\n')
    f.write('ending i iterator: {}'.format(i))
    f.write('\n')
f.close()

print("Code completed.")
