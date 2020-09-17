import os
import sys
dirin = '/linkhome/rech/genimt01/uba22to/RT_CNES/OSE/eval_OSE_tools/OSE_evaluations'
sys.path.append(dirin)  
sys.path.append(dirin+'/scuba/src')
  
import pyinterp
import netCDF4
import numpy as np
import matplotlib.pylab as plt
import glob
import proplot as plot

from yaml import load, Loader
YAML = load(open(str(dirin+'/scuba_alongtrack_DUACS.yaml')), Loader=Loader)
mission_management = dirin+'/scuba/share/MissionManagement.yaml'
import scipy.signal

from scuba.src.mod_constant import *
from scuba.src.mod_segmentation import *
from scuba.src.mod_spectral import *

from src.mod_mission import *
from src.mod_write import *
from src.mod_filter import *
from src.mod_grid import *
from src.mod_alongtrack import *

# function to create recursive paths
def mk_dir_recursive(dir_path):
    if os.path.isdir(dir_path):
        return
    h, t = os.path.split(dir_path)  # head/tail
    if not os.path.isdir(h):
        mk_dir_recursive(h)

    new_path = join_paths(h, t)
    if not os.path.isdir(new_path):
        os.mkdir(new_path)

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

# setup learning configuration
lag        = sys.argv[1]
domain     = sys.argv[2]
wregul     = str2bool(sys.argv[3])
scratchpath = '/gpfsscratch/rech/yrf/uba22to/DINAE_keras/OSE/'+domain
if wregul == True:
    scratchpath = scratchpath+'/resIA_nadir_nadlag_'+lag+'_obs/FP_GENN_wwmissing_wOI_wtrain_wregul'
else:
    scratchpath = scratchpath+'/resIA_nadir_nadlag_'+lag+'_obs/FP_GENN_wwmissing_wOI_wtrain'

# setup c2 alongtrack configuration
residus_files_cmems = '/gpfswork/rech/yrf/uba22to/RAW_DATA/OSE/alongtracks/c2/*'
mission = 'c2'
# study area
lon_min = 295.   # -65
lon_max = 305    # -55
lon_min = 0.
lon_max = 360.
lat_min = 33.
lat_max = 43.
is_circle = False
time_min = '2017-01-01'
time_max = '2017-12-31'
bin_lat_step = 1.
bin_lon_step = 1.
bin_time_step = '1D'

lvar=['DUACS','FP_GENN']
lvar2=['OI','FP-GENN']
output_filename_stat = [ 'stat_OSE_GULFSTREAM_'+lvar[i]+'_'+time_min+'-'+time_max+'_'+mission+'.nc' for i in range(len(lvar)) ]
output_filename_psd = ['psd_OSE_GULFSTREAM_'+lvar[i]+'_'+time_min+'-'+time_max+'_'+mission+'_vxxc.nc' for i in range(len(lvar)) ]

# Read along-track
ds_alongtrack = cmems_alongtrack_dataset_v2(residus_files_cmems, 
                                      lon_min=lon_min, lon_max=lon_max, 
                                      lat_min=lat_min, lat_max=lat_max, 
                                      time_min=time_min, time_max=time_max)

for i in range(len(lvar)):

    print('Processing '+lvar[i]+'...')

    # import Interpolation maps
    grid_files = scratchpath+'/OSE_GULFSTREAM_FPGENN.nc'
    var2add = [lvar2[i]] #  
    var2sub = []

    x_axis, y_axis, z_axis, grid = fpgenn_grid_dataset(grid_files, var2add, var2sub,
                                                       lon_min=lon_min, lon_max=lon_max, 
                                                       lat_min=lat_min, lat_max=lat_max, 
                                                       time_min=time_min, time_max=time_max,
                                                       is_circle=is_circle)


    # Interpolate maps on alongtrack
    map_interp = pyinterp.trivariate(grid, 
                                 ds_alongtrack["longitude"].values, 
                                 ds_alongtrack["latitude"].values,
                                 z_axis.safe_cast(ds_alongtrack.time.values),
                                 bounds_error=False).reshape(ds_alongtrack["longitude"].values.shape)
    ssh_alongtrack = ds_alongtrack["SLA"].values + ds_alongtrack["mdt"].values
    lon_alongtrack = ds_alongtrack["longitude"].values
    lat_alongtrack = ds_alongtrack["latitude"].values
    time_alongtrack = ds_alongtrack["time"].values

    # get mask from map_interp
    msk1 = np.ma.masked_invalid(ssh_alongtrack).mask  
    msk2 = np.ma.masked_invalid(map_interp).mask
    msk = msk1+msk2
    ssh_alongtrack = np.ma.masked_where(msk, ssh_alongtrack).compressed()
    lon_alongtrack = np.ma.masked_where(msk, lon_alongtrack).compressed()
    lat_alongtrack = np.ma.masked_where(msk, lat_alongtrack).compressed()
    time_alongtrack = np.ma.masked_where(msk, time_alongtrack).compressed()
    map_interp = np.ma.masked_where(msk, map_interp).compressed()

    # save statistics in NetCDF file
    ncfile = netCDF4.Dataset(output_filename_stat[i],'w')
    binning = pyinterp.Binning2D(
        pyinterp.Axis(np.arange(0, 360, bin_lon_step), is_circle=True),
        pyinterp.Axis(np.arange(-90, 90 + bin_lat_step, bin_lat_step)))
    # binning alongtrack
    binning.push(lon_alongtrack, lat_alongtrack, ssh_alongtrack, simple=True)
    write_stat(ncfile, 'alongtrack', binning)
    binning.clear()
    # binning map interp
    binning.push(lon_alongtrack, lat_alongtrack, map_interp, simple=True)
    write_stat(ncfile, 'maps', binning)
    binning.clear()
    # binning diff sla-msla
    binning.push(lon_alongtrack, lat_alongtrack, ssh_alongtrack - map_interp, simple=True)
    write_stat(ncfile, 'diff', binning)
    binning.clear()
    # add rmse
    diff2 = (ssh_alongtrack - map_interp)**2
    binning.push(lon_alongtrack, lat_alongtrack, diff2, simple=True)
    var = ncfile.groups['diff'].createVariable('rmse', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T  
    ncfile.close()

    # compute dsp
    # make time vector as days since 1950-01-01
    time_alongtrack = (time_alongtrack - np.datetime64('1950-01-01T00:00:00Z')) / np.timedelta64(1, 'D')

    ref_segment, lon_segment, lat_segment, delta_x, npt, study_segment = compute_segment_alongtrack(ssh_alongtrack,
                                                                                                lon_alongtrack,
                                                                                                lat_alongtrack,
                                                                                                time_alongtrack,
                                                                                                YAML,
                                                                                                map_interp)
    # Power spectrum density reference field
    global_wavenumber, global_psd_ref = scipy.signal.welch(np.asarray(ref_segment).flatten(),
                                                           fs=1.0 / delta_x,
                                                           nperseg=npt,
                                                           scaling='density',
                                                           noverlap=0)
    # Power spectrum density study field
    _, global_psd_study = scipy.signal.welch(np.asarray(study_segment).flatten(),
                                             fs=1.0 / delta_x,
                                             nperseg=npt,
                                             scaling='density',
                                             noverlap=0)
    # Power spectrum density study field
    _, global_psd_diff = scipy.signal.welch(np.asarray(study_segment).flatten()-np.asarray(ref_segment).flatten(),
                                             fs=1.0 / delta_x,
                                             nperseg=npt,
                                             scaling='density',
                                             noverlap=0)

    # Save in NetCDF file
    ds = xr.Dataset({"psd_ref": (["wavenumber"], global_psd_ref),
                 "psd_study": (["wavenumber"], global_psd_study),
                 "psd_diff": (["wavenumber"], global_psd_diff),
                },
                coords={"wavenumber": (["wavenumber"], global_wavenumber),
                       },
               )
    ds.to_netcdf(output_filename_psd[i])

## plot RMSE  
ds_DUACS = xr.open_dataset(output_filename_stat[0], group='diff')
ds_FP_GENN =  xr.open_dataset(output_filename_stat[1], group='diff')
plot.rc['savefig.transparent']=False
fig, ax = plot.subplots(ncols=2, nrows=1, span=False, axwidth='4cm')
ax.format(xlabel='longitude', ylabel='latitude', grid=False, suptitle='RMSE in 1$^{\circ}$x1$^{\circ}$ grid boxes (Cryosat-2 independent)', abc=True, abcloc='ul', xlim=(295, 305), ylim=(33, 43))
ax[0].pcolormesh(ds_DUACS['lon'].values, ds_DUACS['lat'].values, ds_DUACS['rmse'], vmin=0, vmax=0.2)
ax[0].format(title='DUACS OI')
pc = ax[1].pcolormesh(ds_FP_GENN['lon'].values, ds_FP_GENN['lat'].values, ds_FP_GENN['rmse'], vmin=0, vmax=0.2)
ax[1].format(title='FP-GENN')
cbar = fig.colorbar(pc, label='RMSE [m]', ticks=0.02, loc='r', length=0.7, formatter='simple')
plt.savefig(scratchpath+'/RMSE.pdf')

## plot Spectral analysis
ds_DUACS = xr.open_dataset(output_filename_psd[0])
ds_FP_GENN =  xr.open_dataset(output_filename_psd[1])
plot.rc['savefig.transparent']=False
fig, axs = plot.subplots(ncols=2, share=False)
axs.format(xlabel='wavelength [km]', grid=True, abc=True, abcloc='ur', xlim=(50, 500), xscale='log', xlocator=[50, 100, 150, 200, 300, 400, 500])
p1 = axs[0].plot(1./ds_DUACS.wavenumber, ds_DUACS.psd_ref.values, label='Cryosat-2', color='k', lw=3)#, legend='ll', legend_kw={'ncol': 1})
p2 = axs[0].plot(1./ds_DUACS.wavenumber, ds_DUACS.psd_study.values, label='DUACS', lw=6)#, legend='ll', legend_kw={'ncol': 1})
p3 = axs[0].plot(1./ds_FP_GENN.wavenumber, ds_FP_GENN.psd_study.values, label='FP-GENN', lw=2)#, legend='ll', legend_kw={'ncol': 1})
fig.legend([p1, p2, p3], ncols=3, loc='t')
axs[0].format(ylabel='PSD [m$^{2}$cy$^{-1}$km$^{-1}$]', yscale='log')

c2 = plot.scale_luminance('red', 0.5)
axs[1].plot(1./ds_DUACS.wavenumber, 1. - ds_DUACS.psd_diff/ds_DUACS.psd_ref, label='DUACS', lw=6)#, legend='ll', legend_kw={'ncol': 1})
axs[1].plot(1./ds_FP_GENN.wavenumber, 1. - ds_FP_GENN.psd_diff/ds_FP_GENN.psd_ref,  label='FP-GENN', lw=2)#, legend='ll', legend_kw={'ncol': 1})
axs[1].format(ylabel='1.0 - Ratio PSD$_{err}$/PSD$_{c2}$', ylim=(0, 1))
axs[1].hlines(y=0.5, xmin=20, xmax=500, lw=1, color=c2)
plt.savefig(scratchpath+'/DSP.pdf')

