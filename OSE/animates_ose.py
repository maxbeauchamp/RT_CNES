#!/usr/bin/env python

""" display_results.py: script to display maps, taylors diagrams and radially averaged power spectrums from MS-PB-AnDA and MS-VE-DINEOF with three different datasets (nadir / swot / nadirswot) """

__author__ = "Maxime Beauchamp"
__version__ = "0.1"
__date__ = "2019-12-10"
__email__ = "maxime.beauchamp76@gmail.com"

from pb_anda import *
from netCDF4 import Dataset
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
np.random.seed(1)

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

lag        = sys.argv[1]
domain     = sys.argv[2]
wregul     = str2bool(sys.argv[3])

scratchpath = '/gpfsscratch/rech/yrf/uba22to/4DVARNN-DINAE/OSE/'+domain
if wregul == True:
    scratchpath = scratchpath+'/resIA_nadir_nadlag_'+lag+'_obs/GB1_GENN_wwmissing_wOI_wtrain_wregul'
else:
    scratchpath = scratchpath+'/resIA_nadir_nadlag_'+lag+'_obs/GB1_GENN_wwmissing_wOI_wtrain'

# Reload saved GENN result
file_results_GENN=scratchpath+'/saved_path_004_GENN_wwmissing.pickle'
with open(file_results_GENN, 'rb') as handle:
    GT,  Obs_nadir, itrp_GENN, rec_GENN, itrp_OI = pickle.load(handle)

			#*****************#
			# Display results #
			#*****************#

if domain=="OSMOSIS":
    extent     = [-19.5,-11.5,45.,55.]
    extent2    = [-19.25,-11.5,45.25,55.]
    indLat     = 200
    indLon     = 160
elif domain=='GULFSTREAM':
    extent     = [-65.,-55.,33.,43.]
    extent2    = [-64.75,-55.,33.25,43.]
    indLat     = 200
    indLon     = 200
else:
    extent=[-65.,-55.,30.,40.]
    indLat     = 200
    indLon     = 200
lon = np.arange(extent[0],extent[1],1/20)
lat = np.arange(extent[2],extent[3],1/20)
lon = lon[:indLon]
lat = lat[:indLat]

file_validation = "/gpfswork/rech/yrf/uba22to/DATA/OSE/"+domain+"/validation/dataset_nadir_0d.nc"
nc_data_valid   = Dataset(file_validation,'r')
data_valid      = np.copy(nc_data_valid['ssh'][:,:indLat,:indLon])
nc_data_valid.close()


# list of dates
lday=[ datetime.strftime(datetime.strptime("2017-01-01",'%Y-%m-%d')\
                          + timedelta(days=i),"%Y-%m-%d") for i in range(365) ]
lday2 = [ datetime.strptime(lday[i],'%Y-%m-%d') for i in range(len(lday)) ]

## animations (videos) ##
def animate(i):
    print(i)
    # Load data
    OI                          = itrp_OI[i]
    Grad_OI                     = Gradient(OI,2)
    GB_GENN                     = itrp_GENN[i]
    #GB_GENN                     = median_filter(GB_GENN, size=5)
    Grad_GB_GENN                = Gradient(GB_GENN,2)
    for ivar in range(len(var)):
        plot(ax,0,ivar,lon,lat,eval(var[ivar]),title[ivar],\
                extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
    return ax[0][0], ax[0][1]

var=['Grad_OI','Grad_GB_GENN']
title=[r"$\nabla_{OI}$",r"$\nabla_{GB-GENN}$"]
vmax=0.1
vmin=0.
cmap="viridis"
fig, ax = plt.subplots(1,2,figsize=(22,10),squeeze=False,\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, blit=False, frames=np.arange(1,365), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=10, bitrate=5000)
ani.save(scratchpath+"/animation_grads.mp4", writer = writer)
plt.close()


