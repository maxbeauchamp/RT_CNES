#!/usr/bin/env python

from pb_anda import *
from netCDF4 import Dataset
import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from skimage.util.shape import view_as_blocks
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

domain   = sys.argv[1] 
workpath    = "/users/local/m19beauc/4DVARNN-DinAE_xp/OSSE_keras/4Ronan/"+domain
datapath    = '/gpfswork/rech/yrf/uba22to/DATA/'
'''if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
else:
    shutil.rmtree(workpath)
    mk_dir_recursive(workpath)
'''

# 4-plots video individual maps
if domain=="OSMOSIS":
    extent     = [-19.5,-11.5,45.,55.]
    indLat     = 200
    indLon     = 160
elif domain=='GULFSTREAM':
    extent     = [-65.,-55.,33.,43.]
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

# reload GT and obs from original datasets
fileMod=datapath+domain+"/ref/NATL60-CJM165_"+domain+"_ssh_y2013.1y.nc"
nc_data_mod = Dataset(fileMod,'r')
GT  = np.copy(nc_data_mod['ssh'][:,:indLat,:indLon])
print(GT.shape)
nc_data_mod.close()

# Compute GT vorticity
Vorticity_GT=np.zeros(GT.shape)
for i in range(len(GT)):
    Vorticity_GT[i] = cv2.Laplacian(GT[i],cv2.CV_64F)

def animate(i):
    # Load data
    print(i)
    ax[0][0].clear()
    plot(ax,0,0,lon,lat,Vorticity_GT[i],"Vorticity",\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
    return ax[0][0]

dr=range(0,365)
vmax=np.round(np.nanmax([np.abs(Vorticity_GT[dr]),np.abs(Vorticity_GT[dr])]),2)
vmin=-1.0*vmax
cmap="RdBu_r"

# domain 1/20
# animation
dr=range(0,365)
fig, ax = plt.subplots(1,1,figsize=(15,15),squeeze=False,\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
ani = animation.FuncAnimation(fig, animate, frames=dr, blit=False, interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=3, bitrate=5000)
ani.save(workpath+"/animation_vorticity_"+domain+"_1_20.mp4", writer = writer)
plt.close()

# plot mean(|q(t)-q(0)|**2)/var(q(0)
lday=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=151+i),"%Y-%m-%d") for i in range(32) ]
dr=range(151,183)
nmq=[np.mean((Vorticity_GT[x]-Vorticity_GT[dr[0]])**2)/np.var(Vorticity_GT[dr[0]]) for x in dr]
fig, ax = plt.subplots()
plt.plot(lday,nmq)
plt.gcf().autofmt_xdate()
plt.ylabel(r'$\frac{\overline{(q_t-q_0)^2}}{\sigma^2(q_0)}$')
ax.yaxis.label.set_size(16)
plt.savefig(workpath+'/Normalized_mean_vorticity_deviations_MARCH_1_20.png',  bbox_inches='tight')

lday=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=92+i),"%Y-%m-%d") for i in range(32) ]
dr=range(92,124)
nmq=[np.mean((Vorticity_GT[x]-Vorticity_GT[dr[0]])**2)/np.var(Vorticity_GT[dr[0]]) for x in dr]
fig, ax = plt.subplots()
plt.plot(lday,nmq)
plt.gcf().autofmt_xdate()
plt.ylabel(r'$\frac{\overline{(q_t-q_0)^2}}{\sigma^2(q_0)}$')
ax.yaxis.label.set_size(16)
plt.savefig(workpath+'/Normalized_mean_vorticity_deviations_JANUARY_1_20.png', bbox_inches='tight')

# domain 1/5
block_shape = (1,4,4)
Vorticity_GT_s4 = view_as_blocks(Vorticity_GT[:,:,:], block_shape)
Vorticity_GT_s4 = Vorticity_GT_s4.reshape(Vorticity_GT_s4.shape[0], Vorticity_GT_s4.shape[1], Vorticity_GT_s4.shape[2], -1)
Vorticity_GT_s4 = np.mean(Vorticity_GT_s4, axis=3)

index = np.arange(0,len(lon),4)
lon = lon[index]
index = np.arange(0,len(lat),4)
lat = lat[index]
def animate2(i):
    # Load data
    print(i)
    ax[0][0].clear()
    plot(ax,0,0,lon,lat,Vorticity_GT_s4[i],"Vorticity",\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
    return ax[0][0]

dr=range(0,365)
fig, ax = plt.subplots(1,1,figsize=(15,15),squeeze=False,\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
ani = animation.FuncAnimation(fig, animate2, frames=dr, blit=False, interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=3, bitrate=5000)
ani.save(workpath+"/animation_vorticity_"+domain+"_1_5.mp4", writer = writer)
plt.close()

# plot mean(|q(t)-q(0)|**2)/var(q(0)
lday=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=151+i),"%Y-%m-%d") for i in range(32) ]
dr=range(151,183)
nmq=[np.mean((Vorticity_GT_s4[x]-Vorticity_GT_s4[dr[0]])**2)/np.var(Vorticity_GT_s4[dr[0]]) for x in dr]
fig, ax = plt.subplots()
plt.plot(lday,nmq)
plt.gcf().autofmt_xdate()
plt.ylabel(r'$\frac{\overline{(q_t-q_0)^2}}{\sigma^2(q_0)}$')
ax.yaxis.label.set_size(16)
plt.savefig(workpath+'/Normalized_mean_vorticity_deviations_MARCH_1_5.png',  bbox_inches='tight')

lday=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=92+i),"%Y-%m-%d") for i in range(32) ]
dr=range(92,124)
nmq=[np.mean((Vorticity_GT_s4[x]-Vorticity_GT_s4[dr[0]])**2)/np.var(Vorticity_GT_s4[dr[0]]) for x in dr]
fig, ax = plt.subplots()
plt.plot(lday,nmq)
plt.gcf().autofmt_xdate()
plt.ylabel(r'$\frac{\overline{(q_t-q_0)^2}}{\sigma^2(q_0)}$')
ax.yaxis.label.set_size(16)
plt.savefig(workpath+'/Normalized_mean_vorticity_deviations_JANUARY_1_5.png', bbox_inches='tight')


