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

if wregul==True:
    workpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/OSE_keras/OSE/"+domain+"/figs_lag"+lag+"_wregul"
else:
    workpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/OSE_keras/OSE/"+domain+"/figs_lag"+lag
scratchpath = '/users/local/m19beauc/4DVARNN-DinAE_xp/OSE_keras/OSE/'+domain
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
#else:
#    shutil.rmtree(workpath)
#    mk_dir_recursive(workpath)    

# Reload saved GENN result
file_results_GENN_wotrain=scratchpath+'/resIA_nadir_nadlag_'+lag+'_obs/FP_GENN_wwmissing_wOI_wotrain/saved_path_FP_GENN_wwmissing.pickle'
with open(file_results_GENN_wotrain, 'rb') as handle:
    GT, itrp_OI, Obs_nadir, itrp_FP_GENN_wotrain = pickle.load(handle)[:4]
if wregul==True:
    file_results_GENN_wtrain=scratchpath+'/resIA_nadir_nadlag_'+lag+'_obs/FP_GENN_wwmissing_wOI_wtrain_wregul/saved_path_FP_GENN_wwmissing.pickle'
else:
    file_results_GENN_wtrain=scratchpath+'/resIA_nadir_nadlag_'+lag+'_obs/FP_GENN_wwmissing_wOI_wtrain/saved_path_FP_GENN_wwmissing.pickle'
with open(file_results_GENN_wtrain, 'rb') as handle:
    itrp_FP_GENN_wtrain = pickle.load(handle)[3]

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

##*** INIT SSH STATISTICS ***##
## Init variables for temporal analysis (R scores)
nrmse_OI=np.zeros(len(itrp_OI))
nrmse_FP_GENN_wotrain=np.zeros(len(itrp_FP_GENN_wotrain))
nrmse_FP_GENN_wtrain=np.zeros(len(itrp_FP_GENN_wtrain))

## Observations spatial coverage
nadcov=np.zeros(len(data_valid))
nadcov_gt=np.zeros(len(data_valid))

# list of dates
lday=[ datetime.strftime(datetime.strptime("2017-01-01",'%Y-%m-%d')\
                          + timedelta(days=i),"%Y-%m-%d") for i in range(365) ]
lday2 = [ datetime.strptime(lday[i],'%Y-%m-%d') for i in range(len(lday)) ] 

## Spatial analysis

for i in range(0,len(GT)):

    day=lday[i]
    print(day)

    # Load data
    gt                          = data_valid[i,5:indLat,5:indLon]
    OI                          = itrp_OI[i,5:indLat,5:indLon]
    obs 			= Obs_nadir[i,5:indLat,5:indLon]
    FP_GENN_wotrain             = itrp_FP_GENN_wotrain[i,5:indLat,5:indLon]
    FP_GENN_wotrain             = median_filter(FP_GENN_wotrain, size=5)
    itrp_FP_GENN_wotrain[i,5:indLat,5:indLon] = FP_GENN_wotrain
    FP_GENN_wtrain              = itrp_FP_GENN_wtrain[i,5:indLat,5:indLon]
    FP_GENN_wtrain              = median_filter(FP_GENN_wtrain, size=5)
    itrp_FP_GENN_wtrain[i,5:indLat,5:indLon] = FP_GENN_wtrain
    Grad_OI                     = Gradient(OI,2)
    Grad_FP_GENN_wotrain        = Gradient(FP_GENN_wotrain,2)
    Grad_FP_GENN_wtrain         = Gradient(FP_GENN_wtrain,2)
    mask                        = np.ones(gt.shape)
    mask[:,:]                   = np.where(np.isnan(gt),np.nan,1)
    ## Compute spatial coverage
    nadcov[i]	          	= len(np.argwhere(np.isfinite(obs.flatten())))/len(obs.flatten())
    nadcov_gt[i]                = len(np.argwhere(np.isfinite(gt.flatten())))/len(gt.flatten())

    ## Compute NRMSE statistics (i.e. RMSE/stdev(gt))
    nrmse_OI[i]	= (np.sqrt(np.nanmean(((mask*gt-np.nanmean(mask*gt))-(mask*OI-np.nanmean(mask*OI)))**2)))/np.nanstd(mask*gt)
    nrmse_FP_GENN_wotrain[i]    = (np.sqrt(np.nanmean(((mask*gt-np.nanmean(mask*gt))-(mask*FP_GENN_wotrain-np.nanmean(mask*FP_GENN_wotrain)))**2)))/np.nanstd(mask*gt)
    nrmse_FP_GENN_wtrain[i]     = (np.sqrt(np.nanmean(((mask*gt-np.nanmean(mask*gt))-(mask*FP_GENN_wtrain-np.nanmean(mask*FP_GENN_wtrain)))**2)))/np.nanstd(mask*gt)

    if day=='2017-08-04':
        # Display maps
        fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(22,10),
                         subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
        vmin = np.nanmin(obs) ; vmax = np.nanmax(obs)
        cmap="coolwarm"
        plot2(ax1,lon[5:],lat[5:],OI,'OI',\
             extent=extent2,cmap=cmap,vmin=vmin,vmax=vmax)
        plot2(ax2,lon[5:],lat[5:],FP_GENN_wotrain,'FP-GENN (NATL60 2013 model)',\
             extent=extent2,cmap=cmap,vmin=vmin,vmax=vmax)
        plot2(ax3,lon[5:],lat[5:],FP_GENN_wtrain,'FP-GENN (nadir 2017 model)',\
             extent=extent2,cmap=cmap,vmin=vmin,vmax=vmax)
        plt.subplots_adjust(hspace=0.3,wspace=0.5)
        resfile=workpath+"/SSH_maps_"+day+".png"
        plt.savefig(resfile)	# save the figure
        plt.close()			# close the figure

        # Display gradients
        fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(22,10),
                         subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
        vmin = np.nanmin(Grad_FP_GENN_wtrain) ; vmax = np.nanmax(Grad_FP_GENN_wtrain)
        vmin = 0. ; vmax = 0.1 
        cmap="viridis"
        plot2(ax1,lon[5:],lat[5:],Grad_OI,r"$\nabla_{OI}$",\
             extent=extent2,cmap=cmap,vmin=vmin,vmax=vmax)
        plot2(ax2,lon[5:],lat[5:],Grad_FP_GENN_wotrain,r"$\nabla_{FP-GENN}$ (NATL60 2013 model)",\
             extent=extent2,cmap=cmap,vmin=vmin,vmax=vmax)
        plot2(ax3,lon[5:],lat[5:],Grad_FP_GENN_wtrain,r"$\nabla_{FP-GENN}$ (nadir 2017 model)",\
             extent=extent2,cmap=cmap,vmin=vmin,vmax=vmax)
        plt.subplots_adjust(hspace=0.3,wspace=0.5)
        resfile=workpath+"/Grad_SSH_maps_"+day+".png"
        plt.savefig(resfile)        # save the figure
        plt.close()                 # close the figure'''

## animations (videos) ##
if ( lag == "0" and wregul == False ):

    def animate(i):
        print(i)
        # Load data
        OI                          = itrp_OI[i,5:indLat,5:indLon]
        Grad_OI                     = Gradient(OI,2)
        FP_GENN_wtrain              = itrp_FP_GENN_wtrain[i,5:indLat,5:indLon]
        Grad_FP_GENN                = Gradient(FP_GENN_wtrain,2)
        for ivar in range(len(var)):
            plot(ax,0,ivar,lon[5:],lat[5:],eval(var[ivar]),title[ivar],\
                extent=extent2,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
        return ax[0][0], ax[0][1]

    var=['Grad_OI','Grad_FP_GENN']
    title=[r"$\nabla_{OI}$",r"$\nabla_{FP-GENN}$"]
    vmax=0.1
    vmin=0.
    cmap="viridis"
    fig, ax = plt.subplots(1,2,figsize=(22,10),squeeze=False,\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
    plt.subplots_adjust(hspace=0.5)
    ani = animation.FuncAnimation(fig, animate, blit=False, frames=np.arange(1,365), interval=1000, repeat=False)
    writer = animation.FFMpegWriter(fps=10, bitrate=5000)
    ani.save(workpath+"/animation_grads.mp4", writer = writer)
    plt.close()

## SSH score tables 
tab_scores = np.zeros((3,3))
tab_scores[0,0] = np.nanmean(nrmse_OI)
tab_scores[0,1] = np.nanpercentile(nrmse_OI,5)
tab_scores[0,2] = np.nanpercentile(nrmse_OI,95)
tab_scores[1,0] = np.nanmean(nrmse_FP_GENN_wotrain)
tab_scores[1,1] = np.nanpercentile(nrmse_FP_GENN_wotrain,5)
tab_scores[1,2] = np.nanpercentile(nrmse_FP_GENN_wotrain,95)
tab_scores[2,0] = np.nanmean(nrmse_FP_GENN_wtrain)
tab_scores[2,1] = np.nanpercentile(nrmse_FP_GENN_wtrain,5)
tab_scores[2,2] = np.nanpercentile(nrmse_FP_GENN_wtrain,95)
np.savetxt(fname=workpath+"/tab_scores_OSE_SSH.txt",X=tab_scores,fmt='%2.2f')

# Taylor diagram
# apply HPF to visualize Taylor diagrams only for small scales
HR = itrp_OI[:,:indLat,:indLon]
lr = np.copy(HR).reshape(HR.shape[0],-1)
tmp = lr[0,:]
sea_v2 = np.where(~np.isnan(tmp))[0]
lr_no_land = lr[:,sea_v2]
pca = PCA(n_components=.75)
score_global = pca.fit_transform(lr_no_land)
coeff_global = pca.components_.T
mu_global = pca.mean_
DataReconstructed_global = np.dot(score_global, coeff_global.T) +mu_global
lr[:,sea_v2] = DataReconstructed_global
lr = lr.reshape(HR.shape)[:,20:-20,20:-20].flatten()
resfile=workpath+"/Taylor_diagram.png"
label=['GT','OI','FP-GENN (NATL60 2013 model)', 'FP-GENN (nadir 2017 model)']
mask                      = np.ones(data_valid[:,20:-20,20:-20].flatten().shape)
mask[:]                   = np.nan
mask[:]                   = np.where(np.isnan(data_valid[:,20:-20,20:-20].flatten()),np.nan,1)
print(data_valid[:,20:-20,20:-20].shape)
print(itrp_OI[:,20:-20,20:-20].shape)
series={'gt':mask*data_valid[:,20:-20,20:-20].flatten()-lr,
        'OI':mask*itrp_OI[:,20:-20,20:-20].flatten()-lr,
        'FP_GENN_wotrain':mask*itrp_FP_GENN_wotrain[:,20:-20,20:-20].flatten()-lr,
        'FP_GENN_wtrain':mask*itrp_FP_GENN_wtrain[:,20:-20,20:-20].flatten()-lr}
Taylor_diag(series,label,\
            styles=['k','p','p','p',],\
            colors=['k','red','seagreen','steelblue'])
plt.savefig(resfile)
plt.close()

## Plot time series
N = len(lday)
ymax_ = np.ceil(np.nanmax([nrmse_OI, nrmse_FP_GENN_wotrain, nrmse_FP_GENN_wtrain]))
ymax_=.5
# first axis with nRMSE time series
smask  = np.where(np.isfinite(nrmse_OI))[0]
smask2 = np.where(nrmse_OI<=0.4)[0]
smask  = np.intersect1d(smask,smask2)
plt.plot(np.arange(0,365)[smask],nrmse_OI[smask],linestyle='solid',color='red',linewidth=2,label='OI')
plt.plot(np.arange(0,365)[smask],nrmse_FP_GENN_wotrain[smask],linestyle='solid',color='seagreen',linewidth=1,markerSize=2,label='FP-GENN (NATL60 2013 model)')
plt.plot(np.arange(0,365)[smask],nrmse_FP_GENN_wtrain[smask],linestyle='solid',color='steelblue',linewidth=1,markerSize=2,label='FP-GENN (nadir 2017 model)')
# graphical options
plt.ylim(0,ymax_)
plt.ylabel('nRMSE')
plt.xlabel('Time (days)')
plt.xticks([0,31,59,90,120,151,181,212,243,273,304,334],\
           [lday[0],lday[31],lday[59],lday[90],lday[120],lday[151],lday[181],\
           lday[212],lday[243],lday[273],lday[304],lday[334]],\
           rotation=45, ha='right')
plt.margins(x=0)
plt.grid(True,alpha=.3)
plt.legend(loc='upper left',prop=dict(size='small'),frameon=False,bbox_to_anchor=(0,1.02,1,0.2),ncol=3,mode="expand")
# second axis with spatial coverage
axes2 = plt.twinx()
width=0.75
p1 = axes2.bar(range(0,365), nadcov, width,color='r',alpha=0.25)
p2 = axes2.bar(range(0,365), nadcov_gt, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_nRMSE.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()             # close the figure

## Plot averaged RAPS
resfile=workpath+"/results_avg_RAPS.png"
resssh=4
f0, Pf_OI        = avg_raPsd2dv1(itrp_OI,resssh,True)
f1, Pf_FP_GENN   = avg_raPsd2dv1(itrp_FP_GENN_wtrain,resssh,True)
wf0              = 1/f0
wf1              = 1/f1
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf0[1:],Pf_OI[1:],label='OI',color='red',linewidth=2)
ax.plot(wf0[1:],Pf_FP_GENN[1:],label='FP-GENN (nadir 2017 model)',color='steelblue',linewidth=2)
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

