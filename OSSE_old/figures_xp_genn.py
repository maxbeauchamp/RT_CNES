#!/usr/bin/env python

""" display_results.py: script to display maps, taylors diagrams and radially averaged power spectrums from MS-PB-AnDA and MS-VE-DINEOF with three different datasets (nadir / swot / nadirswot) """

__author__ = "Maxime Beauchamp"
__version__ = "0.1"
__date__ = "2019-12-10"
__email__ = "maxime.beauchamp76@gmail.com"

from pb_anda import *
import matplotlib.dates as mdates
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

opt      = sys.argv[1]
type_obs = sys.argv[2]
domain   = sys.argv[3]
workpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/OSSE/"+domain+"/scores_GENN_"+type_obs
scratchpath = '/users/local/m19beauc/4DVARNN-DinAE_xp/OSSE/'+domain
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
#else:
#    shutil.rmtree(workpath)
#    mk_dir_recursive(workpath)    

# Reload AnDA results (for GT and OI)
file_results_AnDA=scratchpath+'/resAnDA_'+opt+'_nadlag_5_'+type_obs+'/saved_path.pickle'
with open(file_results_AnDA, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot = AnDA_ssh_1  
    itrp_dineof_nadirswot = itrp_dineof

# Reload GE-NN results
file_results_1=scratchpath+'/resIA_'+opt+'_nadlag_5_'+type_obs+'/FP_GENN_womissing_wocov/saved_path_019_FP_GENN_womissing.pickle'
with open(file_results_1, 'rb') as handle:
    itrp_FP_GENN_1 = pickle.load(handle)[2]
file_results_2=scratchpath+'/resIA_'+opt+'_nadlag_5_'+type_obs+'/FP_GENN_wmissing_wocov/saved_path_019_FP_GENN_wmissing.pickle'
with open(file_results_2, 'rb') as handle:
    itrp_FP_GENN_2 = pickle.load(handle)[2]
file_results_3=scratchpath+'/resIA_'+opt+'_nadlag_5_'+type_obs+'/FP_GENN_womissing_wOI/saved_path_019_FP_GENN_womissing.pickle'
with open(file_results_3, 'rb') as handle:
    itrp_FP_GENN_3 = pickle.load(handle)[2]
file_results_4=scratchpath+'/resIA_'+opt+'_nadlag_5_'+type_obs+'/FP_GENN_wmissing_wOI/saved_path_019_FP_GENN_wmissing.pickle'
with open(file_results_4, 'rb') as handle:
    itrp_FP_GENN_4 = pickle.load(handle)[2]
file_results_5=scratchpath+'/resIA_'+opt+'_nadlag_5_'+type_obs+'/FP_GENN_wwmissing_wocov/saved_path_014_FP_GENN_wwmissing.pickle'
with open(file_results_5, 'rb') as handle:
    itrp_FP_GENN_5 = pickle.load(handle)[2]
file_results_6=scratchpath+'/resIA_'+opt+'_nadlag_5_'+type_obs+'/FP_GENN_wwmissing_wOI/saved_path_019_FP_GENN_wwmissing.pickle'
with open(file_results_6, 'rb') as handle:
    itrp_FP_GENN_6 = pickle.load(handle)[2]

			#*****************#
			# Display results #
			#*****************#

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


## Init variables for temporal analysis
nrmse_FP_GENN_1=np.zeros(len(itrp_FP_GENN_1))
nrmse_FP_GENN_2=np.zeros(len(itrp_FP_GENN_2))
nrmse_FP_GENN_3=np.zeros(len(itrp_FP_GENN_3))
nrmse_FP_GENN_4=np.zeros(len(itrp_FP_GENN_4))
nrmse_FP_GENN_5=np.zeros(len(itrp_FP_GENN_5))
nrmse_FP_GENN_6=np.zeros(len(itrp_FP_GENN_6))

## Observations spatial coverage
nadswotcov=np.zeros(len(AnDA_ssh_1.GT))

# list of dates
lday1=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=60+i),"%Y-%m-%d") for i in range(20) ]
lday2=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=140+i),"%Y-%m-%d") for i in range(20) ]
lday3=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=220+i),"%Y-%m-%d") for i in range(20) ]
lday4=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=300+i),"%Y-%m-%d") for i in range(20) ]
lday = np.concatenate([lday1,lday2,lday3,lday4])
lday2 = [ datetime.strptime(lday[i],'%Y-%m-%d') for i in range(len(lday)) ] 

## Spatial analysis

for i in range(0,len(AnDA_ssh_1.GT)):

    day=lday[i]
    print(day)
    # Load data
    gt 				= AnDA_ssh_1.GT[i,:indLat,:indLon]
    Grad_gt             	= Gradient(gt,2)
    obs_nadirswot               = AnDA_ssh_1_nadirswot.Obs[i,:indLat,:indLon]
    # nadirswot
    FP_GENN_1                   = itrp_FP_GENN_1[i,:indLat,:indLon]
    Grad_FP_GENN_1              = Gradient(FP_GENN_1,2)
    FP_GENN_2                   = itrp_FP_GENN_2[i,:indLat,:indLon]
    Grad_FP_GENN_2              = Gradient(FP_GENN_2,2)
    FP_GENN_3                   = itrp_FP_GENN_3[i,:indLat,:indLon]
    Grad_FP_GENN_3              = Gradient(FP_GENN_3,2)
    FP_GENN_4                   = itrp_FP_GENN_4[i,:indLat,:indLon]
    Grad_FP_GENN_4              = Gradient(FP_GENN_4,2)
    FP_GENN_5                   = itrp_FP_GENN_5[i,:indLat,:indLon]
    Grad_FP_GENN_5              = Gradient(FP_GENN_5,2)
    FP_GENN_6                   = itrp_FP_GENN_6[i,:indLat,:indLon]
    Grad_FP_GENN_6              = Gradient(FP_GENN_6,2)

    ## Compute spatial coverage
    nadswotcov[i]	= len(np.argwhere(np.isfinite(obs_nadirswot.flatten())))/len(obs_nadirswot.flatten())

    ## Compute NRMSE statistics (i.e. RMSE/stdev(gt) )
    nrmse_FP_GENN_1[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_1-np.nanmean(FP_GENN_1)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_2[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_2-np.nanmean(FP_GENN_2)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_3[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_3-np.nanmean(FP_GENN_3)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_4[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_4-np.nanmean(FP_GENN_4)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_5[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_5-np.nanmean(FP_GENN_5)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_6[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_6-np.nanmean(FP_GENN_6)))**2)))/np.nanstd(gt)

index=list(range(5,16))
index.extend(range(25,36))
index.extend(range(45,56))
index.extend(range(65,76))
tab_scores = np.zeros((6,3))
tab_scores[0,0] = np.nanmean(nrmse_FP_GENN_1[index])
tab_scores[0,1] = np.percentile(nrmse_FP_GENN_1[index],5)
tab_scores[0,2] = np.percentile(nrmse_FP_GENN_1[index],95)
tab_scores[1,0] = np.nanmean(nrmse_FP_GENN_2[index])
tab_scores[1,1] = np.percentile(nrmse_FP_GENN_2[index],5)
tab_scores[1,2] = np.percentile(nrmse_FP_GENN_2[index],95)
tab_scores[2,0] = np.nanmean(nrmse_FP_GENN_3[index])
tab_scores[2,1] = np.percentile(nrmse_FP_GENN_3[index],5)
tab_scores[2,2] = np.percentile(nrmse_FP_GENN_3[index],95)
tab_scores[3,0] = np.nanmean(nrmse_FP_GENN_4[index])
tab_scores[3,1] = np.percentile(nrmse_FP_GENN_4[index],5)
tab_scores[3,2] = np.percentile(nrmse_FP_GENN_4[index],95)
tab_scores[4,0] = np.nanmean(nrmse_FP_GENN_5[index])
tab_scores[4,1] = np.percentile(nrmse_FP_GENN_5[index],5)
tab_scores[4,2] = np.percentile(nrmse_FP_GENN_5[index],95)
tab_scores[5,0] = np.nanmean(nrmse_FP_GENN_6[index])
tab_scores[5,1] = np.percentile(nrmse_FP_GENN_6[index],5)
tab_scores[5,2] = np.percentile(nrmse_FP_GENN_6[index],95)
np.savetxt(fname=workpath+"/tab_scores_"+opt+".txt",X=tab_scores,fmt='%2.2f')


N = len(lday)
print(N)
# first axis with nRMSE time series
plt.plot(range(N),nrmse_FP_GENN_1,linestyle='solid',color='red',linewidth=1,label='GENN_NMNM')
plt.plot(range(N),nrmse_FP_GENN_2,linestyle='dashed',color='red',linewidth=1,label='GENN_MNM')
plt.plot(range(N),nrmse_FP_GENN_5,linestyle='dotted',color='red',linewidth=1,label='GENN_MM')
plt.plot(range(N),nrmse_FP_GENN_3,linestyle='solid',color='blue',linewidth=1,label='GENN_NMNM + OI')
plt.plot(range(N),nrmse_FP_GENN_4,linestyle='dashed',color='blue',linewidth=1,label='GENN_MNM + OI')
plt.plot(range(N),nrmse_FP_GENN_6,linestyle='dotted',color='blue',linewidth=1,label='GENN_MM + OI')
# add vertical bar to divide the 4 periods
plt.axvline(x=19)
plt.axvline(x=39)
plt.axvline(x=59)
# graphical options
plt.ylabel('nRMSE')
plt.xlabel('Time (days)')
plt.xticks([0,20,40,60],\
           [lday[0],lday[20],lday[40],lday[60]],\
           rotation=45, ha='right')
plt.margins(x=0)
plt.grid(True,alpha=.3)
plt.legend(loc='upper left',prop=dict(size='small'),frameon=False,bbox_to_anchor=(0,1.02,1,0.2),ncol=2,mode="expand")
# second axis with spatial coverage
axes2 = plt.twinx()
width=0.75
p1 = axes2.bar(range(N), nadswotcov, width,color='r',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_AnDA_nRMSE_"+opt+".png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()         	# close the figure

