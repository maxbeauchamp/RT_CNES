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

type_obs = sys.argv[1]
workpath = "/home3/scratch/mbeaucha/scores_GENN_"+type_obs
scratchpath = '/home3/scratch/mbeaucha'
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)

# Reload AnDA results (for GT and OI)
file_results_AnDA=scratchpath+'/resAnDA_nadirswot_nadlag_5_'+type_obs+'/saved_path.pickle'
with open(file_results_AnDA, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot = AnDA_ssh_1  
    itrp_dineof_nadirswot = itrp_dineof

# Reload GE-NN results
file_results_1=scratchpath+'/resIA_nadir_nadlag_5_'+type_obs+'/FP_GENN_wwmissing_wocov/saved_path_014_FP_GENN_wwmissing.pickle'
with open(file_results_1, 'rb') as handle:
    itrp_FP_GENN_1 = pickle.load(handle)[2]
file_results_2=scratchpath+'/resIA_nadirswot_nadlag_5_'+type_obs+'/FP_GENN_wwmissing_wocov/saved_path_014_FP_GENN_wwmissing.pickle'
with open(file_results_2, 'rb') as handle:
    itrp_FP_GENN_2 = pickle.load(handle)[2]

lon = np.arange(-65,-55,1/20)
lat = np.arange(30,40,1/20)
indLat  = np.arange(0,200)
indLon  = np.arange(0,200)
lon = lon[indLon]
lat = lat[indLat]
extent_ = [np.min(lon),np.max(lon),np.min(lat),np.max(lat)]

## Init variables for temporal analysis
nrmse_FP_GENN_1=np.zeros(len(itrp_FP_GENN_1))
nrmse_FP_GENN_2=np.zeros(len(itrp_FP_GENN_2))

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
indLon=200
indLat=200

for i in range(0,len(AnDA_ssh_1.GT)):

    day=lday[i]
    print(day)
    # Load data
    gt 				= AnDA_ssh_1.GT[i,:indLon,:indLat]
    obs_nadirswot               = AnDA_ssh_1_nadirswot.Obs[i,:indLon,:indLat]
    # nadirswot
    FP_GENN_1                   = itrp_FP_GENN_1[i,:indLon,:indLat]
    FP_GENN_2                   = itrp_FP_GENN_2[i,:indLon,:indLat]
    ## Compute spatial coverage
    nadswotcov[i]	= len(np.argwhere(np.isfinite(obs_nadirswot.flatten())))/len(obs_nadirswot.flatten())

    ## Compute NRMSE statistics (i.e. RMSE/stdev(gt) )
    nrmse_FP_GENN_1[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_1-np.nanmean(FP_GENN_1)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_2[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_2-np.nanmean(FP_GENN_2)))**2)))/np.nanstd(gt)



N = len(lday)
print(N)
# first axis with nRMSE time

plt.plot(range(N),nrmse_FP_GENN_1,linestyle='solid',color='red',linewidth=1,label='GENN_MM (nadir)')
plt.plot(range(N),nrmse_FP_GENN_2,linestyle='solid',color='blue',linewidth=1,label='GENN_MM (nadirswot)')
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
resfile=workpath+"/TS_AnDA_nRMSE_MM.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()    

