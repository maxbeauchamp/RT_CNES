from tools_OSSE import *
import xarray as xr
import xrft
import logging
import xrft
from dask.diagnostics import ProgressBar
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

## 1) Function for plotting nRMSE
def plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False):

    N   = len(lday)
    GT  = list_data[0]
    ## Compute spatial coverage  
    if any("Obs (nadir+swot)"==s for s in labels_data):
        spatial_coverage_nadirswot = []
        id_obs = np.where(labels_data=="Obs (nadir+swot)")[0][0]
        OBS = list_data[id_obs]
        for j in range(len(lday)):
            spatial_coverage_nadirswot.append(100*len(np.argwhere(np.isfinite(OBS[j].flatten())))/len(OBS[j].flatten()))
    if any("Obs (nadir)"==s for s in labels_data):
        spatial_coverage_nadir = []
        id_obs = np.where(labels_data=="Obs (nadir)")[0][0]
        OBS = list_data[id_obs]
        for j in range(len(lday)):
            spatial_coverage_nadir.append(100*len(np.argwhere(np.isfinite(OBS[j].flatten())))/len(OBS[j].flatten()))
        
    # Compute nRMSE time series
    nRMSE = []
    for i in range(len(labels_data[2:])):
        nRMSE_=[]
        meth_i=list_data[i+2]
        for j in range(len(lday)):
            if gradient == False:
                nRMSE_.append((np.sqrt(np.nanmean(((GT[j]-np.nanmean(GT[j]))-(meth_i[j]-np.nanmean(meth_i[j])))**2)))/np.nanstd(GT[j]))
            else:
                nRMSE_.append((np.sqrt(np.nanmean(((Gradient(GT[j],2)-np.nanmean(Gradient(GT[j],2)))-(Gradient(meth_i[j],2)-np.nanmean(Gradient(meth_i[j],2))))**2)))/np.nanstd(Gradient(GT[j],2)))
        nRMSE.append(nRMSE_)
  
    # plot nRMSE time series
    for i in range(len(labels_data[2:])):
        if gradient == False:
            plt.plot(range(N),nRMSE[i],linestyle=lstyle[i+2],color=colors[i+2],linewidth=lwidth[i+2],label=labels_data[i+2])
        else:
            plt.plot(range(N),nRMSE[i],linestyle=lstyle[i+2],color=colors[i+2],linewidth=lwidth[i+2],label=r"$\nabla_{"+str.split(labels_data[i+2])[0]+"}$ "+str.split(labels_data[i+2])[1])


    # add vertical bar to divide the 4 periods
    plt.axvline(x=19)
    plt.axvline(x=39)
    plt.axvline(x=59)
    # graphical options
    plt.ylim(0,np.ceil(np.max(nRMSE)*100)/100)
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
    if ( (any("Obs (nadir+swot)"==s for s in labels_data)) and any("Obs (nadir)"==s for s in labels_data) ):
        p1 = axes2.bar(range(N), spatial_coverage_nadir, width,color='r',alpha=0.25)
        p2 = axes2.bar(range(N), spatial_coverage_nadirswot-spatial_coverage_nadir, width,bottom=spatial_coverage_nadir,color='g',alpha=0.25)
    elif any("Obs (nadir+swot)"==s for s in labels_data):
        p1 = axes2.bar(range(N), spatial_coverage_nadirswot, width,color='g',alpha=0.25)
    else:
        p1 = axes2.bar(range(N), spatial_coverage_nadir, width,color='r',alpha=0.25)
    axes2.set_ylim(0, 100)
    axes2.set_ylabel('Spatial Coverage (%)')
    axes2.margins(x=0)
    plt.savefig(resfile,bbox_inches="tight")    # save the figure
    plt.close()                                 # close the figure

## 2) Function for plotting Signal-to-Noise ratio
def plot_SNR(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resssh,resfile):

    # select only 10-day windows
    index=list(range(5,16))
    index.extend(range(25,36))
    index.extend(range(45,56))
    index.extend(range(65,76))

    GT  = list_data[0][index]

    # Compute Signal-to-Noise ratio
    SNR = []
    for i in range(len(labels_data[2:])):
        print(labels_data[i+2])
        f, Pf  = avg_err_raPsd2dv1(list_data[i+2][index],GT,resssh,True)
        wf     = 1./f
        SNR.append([wf, Pf])

    # plot Signal-to-Noise ratio
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(labels_data[2:])):
        wf, Pf = SNR[i]
        ax.plot(wf,Pf,linestyle=lstyle[i+2],color=colors[i+2],linewidth=lwidth[i+2],label=labels_data[i+2])
    #plt.axhline(y=0.5, color='r', linestyle='-')
    ax.set_xlabel("Wavenumber", fontweight='bold')
    ax.set_ylabel("Signal-to-noise ratio", fontweight='bold')
    ax.set_xscale('log') ; ax.set_yscale('log')
    plt.legend(loc='best',prop=dict(size='small'),frameon=False)
    plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
    ax.invert_xaxis()
    plt.grid(which='both', linestyle='--')
    plt.savefig(resfile) # save the figure
    plt.close()          # close the figure

## 3) Function for plotting Taylor diagrams
def Taylor_diagram(list_data,labels_data,colors,lstyle,resfile):

    # select only 10-day windows
    index=list(range(5,16))
    index.extend(range(25,36))
    index.extend(range(45,56))
    index.extend(range(65,76))

    # apply High-Pass Filter to visualize Taylor diagrams only for small scales
    HR = list_data[2][index]
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
    lr = lr.reshape(HR.shape).flatten()

    # create a dictionnary with data
    series={}
    for i in np.delete(np.arange(len(list_data)),1):
        series[labels_data[i]] = list_data[i][index].flatten()-lr
    Taylor_diag(series,np.delete(labels_data,1),\
            styles=np.delete(lstyle,1),\
            colors=np.delete(colors,1))
    plt.savefig(resfile)
    plt.close()


## 4) Display individual maps (SSH & Gradients)
def plot_maps(list_data,list_suffix,labels_data,list_day,date,extent,lon,lat,workpath):
    # vmin/vmax based on GT
    iday = np.where(list_day==date)[0][0]
    vmax = np.nanmax(np.abs(list_data[0]))
    vmin = -1.*vmax
    grad_vmax = np.nanmax(np.abs(Gradient(list_data[0][iday],2)))
    grad_vmin = 0
    for i in range(len(list_data)):
        resfile = workpath+"/results_"+list_suffix[i]+'_'+date+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),squeeze=False,
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))   
        plot(ax,0,0,lon,lat,list_data[i][iday],labels_data[i],\
             extent=extent,cmap="coolwarm",vmin=vmin,vmax=vmax)
        plt.savefig(resfile)       # save the figure
        plt.close()                # close the figure

        resfile = workpath+"/results_Grad_"+list_suffix[i]+'_'+date+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),squeeze=False,
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))          
        if len(str.split(labels_data[i]))>1:
            plot(ax,0,0,lon,lat,Gradient(list_data[i][iday],2),r"$\nabla_{"+str.split(labels_data[i])[0]+"}$ "+str.split(labels_data[i])[1],\
                extent=extent,cmap="viridis",vmin=grad_vmin,vmax=grad_vmax)
        else:
            plot(ax,0,0,lon,lat,Gradient(list_data[i][iday],2),r"$\nabla_{"+labels_data[i]+"}$",\
                extent=extent,cmap="viridis",vmin=grad_vmin,vmax=grad_vmax)
        plt.savefig(resfile)       # save the figure
        plt.close()                # close the figure

## 5) Compute R/I/AE-scores
def RIAE_scores(list_data,labels_data,resfile,gradient=False):
    
    GT  = list_data[0]

    # apply High-Pass Filter to visualize Taylor diagrams only for small scales
    HR = list_data[2]
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
    lr = lr.reshape(HR.shape).flatten()

    def Iscore(mask1,gt,itrp):
        return 100*(1-np.nanmean(((mask1*gt-np.nanmean(mask1*gt))-(mask1*itrp-np.nanmean(mask1*itrp)))**2)/np.nanvar(mask1*gt))
    def Rscore(mask1,gt,itrp):
        return 100*(1-np.nanmean(((mask1*gt-np.nanmean(mask1*gt))-(mask1*itrp-np.nanmean(mask1*itrp)))**2)/np.nanvar(mask1*gt))
    def AEscore(gt,itrp):
        return 100*(1-np.nanmean(((gt-np.nanmean(gt))-(itrp-np.nanmean(itrp)))**2)/np.nanvar(gt))

    ## Compute mask nadir / swot  
    if any("Obs (nadir+swot)"==s for s in labels_data):
        spatial_coverage_nadir = []
        id_obs = np.where(labels_data=="Obs (nadir+swot)")[0][0]
        OBS = list_data[id_obs]
        mask1_nadirswot             = np.where(np.isnan(OBS.flatten()),np.nan,1)
        mask2_nadirswot             = np.where(np.isnan(OBS.flatten()),1,np.nan)
    if any("Obs (nadir)"==s for s in labels_data):
        spatial_coverage_nadir = []
        id_obs = np.where(labels_data=="Obs (nadir)")[0][0]
        OBS = list_data[id_obs]
        mask1_nadir                 = np.where(np.isnan(OBS).flatten(),np.nan,1)
        mask2_nadir                 = np.where(np.isnan(OBS.flatten()),1,np.nan)

    # Compute R/I/AE scores
    tab_scores = np.zeros((len(labels_data[2:]),3))
    for i in range(len(labels_data[2:])):
        meth_i=list_data[i+2]
        print(labels_data[i+2])
        if ("nadir+swot" in labels_data[i+2]):
            tab_scores[i,0] = Rscore(mask1_nadirswot,GT.flatten()-lr,meth_i.flatten()-lr)
            tab_scores[i,1] = Iscore(mask2_nadirswot,GT.flatten()-lr,meth_i.flatten()-lr)
        else:
            tab_scores[i,0] = Rscore(mask1_nadir,GT.flatten()-lr,meth_i.flatten()-lr)
            tab_scores[i,1] = Rscore(mask2_nadir,GT.flatten()-lr,meth_i.flatten()-lr)
        if ("rec" in labels_data[i+2]):
            tab_scores[i,2] = AEscore(GT.flatten()-lr,meth_i.flatten()-lr)
        else:
            tab_scores[i,2] = np.nan
    np.savetxt(fname=resfile,X=tab_scores,fmt='%2.2f')

# 5) export NetCDF
def export_NetCDF(list_data,labels_data,list_day,lon,lat,resfile):

    GT  = list_data[0]
    dt64 = [ np.datetime64(datetime.strptime(day,'%Y-%m-%d')) for day in list_day ]
    time_u = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
    mesh_lat, mesh_lon = np.meshgrid(range(len(lat)),range(len(lon)))
    data = xr.Dataset(\
                    data_vars={'longitude': (('lat','lon'),mesh_lon),\
                                'latitude' : (('lat','lon'),mesh_lat),\
                                'Time'     : (('time'),time_u),\
                                'GT'       : (('time','lat','lon'),GT)},\
                    coords={'lon': lon,\
                            'lat': lat,\
                            #'time': range(0,len(time_u))})
                            'time': time_u})
    # add variables
    for i in range(len(labels_data[1:])):
        data[labels_data[i+1]] = (('time','lat','lon'), list_data[i+1])
    # write to file
    data.to_netcdf(resfile)

## 6) New PSD scores
def PSD(ds,time,method):

    # Compute error = SSH_reconstruction - SSH_true
    err = (ds['GT'] - ds[method])
    # rechunk
    err = err.chunk({"lat":1, 'time': err['time'].size, 'lon': err['lon'].size})
    err['time'] = time

    # Rechunk SSH_true
    signal = ds['GT'].chunk({"lat":1, 'time': ds['time'].size, 'lon': ds['lon'].size})
    # make time vector in days units
    signal['time'] = time
    
    # Compute PSD_err and PSD_signal
    psd_err = xrft.power_spectrum(err, dim=['time', 'lon'], detrend='linear', window=True).compute()
    psd_signal = xrft.power_spectrum(signal, dim=['time', 'lon'], detrend='linear', window=True).compute()
        
    # Averaged over latitude
    mean_psd_signal = psd_signal.mean(dim='lat').where((psd_signal.freq_lon > 0.) & (psd_signal.freq_time > 0), drop=True)
    mean_psd_err = psd_err.mean(dim='lat').where((psd_err.freq_lon > 0.) & (psd_err.freq_time > 0), drop=True)
        
    # return PSD-based score
    psd_based_score = (1.0 - mean_psd_err/mean_psd_signal)
    
    # Find the key metrics: shortest temporal & spatial scales resolved based on the 0.5 contour criterion of the PSD_score
    level = [0.5]
    cs = plt.contour(1./psd_based_score.freq_lon.values,1./psd_based_score.freq_time.values, psd_based_score, level)
    x05, y05 = cs.collections[0].get_paths()[0].vertices.T

    shortest_spatial_wavelength_resolved = np.min(x05)
    shortest_temporal_wavelength_resolved = np.min(y05)
    
    return (1.0 - mean_psd_err/mean_psd_signal), np.round(shortest_spatial_wavelength_resolved, 2), np.round(shortest_temporal_wavelength_resolved, 2)

def plot_psd_score(ds_psd,resfile):

    try:
        nb_experiment = len(ds_psd.experiment)
    except:
        nb_experiment = 1

    fig, ax0 =  plt.subplots(1, nb_experiment, sharey=True, figsize=(5*nb_experiment, 5))
    #plt.subplots_adjust(right=0.1, left=0.09)
    for exp in range(nb_experiment):
        try:
            ctitle = ds_psd.experiment.values[exp]
        except:
            ctitle = ''

        if nb_experiment > 1:
            ax = ax0[exp]
            data = (ds_psd.isel(experiment=exp).values)
        else:
            ax = ax0
            data = (ds_psd.values)
        ax.invert_yaxis()
        ax.invert_xaxis()
        c1 = ax.contourf(1./(ds_psd.freq_lon), 1./ds_psd.freq_time, data,
                          levels=np.arange(0,1.1, 0.1), cmap='RdYlGn', extend='both')
        ax.set_xlabel('spatial wavelength (degree_lon)', fontweight='bold', fontsize=18)
        ax0[0].set_ylabel('temporal wavelength (days)', fontweight='bold', fontsize=18)
        #plt.xscale('log')
        ax.set_yscale('log')
        ax.grid(linestyle='--', lw=1, color='w')
        ax.tick_params(axis='both', labelsize=18)
        ax.set_title(f'PSD-based score ({ctitle})', fontweight='bold', fontsize=18)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(ScalarFormatter())
        c2 = ax.contour(1./(ds_psd.freq_lon), 1./ds_psd.freq_time, data, levels=[0.5], linewidths=2, colors='k')

        cbar = fig.colorbar(c1, ax=ax, pad=0.01)
        cbar.add_lines(c2)

    bbox_props = dict(boxstyle="round,pad=0.5", fc="w", ec="k", lw=2)
    ax0[-1].annotate('Resolved scales',
                    xy=(1.2, 0.8),
                    xycoords='axes fraction',
                    xytext=(1.2, 0.55),
                    bbox=bbox_props,
                    arrowprops=
                        dict(facecolor='black', shrink=0.05),
                        horizontalalignment='left',
                        verticalalignment='center')

    ax0[-1].annotate('UN-resolved scales',
                    xy=(1.2, 0.2),
                    xycoords='axes fraction',
                    xytext=(1.2, 0.45),
                    bbox=bbox_props,
                    arrowprops=
                    dict(facecolor='black', shrink=0.05),
                        horizontalalignment='left',
                        verticalalignment='center')

    plt.savefig(resfile)

def plot_psd(ncdf_file,labels_data,list_day,resfile,periods=[[0,20],[20,40],[40,60],[60,80]]):

    ds = xr.open_dataset(ncdf_file)
    dt64 = [ np.datetime64(datetime.strptime(day,'%Y-%m-%d')) for day in list_day ]
    time_u = (dt64 - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 'D')

    # loop over the 4 periods
    list_PSD_all_periods = []
    for i in range(len(periods)):
        list_PSD = []
        ds2  = ds.isel(time=slice(periods[i][0], periods[i][1]))
        time = time_u[periods[i][0]:periods[i][1]]-time_u[periods[i][0]]
        # loop over the methods
        for j in range(len(labels_data[2:])):
            list_PSD.append(PSD(ds2,time,labels_data[j+2])[0])

        ds_psd = xr.concat(list_PSD, dim='experiment')
        ds_psd['experiment'] = labels_data[2:]
        plot_psd_score(ds_psd,resfile+"_period"+str(i)+".png")
  
        list_PSD_all_periods.append(list_PSD)

    # average
    list_PSD = []
    for j in range(len(labels_data[2:])):

        ds_psd = xr.concat([list_PSD_all_periods[i][j] for i in range(len(periods))],dim='period')
        mean   = ds_psd.mean('period')
        print(mean)
        list_PSD.append(mean)

    ds_psd = xr.concat(list_PSD, dim='experiment')
    ds_psd['experiment'] = labels_data[2:]
    plot_psd_score(ds_psd,resfile+"_mean.png")




