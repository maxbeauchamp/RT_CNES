#!/usr/bin/env python

""" AnDA_stat_functions.py: Collection of statistical functions used in AnDA. """

__author__ = "Maxime Beauchamp"
__version__ = "0.1"
__date__ = "2020-12-10"
__email__ = "maxime.beauchamp@imt-atantique.fr"

import os
from os.path import join as join_paths
import sys
import pickle
import numpy as np
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import mpl_toolkits.axisartist.grid_finder as GF
import mpl_toolkits.axisartist.floating_axes as FA
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import einops
import matplotlib.dates as mdates
from datetime import date, datetime, timedelta
from cartopy import crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cv2
np.random.seed(1)

def Gradient(img, order):
    """ calcuate x, y gradient and magnitude """ 
    sobelx = cv2.Sobel(img,cv2.CV_64F,1,0,ksize=3)
    sobelx = sobelx/8.0
    sobely = cv2.Sobel(img,cv2.CV_64F,0,1,ksize=3)
    sobely = sobely/8.0
    sobel_norm = np.sqrt(sobelx*sobelx+sobely*sobely)
    if (order==0):
        return sobelx
    elif (order==1):
        return sobely
    else:
        return sobel_norm

def Imputing_NaN(data, invalid=None):
    """
    Replace the value of invalid 'data' cells (indicated by 'invalid') 
    by the value of the nearest valid data cell
    """
    import scipy.ndimage as nd
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]

def AnDA_RMSE(a,b):
    """ Compute the Root Mean Square Error between 2 n-dimensional vectors. """
    if (a.ndim==1):
        a = a[np.newaxis]
    if (a.ndim>2):
        a = a.reshape(a.shape[0],-1)
    if (b.ndim==1):
        b = b[np.newaxis]    
    if (b.ndim>2):
        b = b.reshape(b.shape[0],-1)
    return np.sqrt(np.nanmean((a-b)**2,1))

    
def AnDA_CRMSE(a,b):
    """ Compute the Root Mean Square Error between 2 n-dimensional vectors. """
    if (a.ndim==1):
        a = a[np.newaxis]
    if (a.ndim>2):
        a = a.reshape(a.shape[0],-1)
    if (b.ndim==1):
        b = b[np.newaxis]
    if (b.ndim>2):
        b = b.reshape(b.shape[0],-1)
    abar = np.nanmean(a)
    bbar = np.nanmean(b)
    return np.sqrt(np.nanmean(((a-abar)-(b-bbar))**2,1))
    
def AnDA_correlate(a,b):
    """ Compute the Correlation between 2 n-dimensional vectors. """
    if (a.ndim==1):
        a = a[np.newaxis]
    if (a.ndim>2):
        a = a.reshape(a.shape[0],-1)
    if (b.ndim==1):
        b = b[np.newaxis] 
    if (b.ndim>2):
        b = b.reshape(b.shape[0],-1)
    a = a - np.nanmean(a,1)[np.newaxis].T
    b = b - np.nanmean(b,1)[np.newaxis].T
    r = np.nansum((a*b),1) / np.sqrt(np.nansum((a*a),1) * np.nansum((b*b),1))
    return r   

def AnDA_stdev(a):
    """ Compute the Correlation between 2 n-dimensional vectors. """
    if (a.ndim==1):
        a = a[np.newaxis]
    if (a.ndim>2):
        a = a.reshape(a.shape[0],-1)
    abar = np.nanmean(a)
    return np.sqrt(np.nanmean((a-abar)**2,1))

def AnDA_Taylor_stats(a,b):
    """ Compute the Taylor Statistics between 2 n-dimensional vectors. """
    abar = np.nanmean(a)
    bbar = np.nanmean(b)
    crmsd = np.sqrt(np.nanmean(((a-abar)-(b-bbar))**2))
    sd    = np.sqrt(np.nanmean((a-abar)**2))
    a = a - np.nanmean(a)
    b = b - np.nanmean(b)
    corr = np.nansum((a*b)) / np.sqrt(np.nansum((a*a)) * np.nansum((b*b)))
    return [crmsd, corr, sd]

def Taylor_diag(series,names,styles,colors):
    """ Taylor Diagram : obs is reference data sample
        in a full diagram (0 --> npi)
        --------------------------------------------------------------------------
        Input: series     - dict with all time series (lists) to analyze  
               series[0]  - is the observation, the reference by default.
    """
    from matplotlib.projections import PolarAxes
    taylor_stats = np.array([AnDA_Taylor_stats(series[i],series[list(series.keys())[0]]) for i in series.keys()])
    crmsd, corr, std = taylor_stats.T
    ref = std[0]
    rlocs = np.concatenate((np.arange(0,10,0.25),[0.95,0.99]))
    str_rlocs = np.concatenate((np.arange(0,10,0.25),[0.95,0.99]))
    tlocs = np.arccos(rlocs)        # Conversion to polar angles
    gl1 = GF.FixedLocator(tlocs)    # Positions
    tf1 = GF.DictFormatter(dict(zip(tlocs, map(str,rlocs))))   
    str_locs2 = np.arange(0,11,0.5)
    tlocs2 =  np.arange(0,11,0.5)      # Conversion to polar angles  
    g22 = GF.FixedLocator(tlocs2)  
    tf2 = GF.DictFormatter(dict(zip(tlocs2, map(str,str_locs2))))
    tr = PolarAxes.PolarTransform()
    smin = 0
    smax = (120/100)*np.max(std)
    ghelper = FA.GridHelperCurveLinear(tr,
                                           extremes=(0,np.pi/2, # 1st quadrant
                                                     smin,smax),
                                           grid_locator1=gl1,
                                           #grid_locator2=g11,
                                           tick_formatter1=tf1,
                                           tick_formatter2=tf2,
                                           )
    fig = plt.figure(figsize=(10,5), dpi=100)
    ax = FA.FloatingSubplot(fig, 111, grid_helper=ghelper)
    fig.add_subplot(ax)
    ax.axis["top"].set_axis_direction("bottom") 
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("Correlation Coefficient")
    ax.axis["left"].set_axis_direction("bottom") 
    ax.axis["left"].label.set_text("Standard Deviation")
    ax.axis["right"].set_axis_direction("top") 
    ax.axis["right"].toggle(ticklabels=True, label=True)
    ax.axis["right"].set_visible(True)
    ax.axis["right"].major_ticklabels.set_axis_direction("bottom")
    ax.axis["right"].label.set_text("Standard Deviation")
    ax.axis["bottom"].set_visible(False) 
    ax.grid(True)
    ax = ax.get_aux_axes(tr)
    t = np.linspace(0, np.pi/2)
    r = np.zeros_like(t) + ref
    ax.plot(t,r, 'k--', label='_')
    rs,ts = np.meshgrid(np.linspace(smin,smax),
                            np.linspace(0,np.pi/2))
    rms = np.sqrt(ref**2 + rs**2 - 2*ref*rs*np.cos(ts))
    CS =ax.contour(ts, rs,rms,cmap=cm.bone)
    plt.clabel(CS, inline=1, fontsize=10)
    ax.plot(np.arccos(0.9999),ref,'k',marker='*',ls='', ms=10)
    aux = range(1,len(corr))
    #colors = plt.matplotlib.cm.jet(np.linspace(0,1,len(corr)))
    for i in reversed(aux):
        ax.plot(np.arccos(corr[i]), std[i],c=colors[i],alpha=0.7,marker=styles[i],label="%s" %names[i])
        #ax.text(np.arccos(corr[i]), std[i],"%s"%i)
    # inset axes....
    '''axins = ax.inset_axes([1.1, 0, 0.3, 0.3])
    def pol2cart(phi, rho):
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(x, y)
    x = np.empty(len(aux))
    y = np.empty(len(aux))
    for i in reversed(aux):
        x[i-1], y[i-1] = pol2cart(np.arccos(corr[i]), std[i]) 
        axins.plot(x[i-1], y[i-1], c=colors[i],alpha=0.7,marker=styles[i],label="%s" %names[i])
        #axins.text(x[i-1], y[i-1],"%s"%i)
    # sub region of the original image
    x1, x2, y1, y2 = np.min(x)-(1/50)*np.min(x), np.max(x)+(1/50)*np.max(x),\
                     np.min(y)-(1)*np.min(y), np.max(y)+(1)*np.max(y)
    axins.set_xlim(x1,x2)
    axins.set_ylim(y1,y2)
    axins.set_xticks([])
    axins.set_yticks([])
    axins.set_xticklabels('')
    axins.set_yticklabels('')'''
    plt.legend(bbox_to_anchor=(1.5, 1),prop=dict(size='small'),loc='best')

def normalise(M):
    """ Normalize the entries of a multidimensional array sum to 1. """
    c = np.sum(M)
    # Set any zeros to one before dividing
    d = c + 1*(c==0)
    M = M/d
    return M

def hanning2d(M, N):
    """
    A 2D hanning window, as per IDL's hanning function.  See numpy.hanning for the 1d description
    """
    
    if N <= 1:
        return np.hanning(M)
    elif M <= 1:
        return np.hanning(N) # scalar unity; don't window if dims are too small
    else:
        return np.outer(np.hanning(M),np.hanning(N))

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(phi, rho)

def avg_raPsd2dv1(img3d,res,hanning):
    """ Computes and plots radially averaged power spectral density mean (power
     spectrum) of an image set img3d along the first dimension.
    """
    N = img3d.shape[0]
    for i in range(N):
        img=img3d[i,:,:]
        f_, Pf_ = raPsd2dv1(img,res,hanning)
        if i==0:
            f, Pf = f_, Pf_
        else:
            f = np.vstack((f,f_))
            Pf= np.vstack((Pf,Pf_))
    Pf = np.mean(Pf,axis=0)
    return f_, Pf

def avg_err_raPsd2dv1(img3d,img3dref,res,hanning):
    """ Computes and plots radially averaged power spectral density error mean (power
     spectrum) of an image set img3d along the first dimension.
    """
    N = img3d.shape[0]
    for i in range(N):
        img1 = img3d[i,:,:]
        img2 = img3dref[i,:,:]
        f_, Pf_ = raPsd2dv1(img1-img2,res,hanning)
        Pf_     = (Pf_/raPsd2dv1(img2,res,hanning)[1])
        if i==0:
            f, Pf = f_, Pf_
        else:
            f = np.vstack((f,f_))
            Pf= np.vstack((Pf,Pf_))
    Pf = np.mean(Pf,axis=0)
    return f_, Pf

def err_raPsd2dv1(img,imgref,res,hanning):
    """ Computes and plots radially averaged power spectral density error (power
     spectrum).
    """
    f_, Pf_ = raPsd2dv1(img-imgref,res,hanning)
    Pf_     = (Pf_/raPsd2dv1(imgref,res,hanning)[1])
    return f_, Pf_
    
def raPsd2dv1(img,res,hanning):
    """ Computes and plots radially averaged power spectral density (power
     spectrum) of image IMG with spatial resolution RES.
    """
    img = img.copy()
    N, M = img.shape
    if hanning:
        img = hanning2d(*img.shape) * img       
    img =  Imputing_NaN(img)     
    imgf = np.fft.fftshift(np.fft.fft2(img))
    imgfp = np.power(np.abs(imgf)/(N*M),2)    
    # Adjust PSD size
    dimDiff = np.abs(N-M)
    dimMax = max(N,M)
    if (N>M):
        if ((dimDiff%2)==0):
            imgfp = np.pad(imgfp,((0,0),(int(dimDiff/2),int(dimDiff/2))),'constant',constant_values=np.nan)
        else:
            imgfp = np.pad(imgfp,((0,0),(int(dimDiff/2),1+int(dimDiff/2))),'constant',constant_values=np.nan)
            
    elif (N<M):
        if ((dimDiff%2)==0):
            imgfp = np.pad(imgfp,((int(dimDiff/2),int(dimDiff/2)),(0,0)),'constant',constant_values=np.nan)
        else:
            imgfp = np.pad(imgfp,((int(dimDiff/2),1+int(dimDiff/2)),(0,0)),'constant',constant_values=np.nan)
    halfDim = int(np.ceil(dimMax/2.))
    X, Y = np.meshgrid(np.arange(-dimMax/2.,dimMax/2.-1+0.00001),np.arange(-dimMax/2.,dimMax/2.-1+0.00001))           
    theta, rho = cart2pol(X, Y)                                              
    rho = np.round(rho+0.5)   
    Pf = np.zeros(halfDim)
    f1 = np.zeros(halfDim)
    for r in range(halfDim):
      Pf[r] = np.nansum(imgfp[rho == (r+1)])
      f1[r] = float(r+1)/dimMax
    f1 = f1/res
    return f1, Pf

def plot(ax,i,j,lon,lat,data,title,extent=[-65,-55,30,40],cmap="coolwarm",gridded=True,vmin=-2,vmax=2,colorbar=True,orientation="horizontal"):
    ax[i][j].set_extent(list(extent))
    if gridded:
        im=ax[i][j].pcolormesh(lon, lat, data, cmap=cmap,\
                          vmin=vmin, vmax=vmax,edgecolors='face', alpha=1, \
                          transform= ccrs.PlateCarree(central_longitude=0.0))
    else:
        im=ax[i][j].scatter(lon, lat, c=data, cmap=cmap, s=1,\
                       vmin=vmin, vmax=vmax,edgecolors='face', alpha=1, \
                       transform= ccrs.PlateCarree(central_longitude=0.0)) 
    im.set_clim(vmin,vmax)
    if colorbar==True:
        clb = plt.colorbar(im, orientation=orientation, extend='both', pad=0.1, ax=ax[i][j])
    ax[i][j].set_title(title, pad=40, fontsize = 15)
    gl = ax[i][j].gridlines(alpha=0.5,draw_labels=True)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_bottom = False
    gl.ylabels_right = False
    gl.xlabel_style = {'fontsize': 10, 'rotation' : 45}
    gl.ylabel_style = {'fontsize': 10}
    ax[i][j].coastlines(resolution='50m')

