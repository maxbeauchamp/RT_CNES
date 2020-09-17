import xarray as xr
import sys
import matplotlib.pylab as plt
from scipy.interpolate import RegularGridInterpolator
import numpy as np

ds = xr.open_dataset(sys.argv[1])
cvar = sys.argv[2]

depth_criteria = 500 # m
var_criteria = 0.01 # m2

ds_bathy = xr.open_dataset('/home/ad/ballarm/scratch/data/MiostA_201909/sad/Rossby_radius.nc')

bathy = np.abs(np.ma.masked_invalid(ds_bathy['H'].values).filled(0.))
lon_bathy = ds_bathy['lon'].values
lat_bathy = ds_bathy['lat'].values
finterp_bathy = RegularGridInterpolator((lon_bathy, lat_bathy), bathy, bounds_error=False, fill_value=None)

lon2d, lat2d = np.meshgrid(ds.lon.values, ds.lat.values)
bathy_interp = finterp_bathy((lon2d, lat2d))

#plt.pcolormesh(np.ma.masked_where(bathy_interp > depth_criteria, bathy_interp))
#plt.show()

data = ds[str(cvar)].values
data = np.ma.masked_outside(data, -10, 10)

# Continental_plateau
data2 = np.ma.masked_where((bathy_interp > depth_criteria) | (bathy_interp <= 0), data)
print(np.ma.mean(np.ma.masked_invalid(data2)))
plt.pcolormesh(data2)
plt.colorbar()
plt.show()

# off shore low variability
data2 = np.ma.masked_where((bathy_interp < depth_criteria) | (ds['variance_ref'].values > var_criteria), data)
print(np.ma.mean(np.ma.masked_invalid(data2)))
plt.pcolormesh(data2)
plt.colorbar()
plt.show()

# off shore high variability
data2 = np.ma.masked_where((bathy_interp < depth_criteria) | (ds['variance_ref'].values <= var_criteria), data)
print(np.ma.mean(np.ma.masked_invalid(data2)))
plt.pcolormesh(data2)
plt.colorbar()
plt.show()
