import numpy as np
from netCDF4 import Dataset
from sys import argv
import matplotlib.pylab as plt

input_file = argv[1]
output_file = argv[2]

# Read input
nc = Dataset(input_file, 'r')
freq = nc.variables['wavenumber'][:]
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
cxy_real = nc.variables['cross_spectrum_real'][:, :, :]
cxy_img = nc.variables['cross_spectrum_imag'][:, :, :]
nc.close()


phase_coherence = np.arctan2(cxy_img, cxy_real)*180./np.pi



nc_out = Dataset(output_file, 'w', format='NETCDF4')
nc_out.createDimension('lat', lat.size)
nc_out.createDimension('lon', lon.size)
nc_out.createDimension('wavenumber', freq.size)

freq_out = nc_out.createVariable('wavenumber', 'f8', ('wavenumber',))
freq_out[:] = np.ma.masked_invalid(freq)

phase_coherence_out = nc_out.createVariable('ratio', 'f8', ('wavenumber', 'lat', 'lon'))
phase_coherence_out[:, :, :] = np.ma.masked_invalid(phase_coherence)

lat_out = nc_out.createVariable('lat', 'f8', ('lat',))
lat_out[:] = lat
lon_out = nc_out.createVariable('lon', 'f8', ('lon',))
lon_out[:] = lon

nc_out.close()
