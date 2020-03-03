import xarray as xr
import numpy as np
dpath = '/opt/skipjack-bigeye-separation/data/'
ncfile = dpath + 'ecco/MXLDEPTH.0001.nc'
mldds = xr.open_dataset(ncfile)
mld = mldds['MXLDEPTH']
mld['i3'] = mld['lon'].isel(i2=0).values
mld['i2'] = mld['lat'].isel(i3=0).values
mld['i1'] = mld['tim'].values
mld = mld.drop(['lat','lon','tim'])
mld = mld.rename({'i1':'time', 'i2':'lat', 'i3':'lon'})
mld['lon'] = np.mod(mld['lon'], 360)
mld = mld.reindex({'lon': np.sort(mld['lon'])}) # or: mld = mld.roll(lon=36)
mld.to_netcdf(dpath+'ecco/MXLDEPTH.0001.correcteddims.nc')