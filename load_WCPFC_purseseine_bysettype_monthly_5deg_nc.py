ncfile='/ltraid2/sleung/data/WCPFC/5by5deg/ByMonth/PurseSeine_0/WCPFC_purseseine_bysettype_monthly_5deg.nc';
fh = nc.Dataset(ncfile, mode='r');

lon = fh.variables['lon'][:];
lat = fh.variables['lat'][:];
time = fh.variables['time'][:]; # monthly data starts jan-1967
timemocatch = pd.date_range('1967-01-01', periods=612, freq='MS')
#np.array(timemocatch)

skj_c_una = fh.variables['skj_c_una'][:];
skj_c_log = fh.variables['skj_c_log'][:];
skj_c_dfad = fh.variables['skj_c_dfad'][:];
skj_c_afad = fh.variables['skj_c_afad'][:];
skj_c_oth = fh.variables['skj_c_oth'][:];

bet_c_una = fh.variables['bet_c_una'][:];
bet_c_log = fh.variables['bet_c_log'][:];
bet_c_dfad = fh.variables['bet_c_dfad'][:];
bet_c_afad = fh.variables['bet_c_afad'][:];
bet_c_oth = fh.variables['bet_c_oth'][:];

sets_una = fh.variables['sets_una'][:];
sets_log = fh.variables['sets_log'][:];
sets_dfad = fh.variables['sets_dfad'][:];
sets_afad = fh.variables['sets_afad'][:];
sets_oth = fh.variables['sets_oth'][:];

# time, y, x
fh.close()

# print(lon)
# print(lat)
# print(time)
# print(np.shape(skj_c_una))
