########################
# Plot setup

# - Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lon, lat);

# - Continent color
contcolor = '0.5';

# - Map proj
mapproj = 'cea'; # choose this more carefully later

# - Lat and lon
wlon = 100; elon = 290;
slat = -70; nlat = 70;
parallels = np.arange(-60,61,30);
meridians = np.arange(100,291,20);
latlonlabsize = 8;
lonlabrot = 45;

# - Colorbar
cbarlabsize = 8;

# - Define common O2 levels
#levsnow = np.linspace(200,850,11,endpoint=True)
#difflevsnow = np.linspace(-100,100,11, endpoint=True)

########################
# Plot seasonal SKJ and BET catch
fig = plt.figure(figsize=(14,9))

iseas=0;
#for isp in np.arange(1,9,2):
for isp in np.arange(1,5):
    s1 = fig.add_subplot(2,4,isp)
    s1.set_title("SKJ catch, season #" + repr(iseas))
    m = Basemap(projection=mapproj,resolution='c',\
            llcrnrlon=wlon, urcrnrlon=elon,\
            llcrnrlat=slat,urcrnrlat=nlat)
    x, y = m(lon2d, lat2d);
    m.drawcoastlines();
    m.fillcontinents(color=contcolor);
    m.drawmapboundary();
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=latlonlabsize);
    merids = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=latlonlabsize);
    for im in merids:
        try:
            merids[im][1][0].set_rotation(lonlabrot)
        except:
            pass
    cs = m.pcolor(x, y, skj_c_tot_seas[iseas,:,:]);
    cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=cbarlabsize)
    cbar.set_label("Catch [metric tonnes]")

    s2 = fig.add_subplot(2,4,isp+4)
    s2.set_title("BET catch, season #" + repr(iseas))
    m = Basemap(projection=mapproj,resolution='c',\
            llcrnrlon=wlon, urcrnrlon=elon,llcrnrlat=slat,urcrnrlat=nlat)
    x, y = m(lon2d, lat2d);
    m.drawcoastlines();
    m.fillcontinents(color=contcolor);
    m.drawmapboundary();
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=latlonlabsize);
    merids = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=latlonlabsize);
    for im in merids:
        try:
            merids[im][1][0].set_rotation(lonlabrot)
        except:
            pass
    cs = m.pcolor(x, y, bet_c_tot_seas[iseas,:,:]);
    cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=cbarlabsize)
    cbar.set_label("Catch [metric tonnes]")

    iseas = iseas + 1;

fig.show()
#fig.savefig('skjbet_totmeancatch_pcolor.png',dpi=200)
