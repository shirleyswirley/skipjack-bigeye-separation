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
