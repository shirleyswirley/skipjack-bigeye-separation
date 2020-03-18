# - Draw coastlines and/or land
ax.coastlines(linewidth=0.4)
#ax.add_feature(cartopy.feature.COASTLINES, facecolor='black')
#ax.add_feature(cartopy.feature.LAND, facecolor='gray')

# - Define gridlines
xlocsnow = np.concatenate(
    [np.linspace(120,180,3),
     np.linspace(-180,-90,4)])
# ...needs 180,-180 repeated to make lat horiz lines go through
ylocsnow = np.linspace(-60, 40, 6)
ylocsnow = np.append(ylocsnow,90)
ylocsnow = np.insert(ylocsnow,0,-90)
# ...2 above lines needed to make lon vert lines go through

# - Draw gridlines
gl = ax.gridlines(
    crs=ccrs.PlateCarree(),
    draw_labels=False, # take care of labels below
    linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlocator = mticker.FixedLocator(xlocsnow)
gl.ylocator = mticker.FixedLocator(ylocsnow)

# - Add lat and lon labels
xticksnow = np.concatenate(
    [np.linspace(120,180,3),
     np.linspace(-150,-90,3)])
# ...leave out one repeated 180,-180
yticksnow = np.linspace(-60, 40, 6)
ax.set_xticks(xticksnow, crs=ccrs.PlateCarree())
ax.set_yticks(yticksnow, crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())

# - Set extent
ax.set_extent([105, 235, -65, 55.25], crs=ccrs.PlateCarree())
