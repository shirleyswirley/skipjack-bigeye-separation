# - Add coastlines
ax.coastlines();

# - Add gridlines
xlocsnow = np.concatenate([np.linspace(120,180,3),np.linspace(-180,-90,4)]) # needs 180,-180 repeated to make y lines go through
xlocsnow = np.append(xlocsnow,-70) # for some reason -90 won't plot w/o this
ylocsnow = np.linspace(-80,80,7)
ylocsnow = np.append(ylocsnow,90); ylocsnow = np.insert(ylocsnow,0,-90) # to make x lines go through
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, # take care of labels below
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.xlocator = mticker.FixedLocator(xlocsnow); gl.ylocator = mticker.FixedLocator(ylocsnow)

# - Add lat and lon labels
xticksnow = np.concatenate([np.linspace(120,180,3),np.linspace(-150,-90,3)]) # get rid of repeated 180,-180
yticksnow = np.linspace(-80,80,9) # get rid of -90,90
ax.set_xticks(xticksnow, crs=ccrs.PlateCarree())
ax.set_yticks(yticksnow, crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())