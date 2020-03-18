# - From https://www.wcpfc.int/convention-text:
# The Convention Area is defined in article 3 of
# the Convention and comprises all
# waters of the Pacific Ocean bounded
# to the south and to the east by a line drawn
# from the south coast of Australia due south
# along the 141° meridian of east
# longitude to its intersection with the 55°
# parallel of south latitude;
# thence due east along the 55° parallel of
# south latitude to its intersection
# with the 150° meridian of east longitude;
# thence due south along the 150°
# meridian of east longitude to its intersection
# with the 60° parallel of south latitude;
# thence due east along the 60° parallel of south latitude
# to its intersection with the 130° meridian of west longitude;
# thence due north along the 130° meridian of west longitude
# to its intersection with the 4° parallel of south latitude;
# thence due west along the 4° parallel of south latitude
# to its intersection with the 150° meridian of west longitude;
# thence due north along the 150° meridian of west longitude.

latptsnow = [-38, -55, -55, -60, -60, -4, -4, 60]
lonptsnow = [141, 141, 150, 150, 180+50, 180+50, 180+30, 180+30]
ax.plot(lonptsnow, latptsnow,
         color='red', linewidth=1,
         transform=ccrs.PlateCarree())

# - Note: normally you would do the following,
# but not needed w/ PlateCarree transform
#lat0, lon0 = -38, 141
#lon0t, lat0t = ccrs.PlateCarree().transform_point(lon0, lat0, ccrs.Geodetic())