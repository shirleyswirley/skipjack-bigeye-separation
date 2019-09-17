# to run this script:
# (https://realpython.com/run-python-scripts/)
#import runpy
#runpy.run_module(mod_name='wcpfc')

# THIS SITE, YES!!
# https://stackoverflow.com/questions/436198/what-is-an-alternative-to-execfile-in-python-3

######################################################
# LOAD LIBRARIES
######################################################
import netCDF4 as nc
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from reverse_colormap import reverse_colormap 
import cmocean
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid

# trying stuff:
# from importlib import reload 
# import os
# import sys
# sys.path.append('/home/test/')
# os.chdir("/ltraid2/sleung/tunaextremes/pyfiles/wcpfc")
# cwd = os.getcwd()
# %run /path/to/filename.py

######################################################
# LOAD MONTHLY 5 BY 5 DEGREE WCPFC CATCH AND EFFORT DATA
######################################################
# - Variables loaded are: lon, lat, time, timemocatch,
# skj_c_una, skj_c_log, skj_c_dfad, skj_c_afad, skj_c_oth,
# bet_c_una, bet_c_log, bet_c_dfad, bet_c_afad, bet_c_oth,
# sets_una, sets_log, sets_dfad, sets_afad, sets_oth
exec(open("load_WCPFC_purseseine_bysettype_monthly_5deg_nc.py").read())

######################################################
# LOAD MONTHLY ENSO TIME SERIES
######################################################
# - Variables loaded/calced are: oni, onien, oniln, onitime
plotensoidx=1
begdate = pd.DatetimeIndex(['1967-1-1']) 
enddate = pd.DatetimeIndex(['2017-12-1'])
onitime = pd.date_range('1967-01-01', periods=612, freq='MS')
exec(open("load_ONI_txt.py").read())

######################################################
# CALCULATE WCPFC TOTAL CATCH, CPUE, EFFORT,
# BET:SKJ CATCH RATIO
######################################################
# - Variables calced are: skj_c_tot, bet_c_tot, sets_tot,
# skj_cp_tot, bet_cp_tot, bettoskj_c_tot
skj_c_tot = skj_c_una + skj_c_log + skj_c_dfad + skj_c_afad + skj_c_oth
bet_c_tot = bet_c_una + bet_c_log + bet_c_dfad + bet_c_afad + bet_c_oth
sets_tot = sets_una + sets_log + sets_dfad + sets_afad + sets_oth
skj_cp_tot = skj_c_tot/sets_tot
bet_cp_tot = bet_c_tot/sets_tot
bettoskj_c_tot = bet_c_tot/skj_c_tot

######################################################
# CALCULATE MEAN AND ENSO COMPOSITES OF TOTAL CATCH,
# CPUE, EFFORT, BET:SKJ CATCH RATIO
######################################################
skj_c_tot_mean = np.mean(skj_c_tot,axis=0)
skj_c_tot_en = np.mean(np.squeeze(skj_c_tot[np.where(onien),:,:]),axis=0)
skj_c_tot_ln = np.mean(np.squeeze(skj_c_tot[np.where(oniln),:,:]),axis=0)

bet_c_tot_mean = np.mean(bet_c_tot,axis=0)
bet_c_tot_en = np.mean(np.squeeze(bet_c_tot[np.where(onien),:,:]),axis=0)
bet_c_tot_ln = np.mean(np.squeeze(bet_c_tot[np.where(oniln),:,:]),axis=0)

skj_cp_tot_mean = np.mean(skj_cp_tot,axis=0)
skj_cp_tot_en = np.mean(np.squeeze(skj_cp_tot[np.where(onien),:,:]),axis=0)
skj_cp_tot_ln = np.mean(np.squeeze(skj_cp_tot[np.where(oniln),:,:]),axis=0)

bet_cp_tot_mean = np.mean(bet_cp_tot,axis=0)
bet_cp_tot_en = np.mean(np.squeeze(bet_cp_tot[np.where(onien),:,:]),axis=0)
bet_cp_tot_ln = np.mean(np.squeeze(bet_cp_tot[np.where(oniln),:,:]),axis=0)

sets_tot_mean = np.mean(sets_tot,axis=0)
sets_tot_en = np.mean(np.squeeze(sets_tot[np.where(onien),:,:]),axis=0)
sets_tot_ln = np.mean(np.squeeze(sets_tot[np.where(oniln),:,:]),axis=0)

bettoskj_c_tot_mean = np.mean(bettoskj_c_tot,axis=0)
bettoskj_c_tot_en = np.mean(np.squeeze(bettoskj_c_tot[np.where(onien),:,:]),axis=0)
bettoskj_c_tot_ln = np.mean(np.squeeze(bettoskj_c_tot[np.where(oniln),:,:]),axis=0)

#test = np.random.rand(3,4,2)
#testlog = [0,1,0]
#test1=test[np.where(testlog),:,:]
#test2=np.squeeze(test1)

######################################################
######################################################
# START PLOTTING HERE
######################################################
######################################################

########################
# Define Basemap plot parameters
# --> parameters used in setup_Basemap.py

# - imshow vs. pcolor vs. pcolormesh
# https://thomas-cokelaer.info/blog/2014/05/matplotlib-difference-between-pcolor-pcolormesh-and-imshow/
# --> use pcolor instead of pcolormesh b/c
# sometimes pcolormesh distorts w/ Basemap? 

# - Define common levels + extend colorbar example
#levsnow = np.linspace(200,850,11,endpoint=True)
#difflevsnow = np.linspace(-100,100,11, endpoint=True)
#cs = m.contourf(x, y, p50yftwoa, levsnow, cmap=cmocean.cm.deep, extend='both')

# - Create 2D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lon, lat)

# - Continent color
contcolor = '0.5'

# - Map proj
mapproj = 'cea' # choose this more carefully later

# - Lat and lon
wlon = 100; elon = 290
slat = -70; nlat = 70
parallels = np.arange(-60,61,30)
meridians = np.arange(100,291,20)
latlonlabsize = 8; lonlabrot = 45

# - Colorbar
cbarlabrot = 45; cbarlabsize = 8

# - Define colormaps
cmocean.cm.balance_r = reverse_colormap(cmocean.cm.balance)
viridis_r = reverse_colormap(plt.get_cmap('viridis'))

######################################################
# PLOT MEAN SKJ CATCH, BET CATCH, BET:SKJ CATCH RATIO
######################################################
fig = plt.figure(figsize=(3,13))

s1 = fig.add_subplot(311)
s1.set_title("Skipjack mean catch");
exec(open("setup_Basemap.py").read())
cs = m.pcolor(x, y, skj_c_tot_mean, cmap='viridis')
#cs = m.pcolor(x, y, skj_c_tot_mean, cmap=cmocean.cm.deep)
exec(open("colorbar_Basemap.py").read())
cbar.set_label("Catch [metric tonnes]")

s2 = fig.add_subplot(312)
s2.set_title("Bigeye mean catch");
exec(open("setup_Basemap.py").read())
cs = m.pcolor(x, y, bet_c_tot_mean)
exec(open("colorbar_Basemap.py").read())
cbar.set_label('Catch [metric tonnes]')

s3 = fig.add_subplot(313)
s3.set_title("Bigeye:skipjack catch");
exec(open("setup_Basemap.py").read())
cs = m.pcolor(x, y, bettoskj_c_tot_mean)
exec(open("colorbar_Basemap.py").read())
cbar.set_label('[unitless]')

fig.show()
#fig.savefig('skjbet_totmeancatch_contourf.png',dpi=200)

########################
# pcolor plot
fig = plt.figure(figsize=(11,4))

s1 = fig.add_subplot(121)
s1.set_title("Skipjack mean catch")
m = Basemap(projection=mapproj,resolution='c',\
        llcrnrlon=wlon, urcrnrlon=elon,llcrnrlat=slat,urcrnrlat=nlat)
x, y = m(lon2d, lat2d)
m.drawcoastlines();
#m.fillcontinents(color=contcolor);
m.drawmapboundary();
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=latlonlabsize);
merids = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=latlonlabsize);
for im in merids:
    try:
        merids[im][1][0].set_rotation(lonlabrot)
    except:
        pass

cs = m.pcolor(x, y, skj_c_tot_mean);
cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=cbarlabsize)
cbar.set_label("Catch [metric tonnes]")

s3 = fig.add_subplot(122)
s3.set_title("Bigeye mean catch")
m = Basemap(projection=mapproj,resolution='c',\
        llcrnrlon=wlon, urcrnrlon=elon,llcrnrlat=slat,urcrnrlat=nlat)
x, y = m(lon2d, lat2d)
m.drawcoastlines()
m.fillcontinents(color=contcolor)
m.drawmapboundary()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=latlonlabsize)
merids = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=latlonlabsize)
for im in merids:
    try:
        merids[im][1][0].set_rotation(lonlabrot)
    except:
        pass

cs = m.pcolor(x, y, bet_c_tot_mean)
cbar = fig.colorbar(cs, orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=cbarlabsize) 
cbar.set_label('Catch [metric tonnes]')

fig.show()
#fig.savefig('skjbet_totmeancatch_pcolor.png',dpi=200)

######################################################
# PLOT CORRELATIONS BTWN MONTHLY SKJ AND BET CATCH
######################################################
betvsskj_c_tot_cc = np.full((np.size(lat),np.size(lon)),fill_value=np.nan)
for ilon in range(np.size(lon)):
    for ilat in range(np.size(lat)):
        betvsskj_c_tot_cc[ilat,ilon] = np.ma.corrcoef(skj_c_tot[:,ilat,ilon],bet_c_tot[:,ilat,ilon])[0,1];

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
# pcolor plot
fig = plt.figure(figsize=(11,4));

plt.title("SKJ vs. BET total purse seine catch, correlation coeff")
m = Basemap(projection=mapproj,resolution='c',\
        llcrnrlon=wlon, urcrnrlon=elon,llcrnrlat=slat,urcrnrlat=nlat)
x, y = m(lon2d, lat2d);
m.drawcoastlines();
#m.fillcontinents(color=contcolor);
m.drawmapboundary();
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=latlonlabsize);
merids = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=latlonlabsize);
for im in merids:
    try:
        merids[im][1][0].set_rotation(lonlabrot)
    except:
        pass

cs = m.pcolor(x, y, betvsskj_c_tot_cc);
cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=cbarlabsize)
cbar.set_label("Correl coeff")

fig.show()
#fig.savefig('skjbet_totmeancatch_pcolor.png',dpi=200)

######################################################
# LOOK AT RELNSHIPS/SCATTERPLOTS BTWN CATCH AND ENV VARIABLES
######################################################
ncfile = '/ltraid2/sleung/data/WOD18/ncfiles/new1955to2018/straddleequator/o2_195501-201807_5deg.nc';
fh = nc.Dataset(ncfile, mode='r');
lonwod = fh.variables['lon'][:];
latwod = fh.variables['lat'][:];
depthwod = fh.variables['depth'][:];
timewod = fh.variables['time'][:]; # monthly data starts jan-1967
timemowod = pd.date_range('1955-01-01', periods=763, freq='MS')

o2 = fh.variables['O2'][:];
fh.close()

ncfile = '/ltraid2/sleung/data/WOD18/ncfiles/new1955to2018/straddleequator/temp_195501-201807_5deg.nc';
fh = nc.Dataset(ncfile, mode='r');
temp = fh.variables['Temp'][:];
fh.close()

ncfile = '/ltraid2/sleung/data/WOD18/ncfiles/new1955to2018/straddleequator/sal_195501-201807_5deg.nc';
fh = nc.Dataset(ncfile, mode='r');
sal = fh.variables['Sal'][:];
fh.close()

# reshape lonwod to 0-360
lonwodlen = np.size(lonwod);
lonwod = np.concatenate((lonwod[np.int(lonwodlen/2):lonwodlen],lonwod[0:np.int(lonwodlen/2)]),axis=0);
lonwod[lonwod<0] = lonwod[lonwod<0]+360;
# reshape o2 to 0-360
masknow = o2.mask;
o2 = np.concatenate((o2[:,:,:,np.int(lonwodlen/2):lonwodlen],\
        o2[:,:,:,0:np.int(lonwodlen/2)]),axis=3)
o2.mask = np.concatenate((masknow[:,:,:,np.int(lonwodlen/2):lonwodlen],\
        masknow[:,:,:,0:np.int(lonwodlen/2)]),axis=3)
# reshape temp to 0-360
masknow = temp.mask;
temp = np.concatenate((temp[:,:,:,np.int(lonwodlen/2):lonwodlen],\
        temp[:,:,:,0:np.int(lonwodlen/2)]),axis=3)
temp.mask = np.concatenate((masknow[:,:,:,np.int(lonwodlen/2):lonwodlen],\
        masknow[:,:,:,0:np.int(lonwodlen/2)]),axis=3)
# reshape sal to 0-360
masknow = sal.mask;
sal = np.concatenate((sal[:,:,:,np.int(lonwodlen/2):lonwodlen],\
        sal[:,:,:,0:np.int(lonwodlen/2)]),axis=3)
sal.mask = np.concatenate((masknow[:,:,:,np.int(lonwodlen/2):lonwodlen],\
        masknow[:,:,:,0:np.int(lonwodlen/2)]),axis=3)

# testing
#a = sp.randn(5,4,2,3)
#a1 = np.concatenate((a[:,:,:,1:3],a[:,:,:,0:1]),axis=3)

# get concurrent times
tidx_wodi = np.int(np.where(timemowod=='1967-01-01')[0]);
tidx_wodf = np.int(np.where(timemowod=='2017-12-01')[0])+1;
#tidx_cati = np.int(np.where(timemocatch=='1967-01-01')[0]);
#tidx_catf = np.int(np.where(timemocatch=='2017-12-01')[0]);

# get spatial idxs to keep to compare to catch data
# (should code this better...)
lonwodc_keepidx = np.arange(np.int(np.where(lonwod==np.min(lon))[0]),np.int(np.where(lonwod==np.max(lon))[0])+1);
latwodc_keepidx = np.arange(np.int(np.where(latwod==np.min(lat))[0]),np.int(np.where(latwod==np.max(lat))[0])+1);

# trim vars to catch area 
latwodc = latwod[latwodc_keepidx];
lonwodc = lonwod[lonwodc_keepidx];
o2c = o2[tidx_wodi:tidx_wodf,:,latwodc_keepidx,:];
o2c = o2c[:,:,:,lonwodc_keepidx];
tempc = temp[tidx_wodi:tidx_wodf,:,latwodc_keepidx,:];
tempc = tempc[:,:,:,lonwodc_keepidx];
salc = sal[tidx_wodi:tidx_wodf,:,latwodc_keepidx,:];
salc = salc[:,:,:,lonwodc_keepidx];

o2100c = np.squeeze(o2c[:,np.where(depthwod==100)[0],:,:]);
sstc = tempc[:,0,:,:];
timemocatch3d = np.full((np.size(timemocatch),np.size(lat),np.size(lon)), fill_value='nan')
for itime in range(np.size(timemocatch)):
    timemocatch3d[itime,:,:] = timemocatch[itime].month;

########################
# Plot setup (look at mean o2100 and sst)

# - Create 2D lat/lon arrays for Basemap
#lon2d, lat2d = np.meshgrid(lon, lat);
lon2d, lat2d = np.meshgrid(lonwod, latwod);

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
# pcolor plot (look at mean o2100 and sst)
fig = plt.figure(figsize=(11,4))

s1 = fig.add_subplot(121)
s1.set_title("Mean O2 at 100 m")
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

cs = m.pcolor(x, y, np.mean(np.squeeze(o2[:,np.where(depthwod==100)[0],:,:]),axis=0));
cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=cbarlabsize)
cbar.set_label("[umol/kg]")

s2 = fig.add_subplot(122)
s2.set_title("Mean SST")
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

cs = m.pcolor(x, y, np.mean(temp[:,0,:,:],axis=0));
cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=cbarlabsize)
cbar.set_label("[degC]")

fig.show()
#fig.savefig('skjbet_totmeancatch_pcolor.png',dpi=200)

########################
# Plot scatterplot btwn SKJ catch and o2100, sst
fig = plt.figure(figsize=(11,4))

s1 = fig.add_subplot(121)
s1.set_title("Monthly SKJ catch vs. O2 at 100 m")
plt.scatter(o2100c,skj_c_tot,c=timemocatch3d);

s2 = fig.add_subplot(122)
s2.set_title("Monthly SKJ catch vs. SST")
plt.scatter(sstc,skj_c_tot,c=timemocatch3d);

fig.show()
#fig.savefig('skjbet_totmeancatch_pcolor.png',dpi=200)

########################
# Plot scatterplot btwn BET catch and o2100, sst
fig = plt.figure(figsize=(11,4))

s1 = fig.add_subplot(121)
s1.set_title("Monthly BET catch vs. O2 at 100 m")
plt.scatter(o2100c,bet_c_tot,c=timemocatch3d);

s2 = fig.add_subplot(122)
s2.set_title("Monthly BET catch vs. SST")
plt.scatter(sstc,bet_c_tot,c=timemocatch3d);

fig.show()
#fig.savefig('skjbet_totmeancatch_pcolor.png',dpi=200)

######################################################
# LOOK AT CLIMATOLOGIES OF SPATIAL CATCH DISTRIBUTIONS
######################################################
skj_c_tot_clim = np.full((12,np.size(lat),np.size(lon)), fill_value=np.nan)
bet_c_tot_clim = np.full((12,np.size(lat),np.size(lon)), fill_value=np.nan)

# this assumes that skj_c_tot, bet_c_tot start in jan
for imonth in np.arange(0,12):
    skj_c_tot_clim[imonth,:,:] = np.mean(skj_c_tot[np.arange(imonth,np.size(timemocatch),12),:,:],axis=0);
    bet_c_tot_clim[imonth,:,:] = np.mean(bet_c_tot[np.arange(imonth,np.size(timemocatch),12),:,:],axis=0);

skj_c_tot_seas = np.full((4,np.size(lat),np.size(lon)), fill_value=np.nan);
bet_c_tot_seas = np.full((4,np.size(lat),np.size(lon)), fill_value=np.nan);

# this assumes that skj_c_tot, bet_c_tot start in jan
skj_c_tot_seas[0,:,:] = np.mean(skj_c_tot_clim[np.array([0,1,11]),:,:],axis=0);
iseas = 1;
for imonth in np.arange(2,9,3):
    skj_c_tot_seas[iseas,:,:] = np.mean(skj_c_tot_clim[imonth:imonth+3,:,:],axis=0);
    bet_c_tot_seas[iseas,:,:] = np.mean(bet_c_tot_clim[imonth:imonth+3,:,:],axis=0);
    iseas = iseas+1;

########################
# Plot climatological SKJ and BET catch

#exec(open("wcpfc_seasonal.py").read(), globals())

########################
# Plot seasonal SKJ and BET catch
exec(open("wcpfc_seasonal.py").read(), globals())
