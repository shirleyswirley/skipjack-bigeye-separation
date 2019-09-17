%----------------------------------------------------
% Load environmental variables from WOD
%----------------------------------------------------
ncfile = '/ltraid2/sleung/data/WOD18/ncfiles/new1955to2018/straddleequator/o2_195501-201807_5deg.nc';
ncid = netcdf.open(ncfile,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'lon');
wodlon = netcdf.getVar(ncid,varid,'double');
varid = netcdf.inqVarID(ncid,'lat');
wodlat = netcdf.getVar(ncid,varid,'double');
varid = netcdf.inqVarID(ncid,'depth');
woddepth = netcdf.getVar(ncid,varid,'double');
wodplat = wodlat-5/2;
wodplon = wodlon-5/2;
timemowod = datetime('1955-01-01') + calmonths(0:762);

varid = netcdf.inqVarID(ncid,'O2');
o2 = netcdf.getVar(ncid,varid,'double');
o2(o2<-1e33)=nan;
netcdf.close(ncid);

ncfile = '/ltraid2/sleung/data/WOD18/ncfiles/new1955to2018/straddleequator/temp_195501-201807_5deg.nc';
ncid = netcdf.open(ncfile,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'Temp');
temp = netcdf.getVar(ncid,varid,'double');
temp(temp<-1e33)=nan;
netcdf.close(ncid);

ncfile = '/ltraid2/sleung/data/WOD18/ncfiles/new1955to2018/straddleequator/sal_195501-201807_5deg.nc';
ncid = netcdf.open(ncfile,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'Sal');
sal = netcdf.getVar(ncid,varid,'double');
sal(sal<-1e33)=nan;
netcdf.close(ncid);

% reshape wodlon to 0-360
wodlonlen = length(wodlon);
wodlon = [wodlon(wodlonlen/2+1:wodlonlen); wodlon(1:wodlonlen/2)];
wodlon(wodlon<0) = wodlon(wodlon<0)+360;
% reshape o2 to 0-360
o2 = cat(1,o2(wodlonlen/2+1:wodlonlen,:,:,:),o2(1:wodlonlen/2,:,:,:));
% reshape temp to 0-360
temp = cat(1,temp(wodlonlen/2+1:wodlonlen,:,:,:),temp(1:wodlonlen/2,:,:,:));
% reshape sal to 0-360
sal = cat(1,sal(wodlonlen/2+1:wodlonlen,:,:,:),sal(1:wodlonlen/2,:,:,:));

% get concurrent times
tidx_wod_wci = find(timemowod==wtimemocatch(1));
tidx_wod_wcf = find(timemowod==wtimemocatch(end));

% get spatial idxs to keep to compare to catch data
% (should code this better...)
wodwlon_keepidx = find(wodlon==min(wcatchlon)):find(wodlon==max(wcatchlon));
wodwlat_keepidx = find(wodlat==min(wcatchlat)):find(wodlat==max(wcatchlat));

% trim vars to catch area 
wo2 = o2(wodwlon_keepidx,wodwlat_keepidx,:,tidx_wod_wci:tidx_wod_wcf);
wtemp = temp(wodwlon_keepidx,wodwlat_keepidx,:,tidx_wod_wci:tidx_wod_wcf);
wsal = sal(wodwlon_keepidx,wodwlat_keepidx,:,tidx_wod_wci:tidx_wod_wcf);

wo2100 = squeeze(wo2(:,:,woddepth==100,:));
wsst = squeeze(wtemp(:,:,1,:));

%----------------------------------------------------
% Load environmental variables from satellite
%----------------------------------------------------
% - Satellite sea surface height anomaly (SSHA)
% 1. Satellite sea surface height data from JSON-2 and others.
% 2. Monthly means, spatial resolution is 0.25 degree by 0.25 degree.
% 3. Located at /ltraid4/sleung/data/aviso/monthly_means_tropicalpacificonly.mat
% with the original netcdf /ltraid4/observations/aviso/gridded/monthly and
% orig matfile in /ltraid4/observations/aviso/gridded/.
% 4. Data source: http://www.aviso.altimetry.fr/en/data/products.html.
% 5. Spans Jan 1993 - Dec 2015.
% 6. Units are meters.
% 7. Variables become sshlon, sshlat, sshtime, ssh (change to ssha for clarity).

% - Satellite sea surface temperature (TS)
% 1. Satellite sea surface temperature data from reanalysis.
% 2. Monthly means, spatial resolution is 1 degree by 1 degree.
% 3. Located in /ltraid4/observations/OAFLUX/ts_*, but not very well-organized.
% 4. Data source: http://oaflux.whoi.edu/heatflux.html).
% 5. Spans Jan 1958 - Dec 2015.
% 6. Units are degC.
% 7. Variables become tslon, tslat, tstime, ts

% - Satellite net primary productivity (NPP)
% 1. Satellite primary productivity data from SeaWiFS satellite
% 2. Monthly means, spatial resolution is 1 degree by 1 degree.
% 3. Located in /ltraid4/sleung/data/SeaWiFS;
% VGPM_1997-2010.nc and cbpm_1997-2010.nc represent primary production
% products from two different algorithms
% 4. Data source: http://orca.science.oregonstate.edu/1080.by.2160.monthly.hdf.vgpm.s.chl.a.sst.php 
% 5. Spans Oct 1997 - March 2010 (VGPM),
% Jan 1997 - Dec 2010 (CbPM, though data only exists from Oct 1997 - May 2010).
% 6. Units are mgC/m^2/day.
% 7. Variables become nppvlon, nppvlat, nppvtime, nppv, nppclon, nppclat, nppctime, nppc

% - Satellite CCMP winds
% 1. Cross-Calibrated Multi-Platform (CCMP) Ocean Surface Wind Vector Analyses 
% 2. Monthly means, spatial res is 1/4 by 1/4 degree
% 3. Located at /ltraid4/sleung/data/ccmpwinds/ccmpwind_monthlymeans_tropicalpacificonly.mat 
% 4. Data source: https://podaac.jpl.nasa.gov/Cross-Calibrated_Multi-Platform_OceanSurfaceWindVectorAnalyses
% 5. Spans Aug 1987 - Dec 2011.
% 6. Units and variables are as follows:
% uwnd = u-wind at 10 m [m/s]
% vwnd = v-wind at 10 m [m/s]
% wspd = wind speed at 10 m [m/s]
% upstr = u-component of pseudostress at 10 m [m^2/s^2]
% vpstr = v-component of pseudostress at 10 m [m^2/s^2]
% 7. Variables become windlon, windlat, windtime,
% and those listed in 6.

%---------
% SSH
%---------
load('/ltraid2/sleung/data/aviso/SSHA_monthlymeans_tropicalpacificonly.mat');
sshlat = latnew; sshlon = lonnew;
sshamnew = permute(sshamnew,[2 1 3]); ssha = sshamnew;
sshtime = datetime('1993-01-01') + calmonths(0:275);

finevar = ssha; finelon = sshlon; finelat = sshlat;
finetime = sshtime; finelondegres = 0.25; finelatdegres = 0.25;
coarselon = wcatchlon; coarselat = wcatchlat;
coarseplon = wcatchplon; coarseplat = wcatchplat; plottocheck = 0;
[coarsevar] = fine2coarsegridbin(finevar,finelon,finelat,finetime,finelondegres,finelatdegres,coarselon,coarselat,coarseplon,coarseplat,plottocheck);
wssha = coarsevar;

tidx_wc_sshi = find(wtimemocatch==sshtime(1));
tidx_wc_sshf = find(wtimemocatch==sshtime(end));

%---------
% TS
%---------
yr1ts = 1958; yr2ts = 2015;
fnamebase = '/ltraid4/observations/OAFLUX/';
for iyear = yr1ts:yr2ts
    fname = [fnamebase 'ts_oaflux_' num2str(iyear) '.nc'];
    ncid = netcdf.open(fname,'nc_nowrite');
    if iyear==yr1ts
        varid = netcdf.inqVarID(ncid,'lon');
        tslon = netcdf.getVar(ncid,varid,'double');
        varid = netcdf.inqVarID(ncid,'lat');
        tslat = netcdf.getVar(ncid,varid,'double');
        ts = nan(length(tslon),length(tslat),12*(yr2ts-yr1ts+1));
        tidx = 1;
    end
    varid = netcdf.inqVarID(ncid,'tmpsf');
    ts(:,:,tidx:tidx+11) = netcdf.getVar(ncid,varid,'double');
    tidx = tidx+12;
end
ts(ts==32766)=nan; % missing value = 32766
ts=0.01.*ts; % scale factor = 0.01
tstime = datetime('1958-01-01') + calmonths(0:695);

%[coarsevar] = fine2coarsegridbin(finevar,finelon,finelat,finetime,finelondegres,finelatdegres,coarselon,coarselat,coarseplon,coarseplat,plottocheck);
[wts] = fine2coarsegridbin(ts,tslon,tslat,tstime,1,1,wcatchlon,wcatchlat,wcatchplon,wcatchplat,1);

tidx_ts_wci = find(tstime==wtimemocatch(1));
tidx_wc_tsf = find(wtimemocatch==tstime(end));

wts = wts(:,:,tidx_ts_wci:end);
tstime = datetime('1967-01-01') + calmonths(0:587);

%---------
% NPP
%---------
% - VGPM
fname = '/ltraid2/sleung/data/SeaWiFS/VGPM_1997-2010.nc';
ncid = netcdf.open(fname,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'lon');
nppvlon = netcdf.getVar(ncid,varid,'double');
varid = netcdf.inqVarID(ncid,'lat');
nppvlat = netcdf.getVar(ncid,varid,'double');
nppvtime = datetime('1997-10-01') + calmonths(0:149);
varid = netcdf.inqVarID(ncid,'NPP');
nppv = netcdf.getVar(ncid,varid,'double');
nppv(nppv==-9999)=nan; % missing value = -9999
% rearrange longitude vector to go from 0 to 360
nppvlon1 = nppvlon;
nppvlon1(nppvlon<0)=nppvlon(nppvlon<0)+360;
minidx = find(nppvlon1==min(abs(nppvlon1)));
% rearrange nppv to go from 0 to 360 longitude
nppvlon1 = [nppvlon1(minidx:end); nppvlon1(1:minidx-1)];
nppv1 = nan(size(nppv));
for itime=1:length(nppvtime)
    nppv1(:,:,itime) = [nppv(minidx:end,:,itime); nppv(1:minidx-1,:,itime)];
end
% isolate the tropical Pacific region you want
nppvlonnew = nppvlon1(550:1800);
nppvlatnew = nppvlat(abs(nppvlat)<30);
nppvnew = nppv1(550:1800,find(abs(nppvlat)<30),:);
clear nppvlon nppvlat nppv nppvlon1 nppv1;
nppvlon = nppvlonnew; nppvlat = nppvlatnew; nppv = nppvnew;

[wnppv] = fine2coarsegridbin(nppv,nppvlon,nppvlat,nppvtime,1,1,wcatchlon,wcatchlat,wcatchplon,wcatchplat,1);

tidx_wc_nppvi = find(wtimemocatch==nppvtime(1));
tidx_wc_nppvf = find(wtimemocatch==nppvtime(end));

% - CbPM
fname = '/ltraid2/sleung/data/SeaWiFS/cbpm_1997-2010.nc';
ncid = netcdf.open(fname,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'lon');
nppclon = netcdf.getVar(ncid,varid,'double');
varid = netcdf.inqVarID(ncid,'lat');
nppclat = netcdf.getVar(ncid,varid,'double');
nppctime = datetime('1997-10-01') + calmonths(0:149);
varid = netcdf.inqVarID(ncid,'npp');
nppc = netcdf.getVar(ncid,varid,'double');
nppc(nppc==-9999)=nan; % missing value = -9999
% rearrange longitude vector to go from 0 to 360
nppclon1 = nppclon;
nppclon1(nppclon<0)=nppclon(nppclon<0)+360;
minidx = find(nppclon1==min(abs(nppclon1)));
% rearrange nppc to go from 0 to 360 longitude
nppclon1 = [nppclon1(minidx:end); nppclon1(1:minidx-1)];
nppc1 = nan(size(nppc));
for itime=1:length(nppctime)
    nppc1(:,:,itime) = [nppc(minidx:end,:,itime); nppc(1:minidx-1,:,itime)];
end
% isolate the tropical Pacific region you want
nppclonnew = nppclon1(550:1800);
nppclatnew = nppclat(abs(nppclat)<30);
nppcnew = nppc1(550:1800,find(abs(nppclat)<30),:);
clear nppclon nppclat nppc nppclon1 nppc1;
nppclon = nppclonnew; nppclat = nppclatnew; nppc = nppcnew;

[wnppc] = fine2coarsegridbin(nppc,nppclon,nppclat,nppctime,1,1,wcatchlon,wcatchlat,wcatchplon,wcatchplat,plottocheck);

tidx_wc_nppci = find(wtimemocatch==nppctime(1));
tidx_wc_nppcf = find(wtimemocatch==nppctime(end));

%---------
% WINDS (many variables)
%---------
load('/ltraid2/sleung/data/ccmpwinds/ccmpwind_monthlymeans_tropicalpacificonly.mat');
windtime = datetime('1987-08-01') + calmonths(0:292);

[wwspd] = fine2coarsegridbin(wspd,windlon,windlat,windtime,1/4,1/4,wcatchlon,wcatchlat,wcatchplon,wcatchplat,plottocheck);
[wupstr] = fine2coarsegridbin(upstr,windlon,windlat,windtime,1/4,1/4,wcatchlon,wcatchlat,wcatchplon,wcatchplat,plottocheck);
[wvpstr] = fine2coarsegridbin(vpstr,windlon,windlat,windtime,1/4,1/4,wcatchlon,wcatchlat,wcatchplon,wcatchplat,plottocheck);
[wuwnd] = fine2coarsegridbin(uwnd,windlon,windlat,windtime,1/4,1/4,wcatchlon,wcatchlat,wcatchplon,wcatchplat,plottocheck);
[wvwnd] = fine2coarsegridbin(vwnd,windlon,windlat,windtime,1/4,1/4,wcatchlon,wcatchlat,wcatchplon,wcatchplat,plottocheck);

tidx_wc_windi = find(wtimemocatch==windtime(1));
tidx_wc_windf = find(wtimemocatch==windtime(end));

%----------------------------------------------------
% Put all env vars of interest into a cell array
% for easier plotting/variable naming 
%----------------------------------------------------
wenvvars = cell(11,1); wenvvars{1} = wo2100;
wenvvars{2} = wsst; wenvvars{3} = wssha;
wenvvars{4} = wnppv; wenvvars{5} = wnppc;
wenvvars{6} = wts; wenvvars{7} = wwspd;
wenvvars{8} = wupstr; wenvvars{9} = wvpstr;
wenvvars{10} = wuwnd; wenvvars{11} = wvwnd;

wenvvarstidxs = cell(11,1);
wenvvarstidxs{1} = [1 length(wtimemocatch)];
wenvvarstidxs{2} = [1 length(wtimemocatch)];
wenvvarstidxs{3} = [tidx_wc_sshi tidx_wc_sshf];
wenvvarstidxs{4} = [tidx_wc_nppvi tidx_wc_nppvf];
wenvvarstidxs{5} = [tidx_wc_nppci tidx_wc_nppcf];
wenvvarstidxs{6} = [1 tidx_wc_tsf];
wenvvarstidxs{7} = [tidx_wc_windi tidx_wc_windf];
wenvvarstidxs{8} = [tidx_wc_windi tidx_wc_windf];
wenvvarstidxs{9} = [tidx_wc_windi tidx_wc_windf];
wenvvarstidxs{10} = [tidx_wc_windi tidx_wc_windf];
wenvvarstidxs{11} = [tidx_wc_windi tidx_wc_windf];

wenvvarsnames = cell(11,1);
wenvvarsnames{1} = 'o2100';
wenvvarsnames{2} = 'sst';
wenvvarsnames{3} = 'ssha';
wenvvarsnames{4} = 'nppv';
wenvvarsnames{5} = 'nppc';
wenvvarsnames{6} = 'ts';
wenvvarsnames{7} = 'wspd';
wenvvarsnames{8} = 'upstr';
wenvvarsnames{9} = 'vpstr';
wenvvarsnames{10} = 'uwnd';
wenvvarsnames{11} = 'vwnd';
% 1 = o2100 from WOD
% 2 = sst from WOD
% 3 = ssha from sat
% 4 = nppv from sat
% 5 = nppc from sat
% 6 = sst (called ts to differentiate from WOD) from sat
% 7 = wspd from sat
% 8 = upstr from sat
% 9 = vpstr from sat
% 10 = uwnd from sat
% 11 = vwnd from sat

%----------------------------------------------------
% Save desired vars out into mat files for easier
% repeated loading/use
%----------------------------------------------------
save('o2tempsalfromWODonWCPFCcatchgrid','wo2','wtemp','wsal','wtimemocatch','wcatchlon','wcatchlat','-v7.3');
save('WODandsatenvvarsonWCPFCcatchgrid','wenvvars','wenvvarstidxs','wenvvarsnames','wtimemocatch','wcatchlon','wcatchlat','-v7.3');
