addpath(genpath('/ltraid2/sleung/matlabroutines'));
addpath(genpath('/ltraid2/sleung/tunaextremes/mfiles/seaarounduscatchecon'));
loadiattc = 0;

%----------------------------------------------------
% Load EEZ outlines for later optional plotting on maps
%----------------------------------------------------
% - Vars are: loneez0pt25, lateez0pt25, eezmap0pt25, nupid, etppid, allpid
wlon = 100;elon = 300; slat = -30;nlat = 30;
[s,a] = shaperead('/ltraid2/sleung/data/EEZs/World_EEZ_v9_20161021_HR_0_360/World_EEZ_v9_2016_HR_0_360','UseGeoCoords',true,'BoundingBox',[wlon,slat;elon,nlat]);
c = struct2cell(a); allctry = c(10,:); allpid = c(1,:);

nuctries = {'Micronesia','Kiribati','Marshall Islands','Nauru','Palau','Papua New Guinea','Solomon Islands','Tuvalu'};
nupid = cell(length(nuctries),1);
for ictry = 1:length(nuctries)
    cidx = find(strcmp(allctry,nuctries{ictry}));
    nupid{ictry} = allpid(cidx);
end

etpctries = {'Peru','Ecuador','Colombia','Panama','Costa Rica','Nicaragua','El Salvador','Guatemala'};
etppid = cell(length(etpctries),1);
for ictry = 1:length(etpctries)
    cidx = find(strcmp(allctry,etpctries{ictry}));
    etppid{ictry} = allpid(cidx);
end

fnamenc = '/ltraid2/sleung/data/EEZs/World_EEZ_v9_20161021_HR_0_360/World_EEZ_v9_2016_HR_0_360_30Nto30S_polyid_0pt25deg.nc';
ncid=netcdf.open(fnamenc,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'z'); % dims: lon, lat
eezmap0pt25 = netcdf.getVar(ncid,varid,'double');,
eezmap0pt25 = eezmap0pt25'; % x = lon, y = lat
netcdf.close(ncid);
numdegreseez = 0.25;
loneez0pt25 = 0:numdegreseez:360;
lateez0pt25 = (-30:numdegreseez:30)';

%----------------------------------------------------
% Define plot parameters for maps
%----------------------------------------------------
ylgnbunow = cbrewer('seq','YlGnBu',11,'linear');
rdylbunow = flipud(cbrewer('div','PuOr',15,'linear'));
rdylbunow(8,:) = [1.0000    1.0000    0.7490];
eezltlinecolor = [0.6 0.6 0.6];
eezdklinecolor = [0.3 0.3 0.3];
eezltlinewidth = 1.5;
eezdklinewidth = 1.5;
mapproj = 'gall-peters';
coastlinewidth = 1;
landcolor = [0.6 0.6 0.6];
wlon = 100;elon = 285;
lonticks = [120 160 200 240 280];
slat = -22.5;nlat = 22.5;
latticks = [-20 -10 0 10 20];
showeezs = 0;

%----------------------------------------------------
% Load monthly enso time series
%----------------------------------------------------
addpath(genpath('/ltraid2/sleung/tunaextremes/mfiles/obsanalysis/scripts/wodenso_anyres_scripts/load'));
plotensoidx=0;

begdate = datetime(1967,1,1);
enddate = datetime(2017,12,1);
wonitime = datetime('1967-01-01') + calmonths(0:611);
loadensoindices;
woni = oni; wonien = onien; woniln = oniln;

if loadiattc==1
begdate = datetime(1958,12,1);
enddate = datetime(2018,12,1);
eonitime = datetime('1958-12-01') + calmonths(0:720);
loadensoindices;
eoni = oni; eonien = onien; eoniln = oniln;
end

%----------------------------------------------------
% Load monthly 5 by 5 degree WCPFC catch and effort data
%----------------------------------------------------
ncfile='/ltraid2/sleung/data/WCPFC/5by5deg/ByMonth/PurseSeine_0/WCPFC_purseseine_bysettype_monthly_5deg.nc';
%ncinfo(ncfile)
%ncdisp(ncfile)
ncid = netcdf.open(ncfile,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'lon');
wcatchlon = netcdf.getVar(ncid,varid,'double');
varid = netcdf.inqVarID(ncid,'lat');
wcatchlat = netcdf.getVar(ncid,varid,'double');
wcatchplat = wcatchlat-5/2;
wcatchplon = wcatchlon-5/2;
wtimemocatch = datetime('1967-01-01') + calmonths(0:611);

varid = netcdf.inqVarID(ncid,'skj_c_una');
skj_c_una = netcdf.getVar(ncid,varid,'double');
skj_c_una(skj_c_una<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'skj_c_log');
skj_c_log = netcdf.getVar(ncid,varid,'double');
skj_c_log(skj_c_log<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'skj_c_dfad');
skj_c_dfad = netcdf.getVar(ncid,varid,'double');
skj_c_dfad(skj_c_dfad<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'skj_c_afad');
skj_c_afad = netcdf.getVar(ncid,varid,'double');
skj_c_afad(skj_c_afad<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'skj_c_oth');
skj_c_oth = netcdf.getVar(ncid,varid,'double');
skj_c_oth(skj_c_oth<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_una');
bet_c_una = netcdf.getVar(ncid,varid,'double');
bet_c_una(bet_c_una<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_log');
bet_c_log = netcdf.getVar(ncid,varid,'double');
bet_c_log(bet_c_log<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_dfad');
bet_c_dfad = netcdf.getVar(ncid,varid,'double');
bet_c_dfad(bet_c_dfad<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_afad');
bet_c_afad = netcdf.getVar(ncid,varid,'double');
bet_c_afad(bet_c_afad<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_oth');
bet_c_oth = netcdf.getVar(ncid,varid,'double');
bet_c_oth(bet_c_oth<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_una');
sets_una = netcdf.getVar(ncid,varid,'double');
sets_una(sets_una<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_log');
sets_log = netcdf.getVar(ncid,varid,'double');
sets_log(sets_log<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_dfad');
sets_dfad = netcdf.getVar(ncid,varid,'double');
sets_dfad(sets_dfad<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_afad');
sets_afad = netcdf.getVar(ncid,varid,'double');
sets_afad(sets_afad<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_oth');
sets_oth = netcdf.getVar(ncid,varid,'double');
sets_oth(sets_oth<-1e33)=nan;

netcdf.close(ncid);

%----------------------------------------------------
% Calculate WCPFC CPUE
%----------------------------------------------------
skj_cp_una = skj_c_una./sets_una;
skj_cp_log = skj_c_log./sets_log;
skj_cp_dfad = skj_c_dfad./sets_dfad;
skj_cp_afad = skj_c_afad./sets_afad;
skj_cp_oth = skj_c_oth./sets_oth;
bet_cp_una = bet_c_una./sets_una;
bet_cp_log = bet_c_log./sets_log;
bet_cp_dfad = bet_c_dfad./sets_dfad;
bet_cp_afad = bet_c_afad./sets_afad;
bet_cp_oth = bet_c_oth./sets_oth;

if loadiattc==1
%----------------------------------------------------
% Load monthly 1 by 1 degree IATTC catch data
%----------------------------------------------------
ncfile='/ltraid2/sleung/data/IATTC/PublicPSTuna/IATTC_purseseine_bysettype_monthly_1deg.nc';
%ncinfo(ncfile)
%ncdisp(ncfile)
ncid = netcdf.open(ncfile,'nc_nowrite');
varid = netcdf.inqVarID(ncid,'lon');
ecatchlon = netcdf.getVar(ncid,varid,'double');
varid = netcdf.inqVarID(ncid,'lat');
ecatchlat = netcdf.getVar(ncid,varid,'double');
ecatchplat = ecatchlat-1/2;
ecatchplon = ecatchlon-1/2;
etimemocatch = datetime('1958-12-01') + calmonths(0:720);

varid = netcdf.inqVarID(ncid,'skj_c_noa');
skj_c_noa = netcdf.getVar(ncid,varid,'double');
skj_c_noa(skj_c_noa<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'skj_c_del');
skj_c_del = netcdf.getVar(ncid,varid,'double');
skj_c_del(skj_c_del<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'skj_c_obj');
skj_c_obj = netcdf.getVar(ncid,varid,'double');
skj_c_obj(skj_c_obj<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_noa');
bet_c_noa = netcdf.getVar(ncid,varid,'double');
bet_c_noa(bet_c_noa<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_del');
bet_c_del = netcdf.getVar(ncid,varid,'double');
bet_c_del(bet_c_del<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'bet_c_obj');
bet_c_obj = netcdf.getVar(ncid,varid,'double');
bet_c_obj(bet_c_obj<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_noa');
sets_noa = netcdf.getVar(ncid,varid,'double');
sets_noa(sets_noa<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_del');
sets_del = netcdf.getVar(ncid,varid,'double');
sets_del(sets_del<-1e33)=nan;

varid = netcdf.inqVarID(ncid,'sets_obj');
sets_obj = netcdf.getVar(ncid,varid,'double');
sets_obj(sets_obj<-1e33)=nan;

netcdf.close(ncid);

%----------------------------------------------------
% Calculate IATTC CPUE
%----------------------------------------------------
skj_cp_noa = skj_c_noa./sets_noa;
skj_cp_del = skj_c_del./sets_del;
skj_cp_obj = skj_c_obj./sets_obj;
bet_cp_noa = bet_c_noa./sets_noa;
bet_cp_del = bet_c_del./sets_del;
bet_cp_obj = bet_c_obj./sets_obj;
end

%----------------------------------------------------
% Calculate WCPFC AND IATTC total catch and CPUE
%----------------------------------------------------
wskj_c_tot = skj_c_una + skj_c_log + skj_c_dfad + skj_c_afad + skj_c_oth;
wbet_c_tot = bet_c_una + bet_c_log + bet_c_dfad + bet_c_afad + bet_c_oth;
wskj_cp_tot = skj_cp_una + skj_cp_log + skj_cp_dfad + skj_cp_afad + skj_cp_oth;
wbet_cp_tot = bet_cp_una + bet_cp_log + bet_cp_dfad + bet_cp_afad + bet_cp_oth;

if loadiattc==1
eskj_c_tot = skj_c_noa + skj_c_del + skj_c_obj;
ebet_c_tot = bet_c_noa + bet_c_del + bet_c_obj;
eskj_cp_tot = skj_cp_noa + skj_cp_del + skj_cp_obj;
ebet_cp_tot = bet_cp_noa + bet_cp_del + bet_cp_obj;
end

%----------------------------------------------------
% Combine data from WCPFC and IATTC
%----------------------------------------------------
% UGH there's an overlap area!!! from 150W-130W, 50S-4S,
% aka 210-230, -50--4

% don't combine them; it's a trap; can't hold
% diff resolution data in one matrix

%----------------------------------------------------
% Calculate mean and enso composites of total catch,
% CPUE, effort, and bet:skj catch ratio
%----------------------------------------------------
wskj_c_tot_mean = nanmean(wskj_c_tot,3);
wbet_c_tot_mean = nanmean(wbet_c_tot,3);
wskj_c_tot_en = nanmean(wskj_c_tot(:,:,logical(wonien)),3);
wbet_c_tot_en = nanmean(wbet_c_tot(:,:,logical(wonien)),3);
wskj_c_tot_ln = nanmean(wskj_c_tot(:,:,logical(woniln)),3);
wbet_c_tot_ln = nanmean(wbet_c_tot(:,:,logical(woniln)),3);

wbettoskj_c_tot = wbet_c_tot./wskj_c_tot;
wbettoskj_c_tot_mean = nanmean(wbettoskj_c_tot,3);
wbettoskj_c_tot_en = nanmean(wbettoskj_c_tot(:,:,logical(wonien)),3);
wbettoskj_c_tot_ln = nanmean(wbettoskj_c_tot(:,:,logical(woniln)),3);

if loadiattc==1
eskj_c_tot_mean = nanmean(eskj_c_tot,3);
ebet_c_tot_mean = nanmean(ebet_c_tot,3);
eskj_c_tot_en = nanmean(eskj_c_tot(:,:,logical(eonien)),3);
ebet_c_tot_en = nanmean(ebet_c_tot(:,:,logical(eonien)),3);
eskj_c_tot_ln = nanmean(eskj_c_tot(:,:,logical(eoniln)),3);
ebet_c_tot_ln = nanmean(ebet_c_tot(:,:,logical(eoniln)),3);

ebettoskj_c_tot = ebet_c_tot./eskj_c_tot;
ebettoskj_c_tot_mean = nanmean(ebettoskj_c_tot,3);
ebettoskj_c_tot_en = nanmean(ebettoskj_c_tot(:,:,logical(eonien)),3);
ebettoskj_c_tot_ln = nanmean(ebettoskj_c_tot(:,:,logical(eoniln)),3);
end

%----------------------------------------------------
% PLOT MEAN ENSO COMPOSITE SKJ AND BET CATCH
%----------------------------------------------------
wplotmeancatch=0; % w for west or WCPFC
if wplotmeancatch==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 18 7],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
climsdefined = 0; cminmax = [0 5000];
unitsnow = '[metric tons]';
axnow=subplot(3,3,1); colormapnow = ylgnbunow;
mapvarnow = wbet_c_tot_mean'; titlenow = 'bet mean'; plotmaps;
axnow=subplot(3,3,2); colormapnow = ylgnbunow;
mapvarnow = wskj_c_tot_mean'; titlenow = 'skj mean'; plotmaps;
unitsnow = '[unitless]';
axnow=subplot(3,3,3); colormapnow = ylgnbunow; 
mapvarnow = wbettoskj_c_tot_mean'; titlenow = 'bet:skj mean'; plotmaps;
unitsnow = '[metric tons]';
axnow=subplot(3,3,4); colormapnow = ylgnbunow;
mapvarnow = wbet_c_tot_en'; titlenow = 'bet en'; plotmaps;
axnow=subplot(3,3,5); colormapnow = ylgnbunow;
mapvarnow = wskj_c_tot_en'; titlenow = 'skj en'; plotmaps;
unitsnow = '[unitless]';
axnow=subplot(3,3,6); colormapnow = ylgnbunow;
mapvarnow = wbettoskj_c_tot_en'; titlenow = 'bet:skj en'; plotmaps;
unitsnow = '[metric tons]';
axnow=subplot(3,3,7); colormapnow = ylgnbunow;
mapvarnow = wbet_c_tot_ln'; titlenow = 'bet ln'; plotmaps;
axnow=subplot(3,3,8); colormapnow = ylgnbunow;
mapvarnow = wskj_c_tot_ln'; titlenow = 'skj ln'; plotmaps;
unitsnow = '[unitless]';
axnow=subplot(3,3,9); colormapnow = ylgnbunow;
mapvarnow = wbettoskj_c_tot_ln'; titlenow = 'bet:skj ln'; plotmaps;
%print('WCPFC_skjbetmeanenlncatch','-dpng');
print('WCPFC_betskjmeanenlncatch_ratios','-dpng');
end

%----------------------------------------------------
% Calculate correlations btwn monthly skj and bet catch
%----------------------------------------------------
[wbetvsskj_c_tot_cc,wbetvsskj_c_tot_rc,~,~,~,~,wbetvsskj_c_tot_ccpval] = tempcorrmapnanwithintandslopeuncertainty(wskj_c_tot,wbet_c_tot);
[wbetvsskj_c_toten_cc,wbetvsskj_c_toten_rc,~,~,~,~,wbetvsskj_c_toten_ccpval] = tempcorrmapnanwithintandslopeuncertainty(wskj_c_tot(:,:,logical(wonien)),wbet_c_tot(:,:,logical(wonien)));
[wbetvsskj_c_totln_cc,wbetvsskj_c_totln_rc,~,~,~,~,wbetvsskj_c_totln_ccpval] = tempcorrmapnanwithintandslopeuncertainty(wskj_c_tot(:,:,logical(woniln)),wbet_c_tot(:,:,logical(woniln)));

%----------------------------------------------------
% Calculate correlations btwn monthly skj catch and bet:skj catch ratio
%----------------------------------------------------
[wbettoskjvsskj_c_tot_cc,wbettoskjvsskj_c_tot_rc,~,~,~,~,wbettoskjvsskj_c_tot_ccpval] = tempcorrmapnanwithintandslopeuncertainty(wskj_c_tot,wbettoskj_c_tot);
[wbettoskjvsskj_c_toten_cc,wbettoskjvsskj_c_toten_rc,~,~,~,~,wbettoskjvsskj_c_toten_ccpval] = tempcorrmapnanwithintandslopeuncertainty(wskj_c_tot(:,:,logical(wonien)),wbettoskj_c_tot(:,:,logical(wonien)));
[wbettoskjvsskj_c_totln_cc,wbettoskjvsskj_c_totln_rc,~,~,~,~,wbettoskjvsskj_c_totln_ccpval] = tempcorrmapnanwithintandslopeuncertainty(wskj_c_tot(:,:,logical(woniln)),wbettoskj_c_tot(:,:,logical(woniln)));

%----------------------------------------------------
% PLOT CORRELATIONS BTWN MONTHLY SKJ AND BET CATCH OR BET:SKJ
%----------------------------------------------------
wplotcatch=0; % w for west or WCPFC
if wplotcatch==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 6 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat; unitsnow = '';
axnow=subplot(3,1,1); colormapnow = rdylbunow; mapvarnow = wbetvsskj_c_tot_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET vs. SKJ tot catch cc'; plotmaps;
axnow=subplot(3,1,2); colormapnow = rdylbunow; mapvarnow = wbetvsskj_c_tot_rc';
climsdefined = 1; cminmax = [-0.17 0.17]; titlenow = 'WCPFC BET vs. SKJ tot catch rc'; plotmaps;
axnow=subplot(3,1,3); colormapnow = rdylbunow; mapvarnow = wbetvsskj_c_tot_ccpval';
climsdefined = 1; cminmax = [0 0.05]; titlenow = 'WCPFC BET vs. SKJ tot catch cc pval'; plotmaps;
%print('WCPFC_betvsskjtotcatch_ccrcpval','-dpng');
print('WCPFC_betvsskjtotcatch_ccrcpval_allmonths','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 6 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat; unitsnow = '';
axnow=subplot(3,1,1); colormapnow = rdylbunow; mapvarnow = wbetvsskj_c_tot_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET vs. SKJ tot catch cc - all months'; plotmaps;
axnow=subplot(3,1,2); colormapnow = rdylbunow; mapvarnow = wbetvsskj_c_toten_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET vs. SKJ tot catch cc - EN months only'; plotmaps;
axnow=subplot(3,1,3); colormapnow = rdylbunow; mapvarnow = wbetvsskj_c_totln_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET vs. SKJ tot catch cc - LN months only'; plotmaps;
print('WCPFC_betvsskjtotcatch_cc_allENLNmonths','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 6 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat; unitsnow = '';
axnow=subplot(3,1,1); colormapnow = rdylbunow; mapvarnow = wbettoskjvsskj_c_tot_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET:SKJ vs. SKJ tot catch cc'; plotmaps;
axnow=subplot(3,1,2); colormapnow = rdylbunow; mapvarnow = wbettoskjvsskj_c_tot_rc';
climsdefined = 1; cminmax = [-1e-5 1e-5]; titlenow = 'WCPFC BET:SKJ vs. SKJ tot catch rc'; plotmaps;
axnow=subplot(3,1,3); colormapnow = rdylbunow; mapvarnow = wbettoskjvsskj_c_tot_ccpval';
climsdefined = 1; cminmax = [0 0.05]; titlenow = 'WCPFC BET:SKJ vs. SKJ tot catch cc pval'; plotmaps;
print('WCPFC_bettoskjvsskjtotcatch_ccrcpval_allmonths','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 6 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat; unitsnow = '';
axnow=subplot(3,1,1); colormapnow = rdylbunow; mapvarnow = wbettoskjvsskj_c_tot_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET:SKJ vs. SKJ tot catch cc - all months'; plotmaps;
axnow=subplot(3,1,2); colormapnow = rdylbunow; mapvarnow = wbettoskjvsskj_c_toten_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET:SKJ vs. SKJ tot catch cc - EN months only'; plotmaps;
axnow=subplot(3,1,3); colormapnow = rdylbunow; mapvarnow = wbettoskjvsskj_c_totln_cc';
climsdefined = 1; cminmax = [-1 1]; titlenow = 'WCPFC BET:SKJ vs. SKJ tot catch cc - LN months only'; plotmaps;
print('WCPFC_bettoskjvsskjtotcatch_cc_allENLNmonths','-dpng');
end

%----------------------------------------------------
% Calclate climatologies/seasons of spatial catch distributions
%----------------------------------------------------
wskj_c_tot_clim = nan(length(wcatchlon),length(wcatchlat),12);
wbet_c_tot_clim = nan(length(wcatchlon),length(wcatchlat),12);
% this assumes that skj_c_tot, bet_c_tot start in jan
for imonth = 1:12
    wskj_c_tot_clim(:,:,imonth) = nanmean(wskj_c_tot(:,:,imonth:12:length(wtimemocatch)),3);
    wbet_c_tot_clim(:,:,imonth) = nanmean(wbet_c_tot(:,:,imonth:12:length(wtimemocatch)),3);
end

% season 1 = winter, 2 = spr, 3 = summ, 4 = fall
wskj_c_tot_seas = nan(length(wcatchlon),length(wcatchlat),4);
wbet_c_tot_seas = nan(length(wcatchlon),length(wcatchlat),4);
wskj_c_tot_seas(:,:,1) = nanmean(wskj_c_tot_clim(:,:,[12 1 2]),3); 
wbet_c_tot_seas(:,:,1) = nanmean(wbet_c_tot_clim(:,:,[12 1 2]),3); 
iseas = 2;
for imonth = [3 6 9]
    wskj_c_tot_seas(:,:,iseas) = nanmean(wskj_c_tot_clim(:,:,imonth:imonth+2),3);
    wbet_c_tot_seas(:,:,iseas) = nanmean(wbet_c_tot_clim(:,:,imonth:imonth+2),3);
    iseas = iseas + 1;
end

%----------------------------------------------------
% PLOT SEASONS OF SPATIAL CATCH DISTRIBUTIONS
%----------------------------------------------------
wplotcatch=0; % w for west or WCPFC
if wplotcatch==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 12 12],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat; unitsnow = '[metric tons]';
colormapnow = ylgnbunow; 
isp = 1;
for iseas = 1:4
    axnow = subplot(4,2,isp); mapvarnow = wskj_c_tot_seas(:,:,iseas)';   
    climsdefined = 1; cminmax = [0 4000];
    titlenow = ['WCPFC SKJ mean seas #' num2str(iseas)]; plotmaps; 
    isp = isp+2;
end
isp = 2;
for iseas = 1:4
    axnow = subplot(4,2,isp); mapvarnow = wbet_c_tot_seas(:,:,iseas)';
    climsdefined = 1; cminmax = [0 300];
    titlenow = ['WCPFC BET mean seas #' num2str(iseas)]; plotmaps;
    isp = isp+2;
end
print('WCPFC_skjbetcatch_diffseasons','-dpng');
end

%----------------------------------------------------
% PLOT CLIMATOLOGIES OF SPATIAL CATCH DISTRIBUTIONS
%----------------------------------------------------
wplotcatch=0; % w for west or WCPFC
if wplotcatch==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 18 12],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat; unitsnow = '[metric tons]';
colormapnow = ylgnbunow;
for imonth = 1:12
    axnow = subplot(4,3,imonth); mapvarnow = wskj_c_tot_clim(:,:,imonth)';
    climsdefined = 1; cminmax = [0 5000];
    titlenow = ['WCPFC SKJ mean month #' num2str(imonth)]; plotmaps;
end
print('WCPFC_skjcatch_climatology','-dpng');
f=figure;
set(f,'color','white','units','inches','position',[0 0 18 12],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat; unitsnow = '[metric tons]';
colormapnow = ylgnbunow;
for imonth = 1:12
    axnow = subplot(4,3,imonth); mapvarnow = wbet_c_tot_clim(:,:,imonth)';
    climsdefined = 1; cminmax = [0 400];
    titlenow = ['WCPFC BET mean month #' num2str(imonth)]; plotmaps;
end
print('WCPFC_betcatch_climatology','-dpng');
end

%----------------------------------------------------
% Load environmental variables from WOD and satellite
% on WCPFC catch grid
%----------------------------------------------------
%loadsaveregridenvvarsfromWODandsatontoWCPFCgrid;
load('o2tempsalfromWODonWCPFCcatchgrid','wo2','wtemp','wsal');
load('WODandsatenvvarsonWCPFCcatchgrid','wenvvars','wenvvarstidxs','wenvvarsnames','wenvvarsunits'); 
%wenvvarsunits = {'umol/kg','degC','m','mgC/m^2/day','mgC/m^2/day','degC','m/s','m^2/s^2','m^2/s^2','m/s','m/s'};
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
% Calculate time vars for plotting
%----------------------------------------------------
wmocatch3d = nan(size(wskj_c_tot)); 
for itime = 1:length(wtimemocatch)
    wmocatch3d(:,:,itime) = month(wtimemocatch(itime));
end

wseascatch3d = nan(size(wskj_c_tot));
for itime = 1:length(wtimemocatch)
    if month(wtimemocatch(itime))==12 | month(wtimemocatch(itime))<=2
        wseascatch3d(:,:,itime) = 1;
    elseif month(wtimemocatch(itime))>=3 & month(wtimemocatch(itime))<=5
        wseascatch3d(:,:,itime) = 2;
    elseif month(wtimemocatch(itime))>=6 & month(wtimemocatch(itime))<=8
        wseascatch3d(:,:,itime) = 3;
    elseif month(wtimemocatch(itime))>=9 & month(wtimemocatch(itime))<=11
        wseascatch3d(:,:,itime) = 4;
    end
end

woni3d = nan(size(wskj_c_tot));
for itime = 1:length(wtimemocatch)
    woni3d(:,:,itime) = woni(itime);
end

%----------------------------------------------------
% SCATTERPLOT BTWN CATCH AND ENV VARS
%----------------------------------------------------
plotscatter=0;
if plotscatter==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 8 12],'resize','on');
for isp = 1:9
    subplot(3,3,isp);
    wseascatch3dnow = wseascatch3d(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    wskj_c_totnow = wskj_c_tot(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    scatter(wenvvars{isp}(:),wskj_c_totnow(:),10,wseascatch3dnow(:));
    title(['skj vs. ' wenvvarsnames{isp}]);
end
%print('WCPFC_skjcatchvsenvvars_scatterplot','-dpng');
print('WCPFC_skjcatchvsenvvars_scatterplot_coloredbyseason','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 14 10],'resize','on');
for isp = 1:9
    subplot(3,3,isp);
    wseascatch3dnow = wseascatch3d(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    wbet_c_totnow = wbet_c_tot(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    scatter(wenvvars{isp}(:),wbet_c_totnow(:),10,wseascatch3dnow(:));
    title(['bet vs. ' wenvvarsnames{isp}]);
end
%print('WCPFC_betcatchvsenvvars_scatterplot','-dpng');
print('WCPFC_betcatchvsenvvars_scatterplot_coloredbyseason','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 8 12],'resize','on');
colormap(rdylbunow);
for isp = 1:9
    subplot(3,3,isp);
    woni3dnow = woni3d(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    wskj_c_totnow = wskj_c_tot(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    scatter(wenvvars{isp}(:),wskj_c_totnow(:),10,woni3dnow(:));
    title(['skj vs. ' wenvvarsnames{isp}]);
end
print('WCPFC_skjcatchvsenvvars_scatterplot_coloredbyoni','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 14 10],'resize','on');
colormap(rdylbunow);
for isp = 1:9
    subplot(3,3,isp);
    woni3dnow = woni3d(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    wbet_c_totnow = wbet_c_tot(:,:,wenvvarstidxs{isp}(1):wenvvarstidxs{isp}(2));
    scatter(wenvvars{isp}(:),wbet_c_totnow(:),10,woni3dnow(:));
    title(['bet vs. ' wenvvarsnames{isp}]);
end
print('WCPFC_betcatchvsenvvars_scatterplot_coloredbyoni','-dpng');
end

plotscatter=0;
if plotscatter==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 8 8],'resize','on');
isp=1;
for ivar = [1 3 4 5 6 7]
    subplot(3,2,isp);
    wseascatch3dnow = wseascatch3d(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    wskj_c_totnow = wskj_c_tot(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    scatter(wenvvars{ivar}(:),wskj_c_totnow(:),10,wseascatch3dnow(:));
    ylabel('Catch [metric tons]'); xlabel(wenvvarsunits{ivar});
    title(['SKJ catch vs. ' wenvvarsnames{ivar}]);
    isp=isp+1;
end
%print('WCPFC_skjcatchvsenvvars_scatterplot_reduced','-dpng');
print('WCPFC_skjcatchvsenvvars_scatterplot_reduced_coloredbyseason','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 8 8],'resize','on');
isp=1;
for ivar = [1 3 4 5 6 7]
    subplot(3,2,isp);
    wseascatch3dnow = wseascatch3d(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    wbet_c_totnow = wbet_c_tot(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));  
    scatter(wenvvars{ivar}(:),wbet_c_totnow(:),10,wseascatch3dnow(:));
    ylabel('Catch [metric tons]'); xlabel(wenvvarsunits{ivar});
    title(['BET catch vs. ' wenvvarsnames{ivar}]);
    isp=isp+1;
end
%print('WCPFC_betcatchvsenvvars_scatterplot_reduced','-dpng');
print('WCPFC_betcatchvsenvvars_scatterplot_reduced_coloredbyseason','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 8 8],'resize','on');
isp=1;
for ivar = [1 3 4 5 6 7]
    subplot(3,2,isp);
    woni3dnow = woni3d(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    wskj_c_totnow = wskj_c_tot(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    scatter(wenvvars{ivar}(:),wskj_c_totnow(:),10,woni3dnow(:));
    ylabel('Catch [metric tons]'); xlabel(wenvvarsunits{ivar});
    colormap(rdylbunow); caxis([-2 2]); %colorbar;
    title(['SKJ catch vs. ' wenvvarsnames{ivar}]);
    isp=isp+1;
end
print('WCPFC_skjcatchvsenvvars_scatterplot_reduced_coloredbyoni','-dpng');
%print('WCPFC_skjcatchvsenvvars_scatterplot_reduced_coloredbyoni_colorbar','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 8 8],'resize','on');
isp=1;
for ivar = [1 3 4 5 6 7]
    subplot(3,2,isp); 
    woni3dnow = woni3d(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    wbet_c_totnow = wbet_c_tot(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    scatter(wenvvars{ivar}(:),wbet_c_totnow(:),10,woni3dnow(:));
    ylabel('Catch [metric tons]'); xlabel(wenvvarsunits{ivar});
    colormap(rdylbunow); caxis([-2 2]);
    title(['BET catch vs. ' wenvvarsnames{ivar}]);
    isp=isp+1;
end 
print('WCPFC_betcatchvsenvvars_scatterplot_reduced_coloredbyoni','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 8 8],'resize','on');
isp=1;
for ivar = [1 3 4 5 6 7]
    subplot(3,2,isp);
    woni3dnow = woni3d(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    wbettoskj_c_totnow = wbettoskj_c_tot(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2));
    scatter(wenvvars{ivar}(:),wbettoskj_c_totnow(:),10,woni3dnow(:));
    ylabel('BET:SKJ catch'); xlabel(wenvvarsunits{ivar});
    colormap(rdylbunow); caxis([-2 2]);
    title(['BET:SKJ vs. ' wenvvarsnames{ivar}]);
    isp=isp+1;
end
print('WCPFC_catchratiovsenvvars_scatterplot_reduced_coloredbyoni','-dpng');
end

%----------------------------------------------------
% Calculate correlations btwn monthly skj/bet catch
% and env vars from WOD and satellite
%----------------------------------------------------
for ivar = 1:length(wenvvars)
    [wskj_c_totvswenvvars_cc{ivar},wskj_c_totvswenvvars_rc{ivar},~,~,~,wskj_c_totvswenvvars_numpts{ivar},wskj_c_totvswenvvars_ccpval{ivar}] = tempcorrmapnanwithintandslopeuncertainty(wenvvars{ivar},wskj_c_tot(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2))); 
end

for ivar = 1:length(wenvvars)
    [wbet_c_totvswenvvars_cc{ivar},wbet_c_totvswenvvars_rc{ivar},~,~,~,wbet_c_totvswenvvars_numpts{ivar},wbet_c_totvswenvvars_ccpval{ivar}] = tempcorrmapnanwithintandslopeuncertainty(wenvvars{ivar},wbet_c_tot(:,:,wenvvarstidxs{ivar}(1):wenvvarstidxs{ivar}(2))); 
end

%----------------------------------------------------
% CORREL COEFF MAPS BTWN CATCH AND ENV VARS
%----------------------------------------------------
plotcorrel=0;
if plotcorrel==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 1; cminmax = [-1 1];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wskj_c_totvswenvvars_cc{isp}'; 
    unitsnow = ''; titlenow = ['cc - skj vs. ' wenvvarsnames{isp}]; plotmaps;
end
print('WCPFC_skjcatchvsenvvars_correlcoeff','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 1; cminmax = [-1 1];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wbet_c_totvswenvvars_cc{isp}'; 
    unitsnow = ''; titlenow = ['cc - bet vs. ' wenvvarsnames{isp}]; plotmaps;
end
print('WCPFC_betcatchvsenvvars_correlcoeff','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 0; cminmax = [-1 1];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wskj_c_totvswenvvars_rc{isp}';
    unitsnow = ''; titlenow = ['rc - skj vs. ' wenvvarsnames{isp}]; plotmaps;
end
print('WCPFC_skjcatchvsenvvars_regresscoeff','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 0; cminmax = [-1 1];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wbet_c_totvswenvvars_rc{isp}';
    unitsnow = ''; titlenow = ['rc - bet vs. ' wenvvarsnames{isp}]; plotmaps;
end
print('WCPFC_betcatchvsenvvars_regresscoeff','-dpng');

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 1; cminmax = [0 0.05];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wskj_c_totvswenvvars_ccpval{isp}';
    unitsnow = ''; titlenow = ['ccpval - skj vs. ' wenvvarsnames{isp}]; plotmaps;
end

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 1; cminmax = [0 0.05];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wbet_c_totvswenvvars_ccpval{isp}';
    unitsnow = ''; titlenow = ['ccpval - bet vs. ' wenvvarsnames{isp}]; plotmaps;
end

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 0; cminmax = [0 0.05];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wskj_c_totvswenvvars_numpts{isp}';
    unitsnow = ''; titlenow = ['numpts - skj vs. ' wenvvarsnames{isp}]; plotmaps;
end

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 0; cminmax = [0 0.05];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = wbet_c_totvswenvvars_numpts{isp}';
    unitsnow = ''; titlenow = ['numpts - bet vs. ' wenvvarsnames{isp}]; plotmaps;
end

end

%----------------------------------------------------
% VALIDATING SST VERSUS TS
%----------------------------------------------------
plotval=0;

if plotval==1
f=figure;
set(f,'color','white','units','inches','position',[0 0 5 5],'resize','on');
sstnow = wenvvars{2}(:,:,wenvvarstidxs{6}(1):wenvvarstidxs{6}(2));
tsnow = wenvvars{6};
scatter(sstnow(:),tsnow(:)); hold on;
plot([-5 35],[-5 35]);
xlabel('In situ SST (degC)');
ylabel('Satellite SST (degC)');

f=figure;
set(f,'color','white','units','inches','position',[0 0 4 4],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = '[degC]'; colormapnow = rdylbunow;
climsdefined = 1; cminmax = [20 30];
axnow = subplot(2,1,1); mapvarnow = nanmean(wenvvars{2},3)';
titlenow = 'Mean in situ SST'; plotmaps;
axnow = subplot(2,1,2); mapvarnow = nanmean(wenvvars{6},3)';
titlenow = 'Mean satellite SST'; plotmaps;

f=figure;
set(f,'color','white','units','inches','position',[0 0 16 8],'resize','off');
plonnow = wcatchplon; platnow = wcatchplat;
unitsnow = ''; colormapnow = rdylbunow;
climsdefined = 0; cminmax = [0 0.05];
for isp = 1:9
    axnow = subplot(3,3,isp); mapvarnow = nanmean(wenvvars{isp},3)';
    unitsnow = ''; titlenow = ['mean - ' wenvvarsnames{isp}]; plotmaps;
end
end

%----------------------------------------------------
% Quotient Analysis
%----------------------------------------------------
disp('pause before quotient analysis')
pause

ts = wenvvars{6};
figure; hist(ts(:));
maxts = max(max(max(ts))); %31.6571
mints = min(min(min(ts))); %-1.0862

wskj_c_totnow = wskj_c_tot(:,:,wenvvarstidxs{6}(1):wenvvarstidxs{6}(2));
tsnow = nan(size(ts));
tsnow(~isnan(wskj_c_totnow)) = ts(~isnan(wskj_c_totnow));
maxtsnow = max(max(max(tsnow))); %31.0324
mintsnow = min(min(min(tsnow))); %1.9568

figure; hist(tsnow(:),20);
[n,x]=hist(tsnow(:),20);
100*n/sum(n)
figure; hist(tsnow(:),50);
[n,x]=hist(tsnow(:),50);
100*n/sum(n)

x = [0:5:20 25:0.5:30]; 
n = histc(tsnow(:),x);
figure; bar(x,n,'histc');
100*n/sum(n)

x = 1.5:0.5:31.5;
n = histc(tsnow(:),x);
figure; bar(x,n,'histc');
100*n/sum(n)

wskj_c_totnowtot = nansum(nansum(nansum(wskj_c_totnow))); 
wskj_c_tot_freq_tsbins = nan(length(n)-1,1); 
for i = 1:(length(n)-1)
    wskj_c_tot_freq_tsbins(i) = ...
        nansum(wskj_c_totnow(tsnow>x(i)&tsnow<x(i+1)))...
        ./wskj_c_totnowtot;
end
wskj_c_tot_ts_qa = wskj_c_tot_freq_tsbins ./ (n(1:end-1)/sum(n(1:end-1)));

figure; ax = gca;
yyaxis left;
bar(x,n,'histc');
yyaxis right;
plot(x(1:end-1)+0.25,wskj_c_tot_ts_qa,'color','r');hold on;
plot(x(1:end-1)+0.25,ones(length(x)-1,1),'color','r');
ax.YColor = 'r';

% without replacement
nruns = 399;
wskj_c_tot_ts_bsqa = nan(length(n)-1,nruns); 
for irun = 1:nruns
    X = randperm(numel(wskj_c_totnow));
    wskj_c_totnowshuffled = reshape(wskj_c_totnow(X),size(wskj_c_totnow)); 

    wskj_c_tot_freq_tsbinsnow = nan(length(n)-1,1);
    for i = 1:(length(n)-1)
        wskj_c_tot_freq_tsbinsnow(i) = ...
            nansum(wskj_c_totnowshuffled(tsnow>x(i)&tsnow<x(i+1)))...
            ./wskj_c_totnowtot;
    end

    wskj_c_tot_ts_bsqa(:,irun) = wskj_c_tot_freq_tsbinsnow ./ (n(1:end-1)/sum(n(1:end-1)));
end
wskj_c_tot_ts_bsqa_2pt5p = nan(length(n)-1,1);
wskj_c_tot_ts_bsqa_97pt5p = nan(length(n)-1,1);
for i = 1:(length(n)-1)
    wskj_c_tot_ts_bsqa_2pt5p(i) = quantile(wskj_c_tot_ts_bsqa(i,:),0.025);
    wskj_c_tot_ts_bsqa_97pt5p(i) = quantile(wskj_c_tot_ts_bsqa(i,:),0.975);
end

f=figure; set(f,'color','white'); ax = gca;
yyaxis left;
bar(x,n,'histc');
yyaxis right;
plot(x(1:end-1)+0.25,wskj_c_tot_ts_qa,'color','r','linewidth',3);hold on;
plot(x(1:end-1)+0.25,ones(length(x)-1,1),'color','r');
plot(x(1:end-1)+0.25,wskj_c_tot_ts_bsqa_2pt5p,'r--');
plot(x(1:end-1)+0.25,wskj_c_tot_ts_bsqa_97pt5p,'r--');
ax.YColor = 'r';

%----------------------------------------------------
% GAMs
%----------------------------------------------------

