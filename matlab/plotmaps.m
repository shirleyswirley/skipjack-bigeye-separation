colormap(axnow,colormapnow); set(axnow,'fontweight','bold');
ax1=axesm('gortho','MapLatLimit',[slat nlat],'MapLonLimit',[wlon elon],...
  'Frame','off','Grid','off','MeridianLabel','off','ParallelLabel','off'); axis off;
m_proj(mapproj,'lon',[wlon elon],'lat',[slat nlat]);
m_pcolor(plonnow,platnow,mapvarnow);
m_grid('xtick',lonticks,'ytick',latticks,'xlabeldir','middle','fontsize',12);
if climsdefined==1; caxis(cminmax); end 
c=colorbar; set(c,'ticklength',0.05,'tickdirection','out'); shading flat;
c.Label.String = unitsnow;
if showeezs==1
    for ireg = 1:length(allpid)
        RR = double(eezmap0pt25==allpid{ireg});
        h1 = m_contour(loneez0pt25,lateez0pt25,RR,[.5 .5],'linecolor',eezltlinecolor,'linewidth',eezltlinewidth);
    end 
    nupidcell = [nupid{:,:}];
    for ireg = 1:length(nupidcell)
        RR = double(eezmap0pt25==nupidcell{ireg});
        h2 = m_contour(loneez0pt25,lateez0pt25,RR,[.5 .5],'linecolor',eezdklinecolor,'linewidth',eezdklinewidth);
    end 
end
m_coast('patch',landcolor,'edgecolor','k','linewidth',coastlinewidth); % must go after shading flat
tightmap;
title(titlenow);
