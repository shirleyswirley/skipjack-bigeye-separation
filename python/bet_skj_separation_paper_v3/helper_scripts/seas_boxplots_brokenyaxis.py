#https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/broken_axis.html
#https://stackoverflow.com/questions/5159065/need-to-add-space-between-subplots-for-x-axis-label-maybe-remove-labelling-of-a

# - Define boxplot x-axis position
xscfactor = 2
xoffsetfactor = 0.18
midpos = np.array(range(len(winnow)))*xscfactor
winpos = midpos - 3*xoffsetfactor
sprpos = midpos - xoffsetfactor
sumpos = midpos + xoffsetfactor
autpos = midpos + 3*xoffsetfactor

# - Define boxplot visual params
widthnow = 0.2
notchnow = False
patchnow = True

# - Plot top boxplot
bpwin = axtopnow.boxplot(winnow, positions=winpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'blue','facecolor':'blue'},
                      medianprops={'linestyle':'-','color':'white'})
bpspr = axtopnow.boxplot(sprnow, positions=sprpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'green','facecolor':'green'},
                      medianprops={'linestyle':'-','color':'white'})
bpsum = axtopnow.boxplot(sumnow, positions=sumpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'red','facecolor':'red'},
                      medianprops={'linestyle':'-','color':'white'})
bpaut = axtopnow.boxplot(autnow, positions=autpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'chocolate','facecolor':'chocolate'},
                      medianprops={'linestyle':'-','color':'white'})
axtopnow.set_ylim([ymintopnow,ymaxnow])

# - Plot bottom boxplot
bpwin = axbotnow.boxplot(winnow, positions=winpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'blue','facecolor':'blue'},
                      medianprops={'linestyle':'-','color':'white'})
bpspr = axbotnow.boxplot(sprnow, positions=sprpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'green','facecolor':'green'},
                      medianprops={'linestyle':'-','color':'white'})
bpsum = axbotnow.boxplot(sumnow, positions=sumpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'red','facecolor':'red'},
                      medianprops={'linestyle':'-','color':'white'})
bpaut = axbotnow.boxplot(autnow, positions=autpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'chocolate','facecolor':'chocolate'},
                      medianprops={'linestyle':'-','color':'white'})
axbotnow.set_ylim([yminnow,ymaxbotnow])

# - Hide the spines btwn axtopnow and axbotnow
axtopnow.spines['bottom'].set_visible(False)
axbotnow.spines['top'].set_visible(False)
axtopnow.xaxis.tick_top()
axtopnow.tick_params(labeltop=False) # don't put tick labels at the top
axbotnow.xaxis.tick_bottom()

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=axtopnow.transAxes, color='k', clip_on=False)
axtopnow.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
axtopnow.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=axbotnow.transAxes)  # switch to the bottom axes
axbotnow.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
axbotnow.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

# - Star significant season comparisons
if plotsignif==1:
    signifpos = signifidxsnow*xscfactor
    ymin, ymax = axbotnow.get_ylim()
    axbotnow.scatter(signifpos, np.full_like(signifpos,ymin),
                  s=35, marker='*', color='black')
    
# - Draw temporary colored lines and use them to create a legend
if plotlegend==1:
    axtopnow.plot([], c='blue', label='Winter')
    axtopnow.plot([], c='green', label='Spring')
    axtopnow.plot([], c='red', label='Summer')
    axtopnow.plot([], c='chocolate', label='Fall')
    axtopnow.legend(frameon=False, loc='upper left', )

# - Turn on the minor TICKS, which are required for the minor GRID
axtopnow.minorticks_on()
axtopnow.grid(which='major', linestyle='--', linewidth='0.5', color='gray', axis='y')
axbotnow.minorticks_on()
axbotnow.grid(which='major', linestyle='--', linewidth='0.5', color='gray', axis='y')

# - Turn off the display of ticks you don't want
axtopnow.tick_params(axis='x', which='minor', top=False, bottom=False)   
axbotnow.tick_params(axis='x', which='minor', top=False, bottom=False)   

# - Add title, labels, ticks
axtopnow.set_title(titlenow)
axtopnow.set_xlim(-2, len(ticks)*2)
axtopnow.set_xticks(range(0, len(ticks) * 2, 2))

axbotnow.set_ylabel(unitsnow)
axbotnow.set_xlim(-2, len(ticks)*2)
axbotnow.set_xticks(range(0, len(ticks) * 2, 2))
axbotnow.set_xticklabels(ticks, rotation='vertical')