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

# - Plot boxplot
bpwin = axnow.boxplot(winnow, positions=winpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'blue','facecolor':'blue'},
                      medianprops={'linestyle':'-','color':'white'})
bpspr = axnow.boxplot(sprnow, positions=sprpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'green','facecolor':'green'},
                      medianprops={'linestyle':'-','color':'white'})
bpsum = axnow.boxplot(sumnow, positions=sumpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'red','facecolor':'red'},
                      medianprops={'linestyle':'-','color':'white'})
bpaut = axnow.boxplot(autnow, positions=autpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      boxprops={'color':'chocolate','facecolor':'chocolate'},
                      medianprops={'linestyle':'-','color':'white'})

# - Change y-axis lims
if setylimnow==1:
    axnow.set_ylim([yminnow,ymaxnow])

# - Star significant season comparisons
if plotsignif==1:
    signifpos = signifidxsnow*xscfactor
    ymin, ymax = axnow.get_ylim()
    axnow.scatter(signifpos, np.full_like(signifpos,ymin),
                  s=35, marker='*', color='black')

# - Draw temporary colored lines and use them to create a legend
axnow.plot([], c='blue', label='Winter')
axnow.plot([], c='green', label='Spring')
axnow.plot([], c='red', label='Summer')
axnow.plot([], c='chocolate', label='Fall')
axnow.legend(frameon=False, loc='upper left', )

# - Turn on the minor TICKS, which are required for the minor GRID
axnow.minorticks_on()
axnow.grid(which='major', linestyle='--', linewidth='0.5', color='gray', axis='y')
#axnow.grid(which='minor', linestyle='--', linewidth='0.5', color='gray', axis='y')

# - Turn off the display of ticks you don't want
axnow.tick_params(axis='x', which='minor', top=False, bottom=False)   

# - Add title, labels, ticks
axnow.set_title(titlenow)
axnow.set_ylabel(unitsnow)
axnow.set_xlim(-2, len(ticks)*2)
axnow.set_xticks(range(0, len(ticks) * 2, 2))
axnow.set_xticklabels(ticks, rotation='vertical')
plt.tight_layout()