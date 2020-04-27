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
meanlinewidth = 1
meanlinecol = 'gold'
meanlinestyle = '-'
signifmarkcolnow = 'black'

# - Plot boxplot
bpwin = axnow.boxplot(winnow, positions=winpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      showmeans=True, meanline=True,
                      meanprops=dict(linewidth=meanlinewidth, color=meanlinecol,
                                   linestyle=meanlinestyle),
                      boxprops={'color':'blue','facecolor':'blue'},
                      medianprops={'linestyle':'-','color':'white'})
bpspr = axnow.boxplot(sprnow, positions=sprpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      showmeans=True, meanline=True,
                      meanprops=dict(linewidth=meanlinewidth, color=meanlinecol,
                                   linestyle=meanlinestyle),
                      boxprops={'color':'green','facecolor':'green'},
                      medianprops={'linestyle':'-','color':'white'})
bpsum = axnow.boxplot(sumnow, positions=sumpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      showmeans=True, meanline=True,
                      meanprops=dict(linewidth=meanlinewidth, color=meanlinecol,
                                   linestyle=meanlinestyle),
                      boxprops={'color':'red','facecolor':'red'},
                      medianprops={'linestyle':'-','color':'white'})
bpaut = axnow.boxplot(autnow, positions=autpos, sym='', whis=0,
                      widths=widthnow, notch=notchnow,
                      patch_artist=patchnow,
                      showmeans=True, meanline=True,
                      meanprops=dict(linewidth=meanlinewidth, color=meanlinecol,
                                   linestyle=meanlinestyle),
                      boxprops={'color':'chocolate','facecolor':'chocolate'},
                      medianprops={'linestyle':'-','color':'white'})

# - Print out medians and means
print(titlenow + ' win medians, means:')
print([item.get_ydata()[1] for item in bpwin['medians']])
print([item.get_ydata()[1] for item in bpwin['means']])
print('spr medians, means:')
print([item.get_ydata()[1] for item in bpspr['medians']])
print([item.get_ydata()[1] for item in bpspr['means']])
print('sum medians, means:')
print([item.get_ydata()[1] for item in bpsum['medians']])
print([item.get_ydata()[1] for item in bpsum['means']])
print('aut medians, means:')
print([item.get_ydata()[1] for item in bpaut['medians']])
print([item.get_ydata()[1] for item in bpaut['means']])

# - Change y-axis lims
if setylimnow==1:
    axnow.set_ylim([yminnow,ymaxnow])

# - Star significant season comparisons
if plotsignif==1:
    signifpos = signifidxsnow*xscfactor
    yplotnow = 0
    axnow.scatter(signifpos, np.full_like(signifpos,yplotnow),
                  s=35, marker='*', color=signifmarkcolnow)

# - Draw temporary colored lines and use them to create a legend
if plotlegend==1:
    axnow.plot([], c='blue', label='Winter')
    axnow.plot([], c='green', label='Spring')
    axnow.plot([], c='red', label='Summer')
    axnow.plot([], c='chocolate', label='Fall')
    axnow.legend(frameon=False, loc='upper left', )

# - Turn on the minor TICKS, which are required for the minor GRID
axnow.minorticks_on()
axnow.grid(which='major', linestyle='--', linewidth='0.5', color='gray', axis='y')

# - Turn off the display of ticks you don't want
axnow.tick_params(axis='x', which='minor', top=False, bottom=False)   

# - Add title, labels, ticks
axnow.set_title(titlenow)
axnow.set_ylabel(unitsnow)
axnow.set_xlim(-2, len(ticks)*2)
axnow.set_xticks(range(0, len(ticks) * 2, 2))
axnow.set_xticklabels(ticks, rotation='vertical')
plt.tight_layout()