#https://matplotlib.org/3.1.0/gallery/statistics/boxplot.html
# - Define boxplot x-axis positions
xscfactor = 2
xoffsetfactor = 0.4
midpos = np.array(range(len(allnow)))*xscfactor 
leftpos = midpos-xoffsetfactor
rightpos = midpos+xoffsetfactor 

# - Define boxplot visual params
widthnow = 0.3
notchnow = False
patchnow = True
meanlinewidth = 1
meanlinecol = 'gold'
meanlinestyle = '-'
signifmarkcolnow = 'black'

# - Plot boxplot
bpl = axnow.boxplot(ennow, positions=leftpos, sym='', whis=0,
                    widths=widthnow, notch=notchnow,
                    patch_artist=patchnow,
                    showmeans=True, meanline=True, 
                    meanprops=dict(linewidth=meanlinewidth, color=meanlinecol,
                                   linestyle=meanlinestyle),
                    boxprops={'color':'red','facecolor':'red'},
                    medianprops={'linestyle':'-','color':'white'})
bpm = axnow.boxplot(allnow, positions=midpos, sym='', whis=0,
                    widths=widthnow, notch=notchnow,
                    patch_artist=patchnow,
                    showmeans=True, meanline=True, 
                    meanprops=dict(linewidth=meanlinewidth, color=meanlinecol,
                                   linestyle=meanlinestyle),
                    boxprops={'color':'black','facecolor':'black'},
                    medianprops={'linestyle':'-','color':'white'})
bpr = axnow.boxplot(lnnow, positions=rightpos, sym='', whis=0,
                    widths=widthnow, notch=notchnow,
                    patch_artist=patchnow,
                    showmeans=True, meanline=True, 
                    meanprops=dict(linewidth=meanlinewidth, color=meanlinecol,
                                   linestyle=meanlinestyle),
                    boxprops={'color':'blue','facecolor':'blue'},
                    medianprops={'linestyle':'-','color':'white'})

# - Print out medians and means
print(titlenow + ' left medians, means:')
print([item.get_ydata()[1] for item in bpl['medians']])
print([item.get_ydata()[1] for item in bpl['means']])
print('middle medians, means:')
print([item.get_ydata()[1] for item in bpm['medians']])
print([item.get_ydata()[1] for item in bpm['means']])
print('right medians, means:')
print([item.get_ydata()[1] for item in bpr['medians']])
print([item.get_ydata()[1] for item in bpr['means']])

# - Change y-axis lims
if setylimnow==1:
    axnow.set_ylim([yminnow,ymaxnow])
    
# - Star significant EN/LN comparisons
if plotsignif==1:
    signifpos = signifidxsnow*xscfactor
    yplotnow = 0
    axnow.scatter(signifpos, np.full_like(signifpos,yplotnow),
        s=35, marker='*', color=signifmarkcolnow)
    
# - Draw temporary red and blue lines and use them to create a legend
if plotlegend==1:
    axnow.plot([], c='red', label='El Niño')
    axnow.plot([], c='black', label='All')
    axnow.plot([], c='blue', label='La Niña')
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