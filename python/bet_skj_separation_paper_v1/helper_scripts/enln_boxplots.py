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

# - Plot boxplot
bpl = axnow.boxplot(ennow, positions=leftpos, sym='', whis=0,
                    widths=widthnow, notch=notchnow,
                    patch_artist=patchnow,
                    boxprops={'color':'red','facecolor':'red'},
                    medianprops={'linestyle':'-','color':'white'})
bpm = axnow.boxplot(allnow, positions=midpos, sym='', whis=0,
                    widths=widthnow, notch=notchnow,
                    patch_artist=patchnow,
                    boxprops={'color':'black','facecolor':'black'},
                    medianprops={'linestyle':'-','color':'white'})
bpr = axnow.boxplot(lnnow, positions=rightpos, sym='', whis=0,
                    widths=widthnow, notch=notchnow,
                    patch_artist=patchnow,
                    boxprops={'color':'blue','facecolor':'blue'},
                    medianprops={'linestyle':'-','color':'white'})

# - Star significant EN/LN comparisons
if plotsignif==1:
    signifpos = signifidxsnow*xscfactor
    ymin, ymax = axnow.get_ylim()
    if signifmarkup==0:
        yplotnow = ymin + (ymax-ymin)/20
    elif signifmarkup==1:
        yplotnow = ymax - (ymax-ymin)/20
    axnow.scatter(signifpos, np.full_like(signifpos,yplotnow),
        s=30, marker='*', color=signifmarkcolnow)
    #medsnow = np.array(list(map(np.median, allnow)))
    #axnow.scatter(signifpos, medsnow[signifidxsnow],
    #              c='black', s=35, marker='*')
    
# - Draw temporary red and blue lines and use them to create a legend
if plotlegend==1:
    axnow.plot([], c='red', label='El Niño')
    axnow.plot([], c='black', label='All')
    axnow.plot([], c='blue', label='La Niña')
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