#https://stackoverflow.com/questions/16592222/matplotlib-group-boxplots

xscfactor = 2
xoffsetfactor = 0.4

midpos = np.array(range(len(allnow)))*xscfactor 
leftpos = midpos-xoffsetfactor
rightpos = midpos+xoffsetfactor 

bpl = axnow.boxplot(ennow, positions=leftpos, sym='', widths=0.3)
bpm = axnow.boxplot(allnow, positions=midpos, sym='', widths=0.3)
bpr = axnow.boxplot(lnnow, positions=rightpos, sym='', widths=0.3)

if plotsignif==1:
    signifpos = signifidxsnow*xscfactor
    ymin, ymax = axnow.get_ylim()
    axnow.scatter(signifpos, np.full_like(signifpos,ymin),
                  s=30, marker='*', color='black')
    #medsnow = np.array(list(map(np.median, allnow)))
    #axnow.scatter(signifpos, medsnow[signifidxsnow],
    #              c='black', s=35, marker='*')
    
set_box_color(bpl, 'red')
set_box_color(bpm, 'black')
set_box_color(bpr, 'blue')

# Draw temporary red and blue lines and use them to create a legend
axnow.plot([], c='red', label='El Nino')
axnow.plot([], c='black', label='All')
axnow.plot([], c='blue', label='La Nina')
axnow.legend(frameon=False, loc='upper left', )

# Turn on the minor TICKS, which are required for the minor GRID
axnow.minorticks_on()
axnow.grid(which='major', linestyle='--', linewidth='0.5', color='gray', axis='y')
#axnow.grid(which='minor', linestyle='--', linewidth='0.5', color='gray', axis='y')

# Turn off the display of ticks you don't want
axnow.tick_params(axis='x', which='minor', top=False, bottom=False)   
        
axnow.set_title(titlenow)
axnow.set_ylabel(unitsnow)
axnow.set_xlim(-2, len(ticks)*2)
axnow.set_xticks(range(0, len(ticks) * 2, 2))
axnow.set_xticklabels(ticks, rotation='vertical')
plt.tight_layout()