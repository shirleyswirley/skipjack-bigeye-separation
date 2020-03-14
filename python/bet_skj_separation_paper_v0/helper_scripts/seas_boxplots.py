#https://stackoverflow.com/questions/16592222/matplotlib-group-boxplots

xscfactor = 2
xoffsetfactor = 0.18

midpos = np.array(range(len(winnow)))*xscfactor
winpos = midpos - 3*xoffsetfactor
sprpos = midpos - xoffsetfactor
sumpos = midpos + xoffsetfactor
autpos = midpos + 3*xoffsetfactor

bpwin = axnow.boxplot(winnow, positions=winpos, sym='', widths=0.2)
bpspr = axnow.boxplot(sprnow, positions=sprpos, sym='', widths=0.2)
bpsum = axnow.boxplot(sumnow, positions=sumpos, sym='', widths=0.2)
bpaut = axnow.boxplot(autnow, positions=autpos, sym='', widths=0.2)

if plotsignif==1:
    signifpos = signifidxsnow*xscfactor
    ymin, ymax = axnow.get_ylim()
    axnow.scatter(signifpos, np.full_like(signifpos,ymin),
                  s=30, marker='*', color='black')

set_box_color(bpwin, 'blue')
set_box_color(bpspr, 'green')
set_box_color(bpsum, 'red')
set_box_color(bpaut, 'chocolate')

# Draw temporary colored lines and use them to create a legend
axnow.plot([], c='blue', label='Winter')
axnow.plot([], c='green', label='Spring')
axnow.plot([], c='red', label='Summer')
axnow.plot([], c='chocolate', label='Fall')
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