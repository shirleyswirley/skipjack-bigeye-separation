#https://stackoverflow.com/questions/16592222/matplotlib-group-boxplots

bpl = axnow.boxplot(ennow, positions=np.array(range(len(ennow)))*2.0-0.4, sym='', widths=0.3)
bpm = axnow.boxplot(allnow, positions=np.array(range(len(allnow)))*2.0, sym='', widths=0.3)
bpr = axnow.boxplot(lnnow, positions=np.array(range(len(lnnow)))*2.0+0.4, sym='', widths=0.3)

set_box_color(bpl, 'red')
set_box_color(bpm, 'black')
set_box_color(bpr, 'blue')

# draw temporary red and blue lines and use them to create a legend
axnow.plot([], c='red', label='El Nino')
axnow.plot([], c='black', label='All')
axnow.plot([], c='blue', label='La Nina')
axnow.legend()

axnow.set_title(titlenow)
axnow.set_ylabel(unitsnow)
axnow.set_xlim(-2, len(ticks)*2)
axnow.set_xticks(range(0, len(ticks) * 2, 2))
axnow.set_xticklabels(ticks, rotation='vertical')
plt.tight_layout()