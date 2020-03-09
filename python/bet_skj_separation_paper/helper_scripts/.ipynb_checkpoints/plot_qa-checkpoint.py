if plotbarorhist=='bar':
    binwidth = bincenters[1]-bincenters[0]
    ax.bar(bincenters-binwidth/6, ivcounts/ivcounts.sum(), binwidth/3, color='silver', label=ivnicename)
    ax.bar(bincenters+binwidth/6, dvcounts/dvcounts.sum(), binwidth/3, color='lightskyblue', label=dvnicename)
elif plotbarorhist=='hist':
    ax.hist(iv, bins=binedges, label='hist', color='white', edgecolor='black')
ax.set_ylabel(ivnicename + ' or ' + dvnicename + '\nFrequency')
ax.set_xlabel(ivnicename + ' [' + ivunits + ']')

quotcolor = 'mediumblue'

ax1 = ax.twinx()
ax1.plot(bincenters, dvquot, color=quotcolor, linestyle='-', label='quot')
ax1.plot(bincenters, dvquot, color=quotcolor, marker='o', markersize=5, label='_nolegend_')

if plotqlimsreplaceT==1: # Bernal et al 2007 PAPER
    ax1.plot(bincenters, qlimsreplaceT.iloc[0], color=quotcolor, linestyle='--', label='w/ repl')
    ax1.plot(bincenters, qlimsreplaceT.iloc[1], color=quotcolor, linestyle='--', label='_nolegend_')
if plotqlimsreplaceF==1: # Bernal et al 2007 CODE
    ax1.plot(bincenters, qlimsreplaceF.iloc[0], 'g--', label='w/o repl')
    ax1.plot(bincenters, qlimsreplaceF.iloc[1], 'g--', label='_nolegend_')

ax1.yaxis.label.set_color(quotcolor)
ax1.tick_params(axis='y', colors=quotcolor)
ax1.spines['right'].set_edgecolor(quotcolor)
ax1.set_ylabel(dvnicename + '/' + ivnicename + '\nQuotient')

# Turn on the minor TICKS, which are required for the minor GRID
ax.minorticks_on()
#ax.grid(which='major', linestyle='--', linewidth='0.5', color='gray', axis='x')
#ax.grid(which='minor', linestyle='--', linewidth='0.5', color='gray', axis='y')
# Turn off the display of ticks you don't want
#ax.tick_params(axis='x', which='minor', top=False, bottom=False)   

if plotlegend==1:
    # - Legend 1
    lines, labels = ax.get_legend_handles_labels()
    lines1, labels1 = ax1.get_legend_handles_labels()
    #ax1.legend(lines + lines1, labels + labels1, loc=2)
    ax1.legend(lines, labels, loc=2, framealpha=0, labelspacing=0, handletextpad=0.2)
    # - Legend 2
    #ax.legend([ivname + ' hist']); ax1.legend(['w/ repl','w/o repl'])
    # - Legend 3
    #fig.legend([ivname + ' hist','quotient','w/ repl','w/o repl'],
    #           bbox_to_anchor=(0.65, 0.1, 0.5, 0.5))