for irow in range(0,len(dfatp)):
    if ~np.isnan(dfatp['ATP'].iloc[irow]):
        if dfatp['ATP'].iloc[irow]==-1:
            if cbfriendly:
                colornow = '#e41a1c'
            else:
                colornow = 'red'
        elif dfatp['ATP'].iloc[irow]==0:
            if cbfriendly:
                colornow = '#dede00'
            else: 
                colornow = 'yellow'
        elif dfatp['ATP'].iloc[irow]==1:
            if cbfriendly:
                colornow = '#4daf4a'
            else:
                colornow = 'green'
        ax.axvspan(dfatp['lbinedges'].iloc[irow], dfatp['rbinedges'].iloc[irow],
                    ymin=ysplits[ispcs], ymax=ysplits[ispcs+1], alpha=1, color=colornow)
        #start, end = ax.get_xlim()
        #ax.xaxis.set_ticks(np.linspace(start, end, 5))
        #ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%0.1f'))
        #ax.tick_params(axis='x', labelrotation=45)
        ax.set_xlabel(ivnicename + ' [' + ivunits + ']')
        ax.set_yticks(ysplits[0:nspcs] + (1/nspcs/2))
        ax.set_yticklabels(spcsnames)
        ax.axhline(1/nspcs, color='black')
        # Turn on the minor TICKS, which are required for the minor GRID
        ax.minorticks_on()
        ax.grid(which='major', linestyle='--', linewidth='0.5', color='gray', axis='x')
        ax.grid(which='minor', linestyle='--', linewidth='0.5', color='gray', axis='x')
        # Turn off the display of ticks you don't want
        ax.tick_params(axis='y', which='minor', left=False, right=False)   