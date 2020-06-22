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
                
        ax.axvspan(dfatp['lbinedges'].iloc[irow],
                   dfatp['rbinedges'].iloc[irow],
                   ymin=ysplits[ispcs], ymax=ysplits[ispcs+1],
                   alpha=1, color=colornow)
        
        ax.set_xlabel(ivnicename + ' [' + ivunits + ']')
        ax.set_yticks(ysplits[0:nspcs] + (1/nspcs/2))
        ax.set_yticklabels(spcsnamesatp)
        ax.axhline(1/nspcs, color='black')
        
        # - Turn on the minor TICKS, which are required for the minor GRID
        #fl = mticker.FixedLocator(np.linspace(ax.get_xlim()[0],
        #                                      ax.get_xlim()[1],
        #                                      len(binedges)))
        #ax.xaxis.set_minor_locator(fl)
        #ax.xaxis.grid(which='minor', color='gray',
        #              linestyle='--', linewidth=0.5)
        ax.minorticks_on()
        ax.grid(which='major', linestyle='--',
                linewidth='0.5', color='gray', axis='x')
        ax.grid(which='minor', linestyle='--',
                linewidth='0.5', color='gray', axis='x')
        
        # - Turn off the display of ticks you don't want
        ax.tick_params(axis='y', which='minor', left=False, right=False)   