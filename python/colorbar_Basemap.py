cbar = plt.colorbar(cs, orientation='vertical', shrink=1)
cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=cbarlabrot)
cbar.ax.tick_params(labelsize=cbarlabsize)
