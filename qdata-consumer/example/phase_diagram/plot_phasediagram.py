import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc, ticker
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable



edge_color = 'dodgerblue'

## Style
rc('text', usetex=True)
font = {'family' : 'normal', 'size' : 20}
rc('font', **font)


fig = plt.figure()

ax = fig.add_subplot(111)#, aspect=1)


gDict = {0:-1, 1:-0.9, 2:-0.8, 3:-0.7, 4:-0.6, 5:-0.5, 6:-0.4, 7:-0.3, 8:-0.2, 9:-0.1, 10:-0.01,
            11:0.1, 12:0.2, 13:0.3, 14:0.4, 15:0.5, 16:0.6, 17:0.7, 18:0.8, 19:0.9, 20:0.99}



edges = np.genfromtxt(f"./edges.txt")
fiedler = np.genfromtxt(f"./phases.txt")[len(gDict):2*len(gDict)]

fiedDict = dict(zip(fiedler[:,0], fiedler[:,1]))

for edge in edges:
    if abs(gDict[edge[0]] - gDict[edge[1]]) < 0.4:
        ax.plot([gDict[edge[0]], gDict[edge[1]]] , [fiedDict[edge[0]], fiedDict[edge[1]]],
                alpha=edge[2], c=edge_color, ls='-', lw=3, marker='')


ax.errorbar(x=np.fromiter(gDict.values(), dtype=float), y=fiedler[:,1],
    c='green', fmt='o', ms=6, mec='black', capsize=4, lw=2)


ax.hlines(y=0, xmin=-1, xmax=+1, ls='--', color='gray', lw=1.5)

ax.tick_params(axis='both', which='both', direction='in',
    bottom=False, top=False, labelbottom=True, labeltop=False,
    left=False, right=False, labelleft=True, labelright=False, labelsize=18)

#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.4))
#ax.yaxis.set_major_locator(ticker.MultipleLocator(0.6))

ax.set_xlabel(r"$g$", labelpad=0)
ax.set_ylabel(r"\textsf{Fiedler value}", labelpad=0)



plt.tight_layout()
plt.savefig("ClusterPhases.pdf")
#plt.show()


