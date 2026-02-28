import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, ticker
import math
import sys

np.set_printoptions(precision=3)

## Style
rc('text', usetex=True)
font = {'family' : 'normal', 'size' : 20}
rc('font', **font)

# n site cluster rank r
n=int(sys.argv[2])
r=int(sys.argv[3])



# threshold for feature highlighting; set to negative to turn off
noise_threshold = 0.4


comp = ['x','y', 'z']
#comp = ['x','y', 'z', 'x2', 'y2', 'z2']
numComps = len(comp)

def phiDim(nspins, rank):
    binomial = math.factorial(nspins) // math.factorial(rank) // math.factorial(nspins - rank)
    return binomial * numComps**rank

def ind_invalid(ind):
    for i in range(ind.size-1):
        for j in range(i+1,ind.size):
            if ind[i] >= ind[j] or ind[i]//numComps == ind[j]//numComps:
                return True
    return False

def advance_ind(ind):
    for ni in range(ind.size):
        nirev = ind.size-1-ni
        idx = ind[nirev]
        if idx+1 == numComps*n:
            l_ind = ind[:nirev]
            maxval = np.amax(l_ind[l_ind < numComps*n-1])//numComps * numComps
            ind[nirev] = maxval
            if maxval//numComps < n-1:
                ind[nirev] += numComps
        else:
            ind[nirev] = idx+1
            return ind
    return ind


def create_feature(ind):
    feat = ['1']*n
    for idx in ind:
        feat[idx//numComps] = comp[idx%numComps]
    return " ".join([c for c in feat])



data = np.genfromtxt(sys.argv[1])
assert data.size == phiDim(n,r)
print(data.size)


## this if-else-block has a problem handling rank 1
if r > 1:
    if n == r:
        feature_indices = np.where(np.abs(data) > noise_threshold)[0].flatten()
        for idx in feature_indices:
            coeff = data[idx]
            feature = " ".join([comp[(idx//(numComps**k)) % numComps] for k in range(r)])
            #print(f"[{idx}]\t\t\t{coeff:.2f}\t\t\t" + feature)
            print(f"[{idx}]\t\t{feature}\t\t{coeff:.2f}")
    else:
        ## NOTE: this is excessive and should be implemented in a more efficient way
        print("Compiling feature list...")
        feature_list = []
        ind = np.arange(0,numComps*r,numComps)
        while len(feature_list) != phiDim(n,r):#!= data.size:
            while ind_invalid(ind):
                ind = advance_ind(ind)
            feature_list.append(create_feature(ind))
            ind = advance_ind(ind)

        print(f"feature list length {len(feature_list)}")

        relevant_feat_inds = np.where(np.abs(data) > noise_threshold)[0].flatten()
        for idx in relevant_feat_inds:
            coeff = data[idx]
            feature = feature_list[idx]
            print(f"[{idx}]\t\t{feature}\t\t{coeff:.2f}")




fig = plt.figure()
ax = fig.add_subplot(111)

ax.tick_params(axis='both', which='both', direction='in',
    bottom=True, top=False, labelbottom=True, labeltop=False,
    left=True, right=False, labelleft=True, labelright=False)

#ax.xaxis.set_major_locator(ticker.MultipleLocator(data.size//6))
#ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))
#ax.set_ylim(ylims)

# Separate noise from true signal
high_idc = np.argwhere(np.abs(data) > noise_threshold).T[0]
low_idc = np.argwhere(np.abs(data) <= noise_threshold).T[0]

ax.errorbar(high_idc, data[high_idc], c='royalblue', fmt='o', ms=6, mec='black', lw=2, alpha=1)
ax.errorbar(low_idc, data[low_idc], c='royalblue', fmt='o', ms=6, alpha=0.15)

ax.hlines(y=0, xmin=0, xmax=data.size, ls='--', lw=1, color='gray')

ax.set_title(sys.argv[1], fontsize=20)
ax.set_xlabel(r"$\mu$")
ax.set_ylabel("Coefficient")

plt.tight_layout()
plt.savefig(f"./features_{sys.argv[1]}.png", dpi=200, bbox_inches='tight')
#plt.show()
