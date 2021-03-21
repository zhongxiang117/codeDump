#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt


# Charge filtration plot in codes illustration for our paper:
#
# Zhong, X.; Velez, C.; Acevedo, O.
#
# Partial Charges Optimized by Genetic Algorithms for Deep Eutectic 
# Solvent Simulations. (Under Review)


RAWDATA = [
    0.59, 0.92, 0.06, 0.37, 0.85, 0.47, 0.50, 0.33, 0.65, 0.82,
    0.88, 0.37, 0.64, 0.73, 0.11, 0.85, 0.71, 0.61, 0.08, 0.33,
    0.14, 0.50, 0.68, 0.57, 0.92, 0.48, 0.63, 0.73, 0.79, 0.75,
    0.84, 0.99, 0.55, 0.35, 0.22, 0.75, 0.98, 0.05, 0.89, 0.35,
    0.63, 0.88, 0.81, 0.14, 0.88, 0.35, 0.58, 0.94, 0.34, 0.16,
    0.07, 0.36, 0.08, 0.52, 0.08, 0.01, 0.91, 0.07, 0.75, 0.29,
    0.98, 0.51, 0.31, 0.55, 0.98, 0.72, 0.76, 0.86, 0.04, 0.22,
    0.09, 0.66, 0.91, 0.94, 0.08, 0.78, 0.33, 0.90, 0.94, 0.02,
    0.71, 0.02, 0.72, 1.00, 0.79, 0.42, 0.89, 0.80, 0.02, 0.65,
    0.56, 0.47, 0.76, 0.57, 0.81, 0.19, 0.31, 0.36, 0.02, 0.21,
]


stepsize = 0.04
rmin = min(RAWDATA)
nm = int((max(RAWDATA)-rmin)/stepsize + 0.5)
# note, since the round-of-error is not considered, so the total
# number of data in bins may be different with that in RAWDATA
data = np.array(RAWDATA)
charges = [rmin+stepsize*n for n in range(nm)]
bins = []
for n in range(nm):
    t = np.logical_and(data>=charges[n], data<charges[n]+stepsize)
    bins.append(np.sum(t))


bmax = max(bins)
mode = rmin + bins.index(bmax)*stepsize


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter([stepsize*n for n in range(len(data))], data, zdir='y')
# to better view, bins is normalized
norbins = [i/bmax for i in bins]
ax.plot(charges, norbins, zdir='x')


ax.set_zlim([0,1])
ax.set_xlabel('stepsize')
ax.set_ylabel('bins')
ax.set_zlabel('charges')


bo = False
try:
    from GAML import function_gen_range as fn
except ImportError:
    bo = True

if bo:
    info = 'Note: gaml is not installed\n'
    info += 'for more info, please install it by\n'
    info += '=> pip3 install gaml'
    ax.text(-0.8, -1.0, 1.8, info, color='black')
else:
    bins = [int(i) for i in bins]
    ndxmin,ndxmax = fn.func_gen_range(bins,0.8)
    chmin = rmin + ndxmin*stepsize
    chmax = rmin + ndxmax*stepsize

    px, py = np.meshgrid(np.linspace(0.0,4.5,10), np.linspace(0.0,1.0,10))
    pzmin = px*0 + chmin
    ax.plot_surface(px,py,pzmin,alpha=0.5,color='cyan')
    ax.text(2.0, 0.5, chmin, 'lower bound')

    pzmax = px*0 + chmax
    ax.plot_surface(px,py,pzmax,alpha=0.5,color='green')
    ax.text(2.0, 0.5, chmax, 'higher bound')

    pzmode = px*0 + mode
    ax.plot_surface(px,py,pzmode,alpha=0.5,color='purple')
    ax.text(2.0, 0.5, mode, 'mode')

    sx, sz = np.meshgrid(np.linspace(0.0,4.5,10), np.linspace(chmin,chmax,10))
    sy = sx*0
    ax.plot_surface(sx,sy,sz,alpha=0.5,color='orange')
    ax.text(2.0, 0.0, (chmin+chmax)/2, '80%')


label = 'Folding stepsize & bins planes'
ax.text(4.0, 0.0, 0.8, label, color='red')


# set view angle, make it easy to see plot
ax.view_init(azim=-60, elev=30)


plt.title('Charge Filtration Plot In 3D Illustration')
plt.show()


