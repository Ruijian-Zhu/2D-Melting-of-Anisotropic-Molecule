# Perform Voronoi analysis for lammps trajectry via Freud library developed by Prof. Sharon Glotzer's group
# The freud package need to be installed a priori
# Detailed illustration can be found in the manual of freud library 
# This code is used to plot the Voronoi diagram for one snapshot 
# Written by Ruijian Zhu (ITP-CAS), Feb 2024

import warnings
import freud
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import matplotlib.cm as cm

N = 2496        #number of COMs
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    box_data = np.genfromtxt("dump_2.92_reorder.lammpstrj", skip_header=5+(N+9)*1500, max_rows=3)
    data = np.genfromtxt("dump_2.92_reorder.lammpstrj", skip_header=9+(N+9)*1500, max_rows=N)
box_data[2] = [0,0]
box = freud.box.Box.from_box(box_data[:,1]-box_data[:,0])
cols = [3, 4, 5]
new_arr = data[:,cols]
# The freud library forces the box to locate at the center
new_arr[:,0] -= box_data[0][1]-box.Lx/2
new_arr[:,1] -= box_data[1][1]-box.Ly/2
voro = freud.locality.Voronoi()
cells = voro.compute((box,new_arr)).polytopes

# Plot the Voronoi diagram
fig = plt.figure(figsize=(8,5))
ax = fig.gca()
voro.plot(ax=ax, cmap="RdBu")
ax.scatter(new_arr[:, 0], new_arr[:, 1], s=0.1, c="k")
ax.axis('off')

# Save the figure
plt.savefig("vor.jpg",dpi=1200,format="jpg",bbox_inches='tight')
plt.savefig("vor.svg",dpi=1200,format="svg",bbox_inches='tight')
plt.show()
