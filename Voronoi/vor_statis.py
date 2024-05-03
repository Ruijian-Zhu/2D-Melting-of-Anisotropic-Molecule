# Perform Voronoi analysis for lammps trajectry via Freud library developed by Prof. Sharon Glotzer's group
# The freud package need to be installed a priori
# Detailed illustration can be found in the manual of freud library 
# This code is used to calculate the number of edge and volume (2D: area) of Voronoi cells for all snapshots 
# Written by Ruijian Zhu (ITP-CAS), Feb 2024

import warnings

import freud
import numpy as np
from matplotlib import pyplot as plt

frame = 2000    # number of snapshots
N = 2496
for i in range(frame):
    count = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # We read the number of particles, the system box, and the particle positions into 3 separate arrays.
        box_data = np.genfromtxt("dump_2.89_reorder.lammpstrj", skip_header=5+(N+9)*i, max_rows=3)
        data = np.genfromtxt("dump_2.89_reorder.lammpstrj", skip_header=9+(N+9)*i, max_rows=N)
    box_data[2] = [0,0]
    box = freud.box.Box.from_box(box_data[:,1]-box_data[:,0])
    cols = [3, 4, 5]
    new_arr = data[:,cols]
    new_arr[:,0] -= box_data[0][1]-box.Lx/2
    new_arr[:,1] -= box_data[1][1]-box.Ly/2
    voro = freud.locality.Voronoi()
    cells = voro.compute((box,new_arr)).polytopes
    filename = "vor" + str(i) + ".txt"  # edges
    # create a numpy array to save polygon edge
    polygon_edge = np.zeros(N)
    for j in range(N):
        polygon_edge[j] = len(cells[j])
    np.savetxt(filename, polygon_edge, delimiter='\t')
    filena = "vol" + str(i) + ".txt"    # areas
    np.savetxt(filena, voro.volumes, delimiter='\t')
    del box_data
    del box
    del data
    del new_arr
    del cells
    del filename
    del filena
    print(i)