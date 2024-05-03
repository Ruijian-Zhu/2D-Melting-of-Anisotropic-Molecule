# 2D-Melting-of-Anisotropic-Molecule

These Cpp programs can be used for analyzing the two-dimensional phase transition of general systems. There are basically four physical quantities are evaluated: translational correlation function, bond-orientational correlation function, body-orientational order parameter (also the distribution of monomer body orientation), and the Voronoi analysis of topological defects. For the first two quantities, we provide two different programs for each, one is for rectangular box while the other is for triclinic box. All the input data should be created from LAMMPS. Below we will describe each program in more detail.

The codes are developed by Ruijian Zhu, from [Yanting Wang's research group] (http://soft-matter.itp.ac.cn/) at Institute of Theoretical Physics, Chinese Academy of Sciences (ITP-CAS). Please contact the developers Ruijian Zhu (zhuruijian@itp.ac.cn) for questions.

### Reorder the COM trajectory

The analysis tool provided by LAMMPS can calculate the COM of each monomer in each snapshot and output as an independent file, whose example is given in the corresponding file. To perform more analysis based on the time evolution of COMs (for most cases, these are also the representative points of unit-cells), we reorganize them into a single LAMMPS trajectory file, which is realized by the corresponding cpp program. The cpp file with a suffix 'tilt' is modified for triclinic box.

### Calculate the translational correlation function

The calculation of translational correlation is separated into two steps: (1) Scannng $\theta$ at a fixed $r$ to find out the crystal axis; (2) Scanning over $r$ along the crystal axis to obtain the correlation function. The first step is performed by theta.cpp while the second step is performed by paircorre.cpp, the file with a suffix 'tilt' is modified for triclinic box. A more detailed description on this procedure can be found in Y.-W. Li and M. P. Ciamarra, Phys. Rev. E 100, 062606 or [arXiv 2302.08305](https://arxiv.org/abs/2302.08305v3).

### Calculate the bond-orientational correlational function

The analysis tool provided by LAMMPS can calculate the bond-orietational order parameter related to each monomer, using the LAMMPS trajectory created in the first step as input, and each snapshot is output as an independent file. The cpp program here use these output files as well as the reorder trajectory as input, calculating the correlation function at different $r$.

### Calculate the body-orientation

This program use both the original trajectory and the COM trajectory as inputs, so as to save the claculation time. 'body_theta.cpp' calculate the body-orientation while 'body_orient.cpp' calculate the corresponding order parameter, both take n-fold rotational symmetry of monomer into account. This is designed for the ball-stick molecule we used in our research work [arXiv 2302.08305](https://arxiv.org/abs/2302.08305v3), but can be easily modified for other anisotropic molecules, e.g., for polygons, the coordinate of each atom can be replaced by the one for the vertex.

### Examples

In the folder 'example', we upload sample documents for both rectangular and triclinic box. There is a pdf file in each folder illustrating the relationship of relative files. It should be noticed that these files are just used as illustration of the input and output of the codes, produced from very short simulation without sufficient relaxation, so the results may be not strict for physics.

