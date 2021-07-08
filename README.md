# Fingerprint
The module fingerprint can be used for the computation of overlap matrix based fingerprints. The fingerprint vector contains the eigenvalues of the overlap matrix which is itself constructed from gaussian type orbitals. More detailed information can be found in [[1]](#1) and [[2]](#2). The fingerprint can be calculated for periodic structures (src/fp_for_crystals.f90) or for clusters (src/fp_for_clusters.f90). For the calculation of non-periodic systems the lattice vectors (alat) have to be set to zero. For periodic systems the matrix alat contains columnwise the lattice vectors.

## Periodic Version
The periodic version (src/fp_for_crystals.f90) is a template program to calculate the fingerprint of an Alanate crystal. If only the fingerprint without the derivative shall be calculated please set only_fp = .true., otherwise also the derivative of the fingerprint will be calculated and written into a file. Compiling this version can be done in the following way:
```bash
mkdir build
cd build
cmake -DDEBUG=OFF -DCRYSTAL=ON .. # compile without debug flags and the gfortran compiler
make
```
If you want to use the intel fortran compiler please use:
```bash
mkdir build
cd build
cmake -DDEBUG=OFF -DINTEL=ON -DCRYSTAL=ON .. # compile without debug flags and the intel compiler
make
```

## Clusters
The cluster version (src/fp_for_clusters.f90) is a template program to calculate the fingerprint of a benzene molecule. The main difference in this version is that the unit cell vectors are set to zero. If only the fingerprint without the derivative shall be calculated please set only_fp = .true., otherwise also the derivative of the fingerprint will be calculated and written to files. Compiling this version can be done in the following way:
```bash
mkdir build
cd build
cmake -DDEBUG=OFF -DCLUSTER=ON .. # compile without debug flags and the gfortran compiler
make
```
If you want to use the intel fortran compiler please use:
```bash
mkdir build
cd build
cmake -DDEBUG=OFF -DINTEL=ON -DCLUSTER=ON .. # compile without debug flags and the intel compiler
make
```
## Fingerprint Distance
In this template program (src/fp_distance.f90) the fingerprint distance of two periodic perovskite crystals is calculated using the hungarian algorithm. The hungarian algorithm is taken from another [GitHub repository](https://github.com/Jonas-Finkler/RMSD-finder) and computes the translation, rotation and permutation of atoms that minimizes the RMSD between two crystalline configurations. Compiling this version can be done in the following way:
```bash
mkdir build
cd build
cmake -DDEBUG=OFF -DDISTANCE=ON .. # compile without debug flags and the gfortran compiler
make
```
If you want to use the intel fortran compiler please use:
```bash
mkdir build
cd build
cmake -DDEBUG=OFF -DINTEL=ON -DDISTANCE=ON .. # compile without debug flags and the intel compiler
make
```

## References

Please cite the following articles if you use this code in an academic setting:

<a id="1">[1]</a> 
Zhu, L. et al. (2016).
A fingerprint based metric for measuring similarities of crystalline structures.
J. Chem. Phys, 144(3)

<a id="2">[2]</a> 
Sadeghi, A. et al. (2013).
Metrics for measuring distances in configuration spaces.
J. Chem. Phys, 139(18)
