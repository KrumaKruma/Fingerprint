# Fingerprint
The module fingerprint can be used for the computation of overlap matrix based fingerprints. The fingerprint vector contains the eigenvalues of the overlap matrix which is itself constructed from gaussian type orbitals. More detailed information can be found in [[1]](#1) and [[2]](#2). The fingerprint can be calculated for periodic structures (src/fp_for_crystals.f90) or for clusters (src/fp_for_clusters.f90). For the calculation of non-periodic systems the lattice vectors (alat) have to be set to zero. For periodic systems the matrix alat contains columnwise the lattice vectors.

## Periodic Version
The periodic version (src/fp_for_crystals.f90) is a template program to calculate the fingerprint of an Alanate crystal. If only the fingerprint without the derivative shall be calculated please set only_fp = .true., otherwise also the derivative of the fingerprint will be calculated and written into a file. Compiling this version can be done in the following way:
```bash
make clean
make
make fp_crystal
cd build
./fp_crystal.x
```
If you want to use the intel fortran compiler please use:
```bash
make clean
make FC=ifort
make fp_crystal FC=ifort
cd build
./fp_crystal.x
```

## Clusters
The cluster version (src/fp_for_clusters.f90) is a template program to calculate the fingerprint of a benzene molecule. The main difference in this version is that the unit cell vectors are set to zero. If only the fingerprint without the derivative shall be calculated please set only_fp = .true., otherwise also the derivative of the fingerprint will be calculated and written to files. Compiling this version can be done in the following way:
```bash
make clean
make
make fp_cluster
cd build
./fp_cluster.x
```
If you want to use the intel fortran compiler please use:
```bash
make clean
make FC=ifort
make fp_cluster FC=ifort
cd build
./fp_cluster.x
```

## Fingerprint Distance
In this template program (src/fp_distance.f90) the fingerprint distance of two periodic perovskite crystals is calculated using the hungarian algorithm. The hungarian algorithm is taken from another [GitHub repository](https://github.com/Jonas-Finkler/RMSD-finder) and computes the translation, rotation and permutation of atoms that minimizes the RMSD between two crystalline configurations. Compiling this version can be done in the following way:
```bash
make clean
make 
make fp_distance
cd build
./fp_distance.x
```
If you want to use the intel fortran compiler please use:
```bash
make clean
make FC=ifort
make fp_distance FC=ifort
cd build
./fp_distance.x
```

## Test the Fingerprint
There are two templates to test the fingerprint: 
Markup: 1. Template for profing rotational invariance by rotating the system randomly
        2. Path integration of the derivative of the fingerprint to test its correctness.
### Test rotational invariance
Compiling and running the test for rotationial invariance can be done the following way:
```bash
make clean
make 
make test_invariance
cd build
./test_invariance.x
```
If you want to use the intel fortran compiler please use:
```bash
make clean
make FC=ifort
make test_invariance FC=ifort
cd build
./test_invariance.x
```
The result for the fingerprint distance should be zero, whereas the result of the scalar product of the coordinates is non zero.
### Test fingerprint derivatives
Compiling and running the test for the fingerprint derivatives can be done the following way:
```bash
make clean
make 
make test_derivative
cd build
./test_derivative.x
```
If you want to use the intel fortran compiler please use:
```bash
make clean
make FC=ifort
make test_derivative FC=ifort
cd build
./test_derivative.x
```
The result is stored in the file Markup :  __Strong text__ or **Path_Integration.dat**. The result can be visualized with gnuplot:
```bash
gnuplot
plot "Path_Integration.dat" u 1:2 w l
replot "Path_Integration.dat" u 1:3 w l
```
If the derivative is correct the two curves should nearly be the same with very small deviations. The deviations should get smaller the more integration steps are used (nint) in the src/test_derivative.f90 file.


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
