# Surface-Data-Augmentation

Laplace-Beltrami eigenfunction Data Augmentation (LB-eigDA) and Chebyshev
polynomial Data Augmentation (C-pDA) 

(C) 2020  Shih-Gu Huang    shihgu@gmail.com
          Anqi Qiu         bieqa@nus.edu.sg
          National University of Singapore

This package contains codes implementing the Laplace-Beltrami
eigenfunction Data Augmentation (LB-eigDA) and Chebyshev polynomial Data
Augmentation (C-pDA) to generate new data on surfaces. The mathematical
detail is given in [1].  


--------------------------------------------------------------------------
REFERENCES:
[1] Huang, S.-G., Chung, M.K., Qiu, A.:
Fast Mesh Data Augmentation via Chebyshev Polynomial of Spectral filtering.
arXiv:2010.02811, 2020.

[2] Hammond, D.K., Vandergheynst, P., Gribonval, R.: Wavelets on graphs via 
spectral graph theory. Applied and Computational Harmonic Analysis 30, 129-150, 2011

[3] Tan, M., Qiu, A.: Spectral Laplace-Beltrami wavelets with applications
in medical images. IEEE Transactions on Medical Imaging 34, 1005-1017, 2015


--------------------------------------------------------------------------
DEMONSTRATIONS: 
- EXAMPLE.m     : This example performs LB-eigDA and C-pDA to generate
                  augmented data from simulated data on a hippocampus
                  surface mesh [1] and generate Figure 2 in [1].


--------------------------------------------------------------------------
FUNCTIONS:
- eigDA.m       : Laplace-Beltrami eigenfunction Data Augmentation (LB-eigDA)

- CpDA.m        : Chebyshev polynomial Data Augmentation (C-pDA). The
                  recurrence relation is modified from the Spectral Graph
                  Wavelet Transform (SGWT) [2].  

- LB/           : Folder containing some functions modified from Spectral
                  Laplace-Beltrami Wavelets [3] for computing LB-operator 

- LB_operator.m : Discretization of Laplace-Beltrami operator

- permutation.m : Generate random permutations for LB-eigDA and C-pDA

- eigen.m       : Compute all or the first K eigenfunctions for LB-eigDA 

- BP_chebyshev.m: Compute Chebyshev expansion coefficients of all
                  bandpass filters in C-pDA 

- figure_trimesh.m : Visualization of surface data

- hippocampus_l.mat: Left hippocampus surface mesh


