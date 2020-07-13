This repository contains the MATLAB files to reproduce the data and figures from "PoincarÃ© maps for multiscale physics discovery and nonlinear Floquet theory" by Jason J. Bramburger and J. Nathan Kutz (Physica D, 2020). Computations use the publicly available SINDy architecture found at https://faculty.washington.edu/kutz/page26/ and should be stored in a folder entitled 'Util'. This repository contains an alternative sparsifyDynamics.m program, entitled sparsifyDynamicsAlt.m, which uses the SVD to solve the overdetermined linear system associated to the model discovery.

The scripts associated to this repository are as follows:

- RC_Section.m: Mapping discovery for the RC circuit model in Section 3.1.

- Hopf_Section.m: Mapping discovery for the Hopf normal form model in Section 3.2.

- Logistic_Section.m: Mapping discovery for the singularly perturbed Logistic model in Section 3.3.

- Brusselator_Section.m: Mapping discovery for the driven Brusselator model in Section 3.4.

- Rossler_Section.m: Mapping discovery for the Rossler model in Section 3.5.

- Spiral_Section.m: Mapping discovery for the coarse-grained dynamics of the dominant PCA modes of spiral waves. Spiral waves are generated from a lambda-omega model while PCA time series data is provided in spiral_data.mat.
