### BRISC library for nearest neighbor Gaussian process (NNGP) approximation

This directory contains C/C++ implementation of NNGP approximation and its R interface (Saha and Datta, 2018). Code is adopted from <https://github.com/ArkajyotiSaha/BRISC>.

To compile the C/C++ code, run the following command in your terminal:
```
R CMD SHLIB BRISC/BRISC.cpp BRISC/BRISC_Predict.cpp BRISC/util.cpp BRISC/init.c BRISC/lbfgs.c
```

Reference:

Saha, A. and Datta, A. (2018). BRISC: bootstrap for rapid inference on spatial covariances. _Stat_, 7(1):e184.