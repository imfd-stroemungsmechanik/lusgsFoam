# lusgsFoam

## Features

- Density-based solver for high-speed, transient, inviscid/viscous flows
- Implicit time advancement allows for large Co-numbers
- Based on the implicit *Lower-Upper Symmetric Gauss-Seidel* (LU-SGS) scheme
- 2nd order spatial accuracy using the using the central schemes of Kurganov and Tadmor
- Accuracy of LUSGS scheme user-selectable by max number of iterations and relative tolerance
- Support for OpenFOAM v2006
