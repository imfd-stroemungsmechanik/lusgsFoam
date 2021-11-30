# lusgsFoam

## Features

- Density-based solver for high-speed, transient, inviscid/viscous flows
- Support for dynamic mesh movement, adaptive mesh refinement and Multiple Reference Frame (MRF) for turbomachinery flows
- Implicit time advancement Based on the implicit *Lower-Upper Symmetric Gauss-Seidel* (LU-SGS) time integration algorithm
- 2nd order spatial accuracy using the using the central schemes of Kurganov and Tadmor
- Convective flux Jacobian computed based on matrix formulation from Blazek (2015)
- Additional performance function object to calculate (corrected) mass flow rate, total pressure ratio, and isentropic efficiency of turbomachines
- Support for OpenFOAM v2106

### Compilation

1. Clone the ```matrix-jacobian``` branch of the ```lusgsFoam``` repository

```bash
git clone -b matrix-jacobian git@gitserv.imfd.tu-freiberg.de:fluid/cfd/lusgsfoam.git
```

2. Use the Allwmake script to compile the `lusgsFoam` solver and the `performance` function object

```bash
./Allwmake
```

