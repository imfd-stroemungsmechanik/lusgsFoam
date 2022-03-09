# lusgsFoam

## Features

- Density-based solver for high-speed, transient, inviscid/viscous flows
- Support for dynamic mesh movement, adaptive mesh refinement and Multiple Reference Frame (MRF) for turbomachinery flows
- Implicit time advancement Based on the implicit *Lower-Upper Symmetric Gauss-Seidel* (LU-SGS) time integration algorithm
- 2nd order spatial accuracy using the using the central schemes of Kurganov and Tadmor
- Convective flux Jacobian computed based on matrix formulation from Blazek (2015)
- Additional performance and loss coefficient function object for turbomachinery postprocessing
- Support for OpenFOAM v2106 and v2112

## Compilation

1. Clone the ```lusgsFoam``` repository

```bash
git clone https://github.com/imfd-stroemungsmechanik/lusgsFoam.git
```

2. Use the Allwmake script to compile the `lusgsFoam` solver and the function objects

```bash
./Allwmake
```

## Validation Cases

Two turbomachienry cases are available for validation including the complete numerical setup:
1. VTI turbine cascade LS89 in ```run/ls89```
2. high-speed centrifugal compressor RadiVer in ```run/radiver``` 
