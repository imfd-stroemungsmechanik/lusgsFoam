# lusgsFoam

## Features

- Density-based solver for high-speed, transient, inviscid/viscous flows
- Support for dynamic mesh movement, adaptive mesh refinement and Multiple Reference Frame (MRF) for turbomachinery flows
- Implicit time advancement Based on the implicit *Lower-Upper Symmetric Gauss-Seidel* (LU-SGS) scheme
- 2nd order spatial accuracy using the using the central schemes of Kurganov and Tadmor
- 1st, 2nd and 3rd order MUSCL reconstruction schemes added (based on blastFoam)
- Additional performance function object to calculate (corrected) mass flow rate, total pressure ratio, and isentropic efficiency of a turbomachine

### Compilation

1. Clone the lusgsFoam repository

```bash
git clone git@gitserv.imfd.tu-freiberg.de:fluid/cfd/lusgsfoam.git
```

2. Use the Allwmake script to compile the `liblusgs` library, the `performance` function object, and the `lusgsFoam` solver

```bash
./Allwmake
```

