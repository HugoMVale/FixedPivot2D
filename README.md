# Description

This code is a Fortran 90 implementation of the *extended fixed pivot method* developed by [Vale & McKenna](https://doi.org/10.1021/ie050179s) to solve [population balance equations (PBE)](https://en.wikipedia.org/wiki/Population_balance_equation) for two-component aggregation processes (also known as bivariate aggregation).

The method is described in detail in [`Vale2005b.pdf`](/Vale2005b.pdf).

# How to use the code

1. Define the desired aggregation/coagulation kernel
- Function `COAGKERNEL` in file [`INPUTS.f90`](/code/INPUTS.f90)

2. Define initial particle size distribution (PSD)
- Function `PSDINI` in file [`INPUTS.f90`](/code/INPUTS.f90)

3. Define end integration time
- Variable: `time_final` in program `FIXEDPIVOT2D`in file [`MAIN.f90`](/code/MAIN.f90)

4. Define grid range
- Variables `xmin`, `xmax`, `ymin`, `ymax` in program `FIXEDPIVOT2D`in file [`MAIN.f90`](/code/MAIN.f90)

5. Define grid type
- Variables `meshx`, `meshy` in program `FIXEDPIVOT2D`in file [`MAIN.f90`](/code/MAIN.f90)
- See subroutine `GRID` in file [`GRID.f90`](/code/GRID.f90) for more information.

6. Define number of grid cells
- Variables `Mx`, `My` in module `GLOBAL` in file [`GLOBAL.f90`](/code/GLOBAL.f90)

7. Update folder path
- Subroutine `OUTPUT` in file [`OUTPUT.f90`](/code/OUTPUT.f90)

8. Build and run
- Do not forget to include [DLSODE or ODEPACK](https://computing.llnl.gov/casc/odepack/), 
e.g. [`libODEPACK.a`](/code/libODEPACK.a) 
