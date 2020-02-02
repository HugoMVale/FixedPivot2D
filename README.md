# Description

This code is a Fortran 90 implementation of the *extended fixed pivot method* developed by [Vale & McKenna](https://doi.org/10.1021/ie050179s) to solve [population balance equations (PBE)](https://en.wikipedia.org/wiki/Population_balance_equation) for two-component aggregation processes (also known as bivariate aggregation).

The method is described in detail in [`Vale2005b.pdf`](/Vale2005b.pdf).

# How to use the code

1. File [`INPUTS.f90`](/code/INPUTS.f90):
- Function `COAGKERNEL`
  - Define the desired aggregation/coagulation kernel.
- Function `PSDINI`
  - Define initial particle size distribution (PSD).

2. File [`MAIN.f90`](/code/MAIN.f90):
- Program `FIXEDPIVOT2D`
  - Variable `time_final`: define end integration time.
  - Variables `xmin`, `xmax`, `ymin`, `ymax`: define grid range.
  - Variables `meshx`, `meshy`: define grid type. See subroutine `GRID` in file [`GRID.f90`](/code/GRID.f90) for more information.

3. File [`GLOBAL.f90`](/code/GLOBAL.f90)
- Module `GLOBAL` 
  - Variables `Mx`, `My`: define number of grid cells.

4. File [`OUTPUT.f90`](/code/OUTPUT.f90)
- Subroutine `OUTPUT`
  - Update folder path.

5. Build and run
- Do not forget to include [DLSODE or ODEPACK](https://computing.llnl.gov/casc/odepack/), 
e.g. [`libODEPACK.a`](/code/libODEPACK.a) 
