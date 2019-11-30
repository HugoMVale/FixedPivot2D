PROGRAM FIXEDPIVOT2D

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This program is a numerical implementation of the extended fixed pivot method to solve the PBE
!for 2-component aggregation, as described in https://doi.org/10.1021/ie050179s.
!
!DISCLAIMER:
!
!Some parts of the code were modified in 2019 to improve its readability and use. During this
!process some bugs may also have been introduced.
!
!KEY STEPS FOR USE:
!
!1. Define the desired coagulation kernel
!   Function COAGKERNEL in file INPUTS.f90
!
!2. Define initial PSD
!   Function PSDINI in file INPUTS.f90
!
!3. Define end integration time
!   Variable: time_final just below
!
!4. Define grid range
!   Variables: xmin, xmax, ymin, ymax just below
!
!5. Define grid type
!   Variables: meshx, meshy just below
!   See subroutine GRID in file GRID.f90 for more information.
!
!6. Define number of grid cells
!   Variables: Mx, My in module GLOBAL in file GLOBAL.f90
!
!7. Update folder path
!   Subroutine OUTPUT in file OUTPUT.f90
!
!8. Build and run
!   DO not forget to include DLSODE or ODEPACK
!   https://computing.llnl.gov/casc/odepack/
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER :: meshx, meshy

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!	*** END-TIME SPECIFICATION ***

time_final = 100.0d0

!	*** GRID SPECIFICATION ***

xmin = 1.0d-3
xmax = 1.0d4

ymin = 1.0d-3
ymax = 1.0d4

! For options, see GRID in GRID.f90
meshx=2
meshy=2


!	*** START OF CALCULATION SEQUENCE ***

!Start-of-simulation message
CALL OUTPUT (1)

!Build grid (different options available)
CALL GRID(meshx,meshy)

!Compute possible particle combinations and weights
CALL COMBINATIONS

!Set initial PSD
CALL AVERAGECELL

!Evaluate the coagulation matrix (for systems where coagulation kernel is time independent)
CALL COAGMATRIX

!Integrate discretized PBE
CALL INTEGRATION

!End-of-simulation message
CALL OUTPUT(4)


END PROGRAM FIXEDPIVOT2D
!##############################################################################################
