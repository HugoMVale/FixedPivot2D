MODULE GLOBAL

IMPLICIT NONE
!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This module is used to define global variables that are shared among the various functions and
!subroutines. This is meant to avoid the use of common blocks.
!
!VARIABLES:
!
!area(i,j)          area of cell Cij
!Bmatrix(k,l,q,r)   matrix of coagulation rates
!dx(i)              width along x of cell Cij
!dy(j)              width along y of cell Cij
!F(i,j)             average number density in cell Cij
!mmij               moment of the number density function
!Mx                 number of cells along x-axis
!My                 number of cells along y-axis
!Np(i,j)            number of particles in cell Cij
!x(i)               x pivot of cell Cij
!y(j)               y pivot of cell Cij
!xright(i)          upper x-boundary of cell Cij
!yright(j)          upper y-boundary of cell Cij
!time_now           current time
!time_final         final solution time
!
!----------------------------------------------------------------------------------------------

DOUBLE PRECISION, PARAMETER :: Pi = 3.14159265d0

INTEGER, PARAMETER :: Mx=40, My=40

DOUBLE PRECISION :: F(Mx,My), Np(Mx,My)
DOUBLE PRECISION :: xmin, xmax, ymin, ymax
DOUBLE PRECISION :: x(0:Mx+1), y(0:My+1), xright(0:Mx), yright(0:My)
DOUBLE PRECISION :: dx(Mx), dy(My), area(Mx,My)
DOUBLE PRECISION :: time_now, time_final
DOUBLE PRECISION :: mm00, mm10, mm01, mm20, mm02, mm11, mm30, mm03, mm21, mm12
DOUBLE PRECISION :: Bmatrix(Mx,My,Mx,My)

INTEGER :: flag_output=0, last=0


TYPE COAG

    INTEGER (KIND=1), DIMENSION (:), POINTER :: k
    INTEGER (KIND=1), DIMENSION (:), POINTER :: q
    INTEGER (KIND=1), DIMENSION (:), POINTER :: l
    INTEGER (KIND=1), DIMENSION (:), POINTER :: r
    DOUBLE PRECISION, DIMENSION (:), POINTER :: eta

END TYPE COAG

TYPE (COAG), DIMENSION (Mx,My) :: cgmx

END MODULE GLOBAL
