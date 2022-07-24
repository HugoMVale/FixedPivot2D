FUNCTION MOMENT(nx,ny)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!Mixed moment of the number density function. See Equation B3 of paper.
!
!VARIABLES:
!
!nx     order of moment with respect to x
!ny     order of moment with respect to y
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT(IN) :: nx, ny
DOUBLE PRECISION :: MOMENT

DOUBLE PRECISION :: sum
INTEGER :: i, j

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

sum = 0.0d0

DO i=1,Mx

    DO j=1,My

        sum = sum + INTEGRAND_MOMENT(i,j,nx,ny)

    END DO

END DO

MOMENT = sum


!##############################################################################################


CONTAINS

    FUNCTION INTEGRAND_MOMENT(i,j,nx,ny)

    IMPLICIT NONE

    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    INTEGER, INTENT (IN) :: i, j, nx, ny
    DOUBLE PRECISION :: INTEGRAND_MOMENT

    !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


    INTEGRAND_MOMENT = Np(i,j)*(x(i)**nx)*(y(j)**ny)


    END FUNCTION INTEGRAND_MOMENT

END FUNCTION MOMENT
!##############################################################################################











