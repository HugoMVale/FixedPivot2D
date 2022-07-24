FUNCTION COAGKERNEL(x1,y1,x2,y2,time)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This is where the coagulation/aggregation kernel is defined. To adapt as desired.
!See for instance Equation 14.
!
!VARIABLES:
!
!B0     aggregation rate factor
!eij    sticking probability between pure particles of i and j
!eff    efficiency
!fxi    mass fraction of component i
!vi     mass of particle i
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

DOUBLE PRECISION, INTENT (IN) :: x1, y1, x2, y2, time
DOUBLE PRECISION :: COAGKERNEL


DOUBLE PRECISION, PARAMETER :: B0=1.0d0
DOUBLE PRECISION :: v1, v2, eff, fx1, fx2, e11=0.8d0, e22=0.1, e12=0.6d0

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

v1 = x1 + y1
v2 = x2 + y2

fx1 = x1/v1
fx2 = x2/v2

eff = 1.0d0 !e11*fx1*fx2 + ((1.0d0-fx1)*fx2 + fx1*(1.0d0-fx2))*e12 + (1.0d0-fx1)*(1.0d0-fx2)*e22


COAGKERNEL = B0 !*eff*(2.0d0 + (v1/v2)**(1.0d0/3.0d0) + (v2/v1)**(1.0d0/3.0d0))

END FUNCTION COAGKERNEL
!##############################################################################################








FUNCTION PSDINI(xf,yf)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!Initial condition for the number density function.
!See, for example, Equations 10 and 11 of paper.
!
!VARIABLES:
!
!mi     mass of component i
!mi0    mean mass of component i
!N0     initial number of particles per unit volume
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

DOUBLE PRECISION, INTENT (IN) :: xf, yf
DOUBLE PRECISION :: PSDINI

DOUBLE PRECISION :: m1, m2, m10, m20, N0

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

m1 = xf
m2 = yf

m10 = 1.0d0
m20 = 1.0d0
N0  = 1.0d0

PSDINI = 4.0d0*N0/(m10*m20)*(m1/m10)*DEXP(-2.0d0*(m1/m10) - (m2/m20))

END FUNCTION PSDINI
!##############################################################################################








FUNCTION EXPO2DMM(m1,m2,m10,m20,N0)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!2D Exponential distribution. Input: mass of each component.
!
!VARIABLES:
!
!mi     mass of component i
!mi0    mean mass of component i
!N0     initial number of particles per unit volume
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

DOUBLE PRECISION, INTENT (IN) :: m1, m2, m10, m20, N0
DOUBLE PRECISION :: EXPO2DMM

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


EXPO2DMM = N0/(m10*m20)*DEXP(-m1/m10 -m2/m20)


END FUNCTION EXPO2DMM
!##############################################################################################
