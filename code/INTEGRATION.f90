SUBROUTINE INTEGRATION

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!Here the particle PBEs are integrated with a popular ODE solver (DLSODE). The variables being
!solved for are the particle densities (not the particle numbers).
!
!VARIABLES:
!
!npoints    number of points where the solution will be saved (nothing to do with time step for
!           integration!)
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!VARIABLES CONCERNING DLSODE

INTEGER, PARAMETER :: neq = Mx*My
INTEGER, PARAMETER :: lrw = 22 + 16*neq + neq**2, liw = 20 + neq
INTEGER :: iopt, istate, itask, itol, mf, meth, miter, iwork(liw)
DOUBLE PRECISION :: atol(neq), rwork(lrw), rtol(neq), FNC(neq), time, time_output
EXTERNAL JAC, PBEBALANCES


INTEGER, PARAMETER :: npoints=20
DOUBLE PRECISION :: aux
INTEGER :: i, j, n

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!	*** DEFINITION OF THE CONDITIONS OF THE INITIAL VALUE PROBLEM ***

time_now = 0.0d0

n=1

DO i=1,Mx

	DO j=1,My

		FNC(n) = F(i,j)

		n=n+1

	END DO

END DO



!	*** RESOLUTION OF THE ODE SYSTEM ***


itol = 4
rtol = 1.0d-4
atol = 1.0d-18 !+ MAXVAL(F)*1.0d-5
itask = 1
istate = 1
iopt = 1
rwork = 0.0d0
iwork = 0
iwork(6) = 5000
meth = 1
miter = 0
mf = 10*meth+miter


DO i=0,npoints

	!IF (i==0) THEN

	!	time_output=0.0d0

	!ELSE

	!	time_output = 1.0d0*(tempo_final/1.0d0)**(DBLE(i-1)/DBLE(npoints-1))

	!END IF


	time_output = DBLE(i)/DBLE(npoints)*time_final


	CALL DLSODE(PBEBALANCES,neq,FNC,time_now,time_output,itol,rtol,atol,itask,istate, &
				iopt,rwork,lrw,iwork,liw,JAC,mf)

	flag_output = 1

	IF (i==npoints) last=1

	CALL PBEBALANCES(neq,time,FNC,aux)

	flag_output = 0

END DO


END SUBROUTINE INTEGRATION
!##############################################################################################









SUBROUTINE PBEBALANCES(neq,time,FNC,DFNC)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!Right hand side of the population balance equations.
!Depending on flag values, calculation of moments and output messages may also be activated.
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTERFACE

	FUNCTION RATECOAG(i,j,time)
		INTEGER, INTENT (IN) :: i, j
		DOUBLE PRECISION, INTENT (IN) :: time
		DOUBLE PRECISION :: RATECOAG
	END FUNCTION RATECOAG

	FUNCTION MOMENT(nx,ny)
		INTEGER, INTENT(IN) :: nx, ny
		DOUBLE PRECISION :: MOMENT
	END FUNCTION MOMENT

END INTERFACE



INTEGER, INTENT (IN) :: neq
DOUBLE PRECISION, INTENT (IN) :: time, FNC(neq)
DOUBLE PRECISION, INTENT (OUT):: DFNC(neq)


INTEGER :: i, j, n

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


n=1

DO i=1,Mx

	DO j=1,My

		F(i,j)= FNC(n)
		Np(i,j) = F(i,j)*area(i,j)

		n=n+1

	END DO

END DO





!	*** POPULATION BALANCE EQUATIONS ***



IF (flag_output==0) THEN

	n=1

	DO i=1,Mx

		DO j=1,My

			DFNC(n) = RATECOAG(i,j,time)/area(i,j)

			n=n+1

		END DO

	END DO

END IF




!	*** OUTPUT ***


!CALL OUTPUT (2)

IF (flag_output==1) THEN

	mm00 = MOMENT(0,0)
	mm10 = MOMENT(1,0)
	mm01 = MOMENT(0,1)
	mm20 = MOMENT(2,0)
	mm02 = MOMENT(0,2)
	mm11 = MOMENT(1,1)
	mm30 = MOMENT(3,0)
	mm03 = MOMENT(0,3)
	mm21 = MOMENT(2,1)
	mm12 = MOMENT(1,2)
	!mm40 = MOMENT(4,0)
	!mm22 = MOMENT(2,2)
	!mm50 = MOMENT(5,0)
	!mm32 = MOMENT(3,2)
	!mm40 = MOMENT(4,0)


	CALL OUTPUT(3)

END IF



END SUBROUTINE PBEBALANCES
!##############################################################################################










FUNCTION RATECOAG(i,j,time)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This function computes the net rate of coagulation, i.e. birth minus death.
!See Equation 8 of paper.
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTERFACE

	FUNCTION DELTAK(i,j)
		INTEGER, INTENT (IN) :: i,j
		DOUBLE PRECISION :: DELTAK
	END FUNCTION DELTAK

END INTERFACE

INTEGER, INTENT (IN) :: i, j
DOUBLE PRECISION, INTENT (IN) :: time
DOUBLE PRECISION :: RATECOAG

DOUBLE PRECISION :: birth, death, eta
INTEGER :: k, l, q, r, n

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!	***	BIRTH ***

birth = 0.0d0
n=1

DO
	k = cgmx(i,j)%k(n)

	IF (k==0) EXIT

	q = cgmx(i,j)%q(n)
	l = cgmx(i,j)%l(n)
	r = cgmx(i,j)%r(n)

	eta = cgmx(i,j)%eta(n)

	birth = birth + (1.0d0 - 0.5d0*DELTAK(k,q)*DELTAK(l,r))*eta*Bmatrix(k,l,q,r)*Np(k,l)*Np(q,r)

	n=n+1

END DO




!	*** DEATH ***

death=0.0d0


DO k=1,Mx

	DO l=1,My

		death = death + Bmatrix(i,j,k,l)*Np(k,l)

	END DO

END DO

death = death*Np(i,j)



!	*** TOTAL ***

RATECOAG = birth - death


END FUNCTION RATECOAG
!##############################################################################################









FUNCTION DELTAK(i,j)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!Delta Kronecker.
!
!----------------------------------------------------------------------------------------------

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: i, j
DOUBLE PRECISION :: DELTAK

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


IF (i==j) THEN

    DELTAK = 1.0d0

ELSE

    DELTAK = 0.0d0

END IF


END FUNCTION DELTAK
!##############################################################################################









SUBROUTINE JAC

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine is not used. However, it must be declared to run DLSODE.
!
!----------------------------------------------------------------------------------------------

END SUBROUTINE JAC
!##############################################################################################
