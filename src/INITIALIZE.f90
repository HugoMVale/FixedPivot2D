SUBROUTINE AVERAGECELL

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!The average cell value, F(i,j) is computed by definition, i.e., by evaluating the
!integral of F over the cell domain. The multiple integral is computed as an iterated
!integral; Simpson's 1/3 rule is applied in both directions.
!
!VARIABLES:
!
!xf     quadrature points along x-axis
!yf     quadrature points along y-axis
!z      function values at the quadrature points
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTERFACE

    FUNCTION PSDINI(xf,yf)
        DOUBLE PRECISION, INTENT (IN) :: xf, yf
        DOUBLE PRECISION :: PSDINI
    END FUNCTION PSDINI

END INTERFACE


INTEGER, PARAMETER :: npoints=9
DOUBLE PRECISION :: xf(npoints), yf(npoints), z(npoints)
INTEGER :: i, j, n

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


DO i=1,Mx

    DO j=1,My

        xf(1)=xright(i-1)
        yf(1)=yright(j-1)

        xf(2)=x(i)
        yf(2)=yright(j-1)

        xf(3)=xright(i)
        yf(3)=yright(j-1)

        xf(4)=xright(i-1)
        yf(4)=y(j)

        xf(5)=x(i)
        yf(5)=y(j)

        xf(6)=xright(i)
        yf(6)=y(j)

        xf(7)=xright(i-1)
        yf(7)=yright(j)

        xf(8)=x(i)
        yf(8)=yright(j)

        xf(9)=xright(i)
        yf(9)=yright(j)

        DO n=1,npoints

            z(n) = PSDINI(xf(n),yf(n))

        END DO


        F(i,j)= 1.0d0/36.0d0*((z(1) + 4.0d0*z(2) + z(3)) &

                + 4.0d0*(z(4) + 4.0d0*z(5) + z(6)) + (z(7) + 4.0d0*z(8) + z(9)))

    END DO

END DO


END SUBROUTINE AVERAGECELL
!##############################################################################################








SUBROUTINE COAGMATRIX

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!Here we (pre)compute coagulation kernel for all possible particle combinations and store result
!in a coagulation matrix. Only makes sense if the coagulation kernel is time-independent.
!
!See Equation 9 of paper.
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTERFACE

	FUNCTION COAGKERNEL(x1,y1,x2,y2,time)
		DOUBLE PRECISION, INTENT (IN) :: x1, y1, x2, y2, time
		DOUBLE PRECISION :: COAGKERNEL
	END FUNCTION COAGKERNEL

END INTERFACE


INTEGER :: k, q, l, r

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!	*** COMPUTE MATRIX WITH B VALUES ***

PRINT *
PRINT *, "Computing B values..."
PRINT *


Bmatrix = 0.0d0

DO k=1,Mx

	DO q=1,Mx

		DO l=1,My

			DO r=1,My

				Bmatrix(k,l,q,r)=COAGKERNEL(x(k),y(l),x(q),y(r),0.0d0)

			END DO

		END DO

	END DO

END DO


END SUBROUTINE COAGMATRIX
!##############################################################################################
