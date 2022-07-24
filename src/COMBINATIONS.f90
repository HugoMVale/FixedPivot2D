SUBROUTINE COMBINATIONS

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine is one of the core components of the algorithm. The possible (non-repeated)
!combinations of cells Ckl and Cqr leading to the formation of a new particle in region Rij
!are computed and stored. The values of eta (i.e., the weights) are also computed.
!
!This information occupies a lot of memory and, thus, is stored in a variable with a special
!data structure. It is hard to accurately estimate the necessary storage space before computing
!the combinations. Increase "dim_cgmx" if required.
!
!See Equations 4 and 8 of paper.
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, PARAMETER :: dim_cgmx=10*Mx*My

DOUBLE PRECISION :: xx(Mx,Mx,My,My), yy(Mx,Mx,My,My), xf, yf
INTEGER :: coagmatrix(4,dim_cgmx)
INTEGER :: i, j, k, l, q, r, n, s

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


DO k=1,Mx

    DO q=1,Mx

        DO l=1,My

            DO r=1,My

				xx(k,q,l,r) = x(k) + x(q)

				yy(k,q,l,r) = y(l) + y(r)

                !Uncomment if used for systems with total mass and composition as internal coordinates
				!yy(k,q,l,r)= (x(k)*y(l) + x(q)*y(r))/(x(k) + x(q))

			END DO

		END DO

	END DO

END DO






DO j=1,My

	DO i=1,Mx

		coagmatrix = 0

		n=1

		DO l=1,j

			DO k=1,i

				DO r=l,j ! My !j

					IF (r==l) THEN

						s=k

					ELSE

						s=1

					END IF

					DO q=s,i

						xf = xx(k,q,l,r)
						yf = yy(k,q,l,r)

						IF (x(i-1)<xf .AND. (xf<x(i+1) .AND. (y(j-1)<yf .AND. yf<y(j+1)))) THEN

							coagmatrix(1,n) = k
							coagmatrix(2,n) = q
							coagmatrix(3,n) = l
							coagmatrix(4,n) = r

							n=n+1

						END IF

						IF (n>dim_cgmx)	THEN

							PRINT *, "dim_cgmax too small."
							STOP

						END IF

					END DO

				END DO

			END DO

		END DO

		coagmatrix(:,n) = 0

		ALLOCATE (cgmx(i,j)%k(1:n))
		ALLOCATE (cgmx(i,j)%q(1:n))
		ALLOCATE (cgmx(i,j)%l(1:n))
		ALLOCATE (cgmx(i,j)%r(1:n))
		ALLOCATE (cgmx(i,j)%eta(1:n))

		cgmx(i,j)%k(1:n) = coagmatrix(1,1:n)
		cgmx(i,j)%q(1:n) = coagmatrix(2,1:n)
		cgmx(i,j)%l(1:n) = coagmatrix(3,1:n)
		cgmx(i,j)%r(1:n) = coagmatrix(4,1:n)

	END DO

END DO





DO i=1,Mx

	DO j=1,My

		n=1

		DO

			k = cgmx(i,j)%k(n)

			IF (k==0) EXIT

			q = cgmx(i,j)%q(n)
			l = cgmx(i,j)%l(n)
			r = cgmx(i,j)%r(n)

			xf = xx(k,q,l,r)
			yf = yy(k,q,l,r)


			IF (x(i-1)<=xf .AND. (xf<=x(i) .AND. (y(j)<=yf .AND. yf<=y(j+1)))) THEN

				cgmx(i,j)%eta(n) = WEIGHT(2,x(i-1),x(i),y(j),y(j+1))

			ELSE IF (x(i)<=xf .AND. (xf<=x(i+1) .AND. (y(j)<=yf .AND. yf<=y(j+1)))) THEN

				cgmx(i,j)%eta(n) = WEIGHT(1,x(i),x(i+1),y(j),y(j+1))

			ELSE IF (x(i-1)<=xf .AND. (xf<=x(i) .AND. (y(j-1)<=yf .AND. yf<=y(j)))) THEN

				cgmx(i,j)%eta(n) = WEIGHT(4,x(i-1),x(i),y(j-1),y(j))

			ELSE IF (x(i)<=xf .AND. (xf<=x(i+1) .AND. (y(j-1)<=yf .AND. yf<=y(j)))) THEN

				cgmx(i,j)%eta(n) = WEIGHT(3,x(i),x(i+1),y(j-1),y(j))


			END IF


			n=n+1

		END DO

	END DO

END DO


!##############################################################################################


CONTAINS

FUNCTION WEIGHT(i,x1,x2,y1,y2)

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: i
DOUBLE PRECISION, INTENT (IN) :: x1, x2, y1, y2
DOUBLE PRECISION :: WEIGHT

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

WEIGHT = 0.0d0

SELECT CASE (i)

	CASE (1)

		WEIGHT=(xf-x2)*(yf-y2)/((x1-x2)*(y1-y2))

	CASE (2)

		WEIGHT=(xf-x1)*(yf-y2)/((x1-x2)*(y2-y1))

	CASE (3)

		WEIGHT=(xf-x2)*(yf-y1)/((x2-x1)*(y1-y2))

	CASE (4)

		WEIGHT=(xf-x1)*(yf-y1)/((x1-x2)*(y1-y2))

END SELECT

END FUNCTION WEIGHT


END SUBROUTINE COMBINATIONS
!##############################################################################################
