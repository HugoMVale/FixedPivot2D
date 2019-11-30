SUBROUTINE GRID(meshx,meshy)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!The construction of the grid is done here. The 2 axes are handled independently. For each axis,
!one may select between:
!   1: linear mesh
!   2: geometric mesh
!   3: dirac (on first cell) + geometric mesh
!
!VARIABLES:
!
!meshx  mesh type for x-axis
!meshy  mesh type for y-axis
!Rx     mesh ratio or width for x-axis
!Ry     mesh ratio or width for y-axis
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT(IN) :: meshx, meshy

DOUBLE PRECISION :: Rx, Ry
INTEGER :: i, j

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!  *** GRID ***


SELECT CASE (meshx)

    CASE (1)

    !	* LINEAR MESH

    Rx = (xmax - xmin)/DBLE(Mx)

    DO i=0,Mx

        xright(i) = xmin + Rx*DBLE(i)

    END DO



    CASE(2)

    !	* GEOMETRIC MESH

    Rx = (xmax/xmin)**(1.0d0/DBLE(Mx))

    DO i=0,Mx

        xright(i) = xmin*Rx**(i)

    END DO



    CASE (3)

    !	* DIRAC + GEOMETRIC MESH

    Rx = (xmax/xmin)**(1.0d0/DBLE(Mx-1))

    xright(0) = 0.0d0

    DO i=1,Mx

        xright(i) = xmin*Rx**(i-1)

    END DO


END SELECT



SELECT CASE (meshy)

    CASE (1)

    !	* LINEAR MESH

    Ry = (ymax - ymin)/DBLE(My)

    DO i=0,My

        yright(i) = ymin + Ry*DBLE(i)

    END DO



    CASE (2)

    !	* GEOMETRIC MESH

    Ry = (ymax/ymin)**(1.0d0/DBLE(My))

    DO i=0,My

        yright(i) = ymin*Ry**(i)

    END DO



    CASE (3)

    !	* DIRAC + GEOMETRIC MESH

    Ry = (ymax/ymin)**(1.0d0/DBLE(My-1))

    yright(0) = 0.0d0

    DO i=1,My

        yright(i) = ymin*Ry**(i-1)

    END DO


END SELECT




!	*** CELL VALUES and SIZES ***


SELECT CASE (meshx)

    CASE (1,2)

        DO i=1,Mx

            x(i) = (xright(i) + xright(i-1))/2.0d0
            dx(i) = xright(i) - xright(i-1)

        END DO

    CASE (3)

        x(1) = xright(0)
        dx(1) = xright(1)

        DO i=2,Mx

            x(i) = (xright(i) + xright(i-1))/2.0d0
            dx(i) = xright(i) - xright(i-1)

        END DO

END SELECT



SELECT CASE (meshy)

    CASE (1,2)

        DO i=1,My

            y(i) = (yright(i) + yright(i-1))/2.0d0
            dy(i) = yright(i) - yright(i-1)

        END DO

    CASE (3)

        y(1) = yright(0)
        dy(1) = yright(1)

        DO i=2,My

            y(i) = (yright(i) + yright(i-1))/2.0d0
            dy(i) = yright(i) - yright(i-1)

        END DO

END SELECT



DO i=1,Mx

    DO j=1,My

        area(i,j) = dx(i)*dy(j)

    END DO

END DO

x(0) = xright(0)
y(0) = yright(0)
x(Mx+1) = xright(Mx)
y(My+1) = yright(My)


END SUBROUTINE GRID
!##############################################################################################
