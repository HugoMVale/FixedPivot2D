SUBROUTINE OUTPUT(message)

!----------------------------------------------------------------------------------------------
!
!DESCRIPTION:
!
!This subroutine manages all outputs to console and files.
!The file paths must be adapted depending on program location.
!
!VARIABLES:
!
!Fdiag(i)   average number density function along grid diagonal, i.e. in cell Cij
!nevals     number of function calls
!
!----------------------------------------------------------------------------------------------

USE GLOBAL

IMPLICIT NONE

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

INTEGER, INTENT (IN) :: message

DOUBLE PRECISION :: Fdiag(1:Mx)
INTEGER :: i, nevals=0

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


SELECT CASE (message)

! Start of simulation
CASE (1)

    PRINT *, "RUNNING SIMULATION..."
    PRINT *
    PRINT *, FDATE()
    PRINT *
    PRINT *

    OPEN (UNIT=1, FILE="C:\FortranProjects\FixedPivot2D\code\F.dat", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    OPEN (UNIT=2, FILE="C:\FortranProjects\FixedPivot2D\code\Moments.dat", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")
    OPEN (UNIT=3, FILE="C:\FortranProjects\FixedPivot2D\code\Fdiag.dat", STATUS="OLD", ACTION="WRITE", POSITION="REWIND")

    PRINT '(2X,A6,4(A12))', "Time", "M_00", "M_10", "M_01", "M_11"

    !Header file 1

    WRITE (1,'(1X,A15)', ADVANCE="NO") "x(i)"

    DO i=1,Mx

        IF (i<Mx) THEN
            WRITE (1,'(A15))', ADVANCE="NO") "y(j)/F(i,j)"
        ELSE
            WRITE (1,'(A15))') "y(j)/F(i,j)"
        END IF

    END DO

    !Header file 2

    WRITE (2,'(1X,11(A15))') "time", "mm00", "mm10", "mm01", "mm20", "mm02", "mm11", "mm30", "mm03", "mm12", "mm21"

    !Header file 3

    WRITE (3,'(1X,A15)', ADVANCE="NO") "time"

    DO i=1,Mx

        IF (i<Mx) THEN
            WRITE (3,'(A15))', ADVANCE="NO") "Fdiag"
        ELSE
            WRITE (3,'(A15))') "Fdiag"
        END IF

    END DO

! All intermediate function calls
CASE (2)

    nevals = nevals + 1

    PRINT '(1X,A,I7)',"dF/dt Evals.: ", nevals
    PRINT *

! Predefined time steps
CASE (3)

    !CALL BEEPQQ(1000,150)

    PRINT '(2X,F6.2,4(E12.4))', time_now, mm00, mm10, mm01, mm11
    PRINT *


    DO i=1,Mx

        Fdiag(i) = F(i,i)

    END DO

    WRITE (2,'(1X,11(E15.5))') time_now, mm00, mm10, mm01, mm20, mm02, mm11, mm30, mm03, mm12, mm21
    WRITE (3,'(1X,51(E15.5))') time_now, (Fdiag(1:Mx) + 1.0d-50)



    IF (last==1) THEN

        WRITE (1,'(1X,51(E15.5))') 0.0d0, y(1:My)

        DO i=1,Mx

            WRITE (1,'(1X,51(E15.5))') x(i), (DABS(F(i,:))+1.0d-50)

        END DO

    END IF

! End of simulation
CASE (4)

    CLOSE (1)
    CLOSE (2)

    PRINT *
    PRINT *," *** END OF SIMULATION ***"
    PRINT *
    PRINT *, FDATE()
    PRINT *

END SELECT



END SUBROUTINE OUTPUT
!##############################################################################################
