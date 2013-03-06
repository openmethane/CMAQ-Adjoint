
C***********************************************************************
C   Portions of Models-3/CMAQ software were developed or based on      *
C   information from various groups: Federal Government employees,     *
C   contractors working on a United States Government contract, and    *
C   non-Federal sources (including research institutions).  These      *
C   research institutions have given the Government permission to      *
C   use, prepare derivative works, and distribute copies of their      *
C   work in Models-3/CMAQ to the public and to permit others to do     *
C   so.  EPA therefore grants similar permissions for use of the       *
C   Models-3/CMAQ software, but users are requested to provide copies  *
C   of derivative works to the Government without restrictions as to   *
C   use by others.  Users are responsible for acquiring their own      *
C   copies of commercial software associated with Models-3/CMAQ and    *
C   for complying with vendor requirements.  Software copyrights by    *
C   the MCNC Environmental Modeling Center are used with their         *
C   permissions subject to the above restrictions.                     *
C***********************************************************************

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SCIPROC_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS, XFIRST)

C-----------------------------------------------------------------------
C Function:
C   Science process adjoint for sensitivity analysis

C Revision History:
C   Apr 2007 Kumaresh Singh and Adrian Sandu at Virginia Tech: created
C   Aug 2010 Peter Percell at UH-IMAQS: rewrote for adjoint/4dvar
C-----------------------------------------------------------------------

      IMPLICIT NONE

C Include files:
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

C Arguments:
      REAL, POINTER :: LGRID(:,:,:,:) ! adjoint variable
      INTEGER, INTENT(INOUT) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(INOUT) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                      !   TSTEP(1) = local output step
                                      !   TSTEP(2) = sciproc sync. step (chem)
      INTEGER, INTENT(IN) :: NADVSS(:) ! number of advection sub-steps per sync step
      LOGICAL, INTENT(INOUT) :: XFIRST

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      CHARACTER(16) :: PNAME = 'SCIPROC_BWD'
      CHARACTER(36) :: NMSG = 'After NEXTIME: returned JDATE, JTIME'

      INTERFACE
         SUBROUTINE AERO_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE AERO_BWD
         SUBROUTINE CHEM_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE CHEM_BWD
         SUBROUTINE CLDPROC_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE CLDPROC_BWD
         SUBROUTINE HDIFF_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE HDIFF_BWD
         SUBROUTINE ZADV_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE ZADV_BWD
         SUBROUTINE XADV_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
            INTEGER, INTENT(IN) :: NADVSS(:)
         END SUBROUTINE XADV_BWD
         SUBROUTINE YADV_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
            INTEGER, INTENT(IN) :: NADVSS(:)
         END SUBROUTINE YADV_BWD
         SUBROUTINE VDIFF_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:,:,:,:)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE VDIFF_BWD
c$$$         SUBROUTINE ADJADV_BWD(LGRID, JDATE, JTIME, TSTEP)
c$$$            IMPLICIT NONE
c$$$            REAL, POINTER       :: LGRID(:,:,:,:)
c$$$            INTEGER, INTENT(IN) :: JDATE, JTIME
c$$$            INTEGER, INTENT(IN) :: TSTEP(2)
c$$$         END SUBROUTINE ADJADV_BWD
c$$$         SUBROUTINE COUPLE(LGRID, JDATE, JTIME, TSTEP)
c$$$            IMPLICIT NONE
c$$$            REAL, POINTER       :: LGRID(:,:,:,:)
c$$$            INTEGER, INTENT(IN) :: JDATE, JTIME
c$$$            INTEGER, INTENT(IN) :: TSTEP(2)
c$$$         END SUBROUTINE COUPLE
c$$$         SUBROUTINE DECOUPLE(LGRID, JDATE, JTIME, TSTEP)
c$$$            IMPLICIT NONE
c$$$            REAL, POINTER       :: LGRID(:,:,:,:)
c$$$            INTEGER, INTENT(IN) :: JDATE, JTIME
c$$$            INTEGER, INTENT(IN) :: TSTEP(2)
c$$$         END SUBROUTINE DECOUPLE
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
      END IF

      CALL NEXTIME(JDATE, JTIME, -TSTEP(2))
      WRITE(LOGDEV, '(/ 5X, A, I8, I7.6)') NMSG, JDATE, JTIME

      CALL AERO_BWD(LGRID, JDATE, JTIME, TSTEP)

      CALL CHEM_BWD(LGRID, JDATE, JTIME, TSTEP)

      CALL CLDPROC_BWD(LGRID, JDATE, JTIME, TSTEP)

c$$$      CALL DECOUPLE(LGRID, JDATE, JTIME, TSTEP)

      CALL HDIFF_BWD(LGRID, JDATE, JTIME, TSTEP)

c$$$      CALL ADJADV_BWD(LGRID, JDATE, JTIME, TSTEP)

      CALL ZADV_BWD(LGRID, JDATE, JTIME, TSTEP)

      IF ( XFIRST ) THEN
         XFIRST = .FALSE.

         CALL XADV_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS)

         CALL YADV_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS)
      ELSE
         XFIRST = .TRUE.

         CALL YADV_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS)

         CALL XADV_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS)
      END IF

c$$$      CALL COUPLE(LGRID, JDATE, JTIME, TSTEP)

      CALL VDIFF_BWD(LGRID, JDATE, JTIME, TSTEP)

      RETURN
      END
