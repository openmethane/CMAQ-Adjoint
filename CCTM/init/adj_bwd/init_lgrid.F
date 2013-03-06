
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
      SUBROUTINE INIT_LGRID(LGRID, SDATE, STIME, EDATE, ETIME, TSTEP, NSTEPS)

C-----------------------------------------------------------------------
C Function:
C   Determine the simulation time period and time stepping constants;
C   initialize the LGRID array.

C Revision History:
C   Aug 2010 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C-----------------------------------------------------------------------

      USE ADJOINT_FILES         ! adjoint file data and functions

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_FILES_ID    ! file name parameters

C External Functions (not already declared by IODECL3.EXT):
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC

C Arguments:
      REAL, POINTER :: LGRID(:, :, :, :)
      INTEGER, INTENT(OUT) :: SDATE ! start date, format YYYYDDD
      INTEGER, INTENT(OUT) :: STIME ! start time, format HHMMSS
      INTEGER, INTENT(OUT) :: EDATE ! end date, format YYYYDDD
      INTEGER, INTENT(OUT) :: ETIME ! end time, format HHMMSS
      INTEGER, INTENT(OUT) :: TSTEP(2) ! time step vector (HHMMSS)
                                       !   TSTEP(1) = local output step
                                       !   TSTEP(2) = sciproc sync. step (chem)
      INTEGER, INTENT(OUT) :: NSTEPS ! run duration: number of output time steps

C Local Variables:
      CHARACTER(16) :: PNAME = 'INIT_LGRID'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: JDATE          ! current model date, coded YYYYDDD
      INTEGER :: JTIME          ! current model time, coded HHMMSS

      INTERFACE
         SUBROUTINE DEFINE_RECEPTOR(LGRID)
            IMPLICIT NONE
            REAL, POINTER :: LGRID(:, :, :, :)
         END SUBROUTINE DEFINE_RECEPTOR
      END INTERFACE

C-----------------------------------------------------------------------

      JDATE = 0
      JTIME = 0

C Get file characteristics of CTM_CONC_FWD
      IF ( .NOT. OPEN3(CTM_CONC_FWD, FSREAD3, PNAME) ) THEN
         XMSG = 'Could not open ' // TRIM(CTM_CONC_FWD) //
     &        ' file as OLD:READ-ONLY'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      IF ( .NOT. DESC3(CTM_CONC_FWD) ) THEN
         XMSG = 'Could not get ' // TRIM(CTM_CONC_FWD) // ' file description'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      IF ( .NOT. CLOSE3(CTM_CONC_FWD) ) THEN
         XMSG = 'Could not close ' // TRIM(CTM_CONC_FWD)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

C Set the output time step, TSTEP(1), and number of output steps, NSTEPS, to
C the values used for the forward simulation
      TSTEP(1) = TSTEP3D
      NSTEPS = MXREC3D - 1

C Set start date and time and end date and time
      SDATE = SDATE3D
      STIME = STIME3D

      EDATE = SDATE3D
      ETIME = STIME3D
      CALL NEXTIME(EDATE, ETIME, SEC2TIME(NSTEPS * TIME2SEC(TSTEP(1))))

C Initialize LGRID
      IF ( .NOT. OPEN3(INIT_LGRID_1, FSREAD3, PNAME) ) THEN
         XMSG = 'Could not open ' // TRIM(INIT_LGRID_1) //
     &        '; initializing LGRID using DEFINE_RECEPTOR'
         CALL M3WARN(PNAME, EDATE, ETIME, XMSG)

         CALL DEFINE_RECEPTOR(LGRID)
      ELSE
         CALL RD_ADJ_FILE(LGRID, EDATE, ETIME, INIT_LGRID_1)
      END IF

      END SUBROUTINE INIT_LGRID
