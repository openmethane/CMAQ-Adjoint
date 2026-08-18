
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

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_fwd/driver_fwd.F,v 1.5 2013/05/15 20:14:11 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)driver_fwd.F  1.2 /project/mod3/CMAQ/src/driver/ctm/SCCS/s.driver_fwd.F 07 Jul 1997 09:20:18

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      PROGRAM DRIVER_FWD

C-----------------------------------------------------------------------
C Function:
C   Do forward simulation and write checkpoint states

C Revision History:
C   Apr 2007 Kumaresh Singh and Adrian Sandu at Virginia Tech: created
C   Aug 2010 by Peter Percell at UH-IMAQS: rewrote for adjoint/4dvar
C   May 2013 by Peter Percell at UH-IMAQS:
C   - support both constant and variable synchronization steps, depending on
C     the choice of checkpointing method
C   - changed the calling argument NADVSS back to the traditional ASTEP
C   - restored output of ACONC
C   - changed XFIRST to an array indexed by layer
C-----------------------------------------------------------------------

      USE PCGRID_DEFN           ! inherits HGRD_DEFN and CGRID_SPCS
      USE VGRD_DEFN             ! vertical layer specifications
      USE WVEL_DEFN             ! derived vertical velocity component
      USE STD_CONC              ! standard CONC
      USE AVG_CONC              ! integral average CONC
      USE ADJOINT_FILES         ! adjoint file data and functions
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines

      USE ADJ_FDM_TEST          ! FDM test

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Load the mechanism COMMON ... (next two INCLUDE files)
      INCLUDE SUBST_RXCMMN      ! chemical mechamism reactions COMMON
      INCLUDE SUBST_RXDATA      ! chemical mechamism reactions DATA

C External Functions:
      INTEGER, EXTERNAL :: GETEFILE

C Local Variables:
      INTEGER, SAVE :: LOGDEV   ! FORTRAN unit number for log file

      CHARACTER(16) :: PNAME = 'DRIVER_FWD'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(2) :: COLROW = 'CR' ! col/row arg list order

      INTEGER, ALLOCATABLE, SAVE :: ASTEP(:) ! advection time steps
      LOGICAL, ALLOCATABLE, SAVE :: XFIRST(:)
      INTEGER :: ALLOCSTAT

      REAL, POINTER :: CGRID(:, :, :, :)

      REAL :: CLOCK             ! Wall-clock time (sec) at initialization.
      INTEGER :: PAR_ERR        ! Error code from parallel initialization.

      INTEGER :: JDATE          ! current model date, coded YYYYDDD
      INTEGER :: JTIME          ! current model time, coded HHMMSS

      INTEGER :: SDATE          ! start date, format YYYYDDD
      INTEGER :: STIME          ! start time, format HHMMSS
      INTEGER :: TSTEP(2)       ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
      INTEGER :: NSTEPS         ! run duration: number of output time steps

      INTEGER :: NREPS          ! number of sync steps per output step

      INTEGER :: IUNIT

      INTEGER :: STEP_NUM       ! current output time step number
      INTEGER :: IREP           ! model step number within this output step

      INTEGER :: C, R, L, K, S, V ! loop counters

      REAL, ALLOCATABLE :: AGRID( :,:,:,: )
      REAL :: DIVFAC            ! trapezoidal average factor
      INTEGER A_NLYS

!slz ptb/perturbation
!slz      REAL    :: PTB  !slz
!slz      INTEGER :: PTBV  !slz
!slz      REAL, EXTERNAL :: ENVREAL  !slz
!slz      INTEGER, EXTERNAL :: ENVINT  !slz
!slz      INTEGER :: STATUS  !slz

      INTERFACE
         SUBROUTINE INITSCEN(CGRID, SDATE, STIME, TSTEP, NSTEPS)
            IMPLICIT NONE
            REAL, POINTER        :: CGRID(:, :, :, :)
            INTEGER, INTENT(OUT) :: SDATE, STIME
            INTEGER, INTENT(OUT) :: TSTEP(2)
            INTEGER, INTENT(OUT) :: NSTEPS
         END SUBROUTINE INITSCEN
         SUBROUTINE SYNCSTEP(JDATE, JTIME, TSTEP, NSTEPS, NREPS)
            IMPLICIT NONE
            INTEGER, INTENT(IN)    :: JDATE, JTIME
            INTEGER, INTENT(INOUT) :: TSTEP(2)
            INTEGER, INTENT(IN)    :: NSTEPS
            INTEGER, INTENT(OUT)   :: NREPS
         END SUBROUTINE SYNCSTEP
         SUBROUTINE ADVSTEP(JDATE, JTIME, TSTEP, ASTEP)
            IMPLICIT NONE
            INTEGER, INTENT(IN)  :: JDATE, JTIME
            INTEGER, INTENT(IN)  :: TSTEP(2)
            INTEGER, INTENT(OUT) :: ASTEP(:)
         END SUBROUTINE ADVSTEP
         SUBROUTINE SCIPROC(CGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)
            IMPLICIT NONE
            REAL, POINTER          :: CGRID(:, :, :, :)
            INTEGER, INTENT(INOUT) :: JDATE, JTIME
            INTEGER, INTENT(IN)    :: TSTEP(2)
            INTEGER, INTENT(IN)    :: ASTEP(:)
            LOGICAL, INTENT(INOUT) :: XFIRST(:)
         END SUBROUTINE SCIPROC
         SUBROUTINE WR_ACONC(AGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL,    INTENT(IN) :: AGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME, TSTEP
         END SUBROUTINE WR_ACONC
         SUBROUTINE WR_CGRID(CGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: CGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME, TSTEP
         END SUBROUTINE WR_CGRID
      END INTERFACE

C-----------------------------------------------------------------------

      JDATE = 0
      JTIME = 0

C Get number of species, and starting indices for CGRID array.
      CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

C Start up processor communication and define horizontal domain decomposition
C and vertical layer structure
      CALL PAR_INIT(COLROW, NSPCSD, CLOCK, PAR_ERR)
      IF ( PAR_ERR /= 0 ) THEN
         XMSG = 'Error in PAR_INIT'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      LOGDEV = INIT3()

      IF ( NSPCSD .GT. MXVARS3 ) THEN
         WRITE(XMSG,'(5X, A, I5, A)') 'The number of variables,', NSPCSD,
     &        ' to be written to the State CGRID File'
         WRITE(LOGDEV, '(A)') XMSG
         WRITE(XMSG,'(5X, A, I5)') 'exceeds the I/O-API limit:', MXVARS3
         WRITE(LOGDEV, '(A)') XMSG
         XMSG = 'Recompile with an I/O-API lib having a larger MXVARS3'
         WRITE(LOGDEV, '(5X, A)') XMSG
         CALL M3EXIT(PNAME, JDATE, JTIME, ' ', XSTAT1)
      END IF

C Initialize PCGRID
      IF ( .NOT. PCGRID_INIT() ) THEN
         XMSG = 'Failure defining horizontal domain'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
      END IF

      CGRID => PCGRID(1:NCOLS, 1:NROWS, :, :)

C Initalize CONC definitions (in STD_CONC F90 module)
      CALL CONC_DEFN()

C Get avg CONC definitions, species and layer pointers (in AVG_CONC F90 module)
      CALL A_CONC_DEFN()
      A_NLYS = ACONC_ELEV - ACONC_BLEV + 1

C Allocate memory
      ALLOCATE(
     &     ASTEP(1:NLAYS),
     &     XFIRST(1:NLAYS),
     &     AGRID(NCOLS, NROWS, A_NLYS, N_ASPCS),
     &     STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating ASTEP, XFIRST or AGRID'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF
      AGRID = 0.0

C Initialize optional derived vertical velocity writes to conc file
      IF ( .NOT. WVEL_INIT() ) THEN
         XMSG = 'Failure initializing derived vertical velocity writes'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2 )
      END IF

C Initialize conc field: Copy IC's to CONC file as step 0
C Convention: the input file concentration units are always ppmV.
      CALL INITSCEN(CGRID, SDATE, STIME, TSTEP, NSTEPS)

C Verify input file header consistency with COORD.EXT and run duration
      CALL FLCHECK(SDATE, STIME, TSTEP(1))

C If the sync step must be constant for this run, get a synchronization step,
C TSTEP(2), and number of repetitions per output step, NREPS, for the whole run
C (Note: For IOAPI compatibility, TSTEP(2) and NREPS must be constant for a run
C because checkpoints are written every sync step and the IOAPI tools require a
C constant output time step for files.)
      IF ( CONST_SYNC_STEP )
     &     CALL SYNCSTEP(SDATE, STIME, TSTEP, NSTEPS, NREPS)

C Read XFIRST from the file CTM_XFIRST_IN
      IUNIT = GETEFILE(CTM_XFIRST_IN, .TRUE., .TRUE., PNAME)
      IF ( IUNIT > 0 ) THEN
         READ(IUNIT, '(L1)') XFIRST(1:NLAYS)
         CLOSE(IUNIT)
      ELSE
         XFIRST(1:NLAYS) = .TRUE.
      END IF

C Set the current date and time to the start date and time
      JDATE = SDATE
      JTIME = STIME

!slz ptb/perturbation
      CALL FDM_INIT(JDATE, JTIME)
      IF(FDM_TEST) CREATE_CHK=.FALSE.

C Do forward simulation and write checkpoint states
      DO STEP_NUM = 1, NSTEPS   ! output time step loop

C If the sync step does not need to be constant for this run, get a
C synchronization step, TSTEP(2), and number of repetitions per output
C step, NREPS, just for this output step
         IF ( .NOT. CONST_SYNC_STEP )
     &        CALL SYNCSTEP(JDATE, JTIME, TSTEP, 1, NREPS)

C Get the advection time steps, ASTEP(L)
         CALL ADVSTEP(JDATE, JTIME, TSTEP, ASTEP)

         DO V = 1, N_ASPCS
            S = AVG_CONC_MAP(V)
            L = 0
            DO K = ACONC_BLEV, ACONC_ELEV
               L = L + 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     AGRID(C, R, L, V) = CGRID(C, R, K, S)
                  END DO
               END DO
            END DO
         END DO

         DIVFAC = 0.5 / FLOAT(NREPS)
         DO IREP = 1, NREPS - 1
!slz perturbation
            IF ( FDM_TEST.AND.PTB_TIME.EQ.JTIME ) CALL FDM_PTB(CGRID, JDATE, JTIME)
            CALL SCIPROC(CGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST) ! science processes

            DO V = 1, N_ASPCS
               S = AVG_CONC_MAP(V)
               L = 0
               DO K = ACONC_BLEV, ACONC_ELEV
                  L = L + 1
                  DO R = 1, MY_NROWS
                     DO C = 1, MY_NCOLS
                        AGRID(C, R, L, V) = AGRID(C, R, L, V)
     &                       + 2.0 * CGRID(C, R, K, S)
                     END DO
                  END DO
               END DO
            END DO
         END DO

         IF ( FDM_TEST.AND.PTB_TIME.EQ.JTIME ) CALL FDM_PTB(CGRID, JDATE, JTIME)
         CALL SCIPROC(CGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)

         DO V = 1, N_ASPCS
            S = AVG_CONC_MAP(V)
            L = 0
            DO K = ACONC_BLEV, ACONC_ELEV
               L = L + 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     AGRID(C, R, L, V) = DIVFAC * ( AGRID(C, R, L, V)
     &                                              + CGRID(C, R, K, S) )
                  END DO
               END DO
            END DO
         END DO

         DO V = 1, N_CSPCS
            S = CONC_MAP(V)
            L = 0
            DO K = CONC_BLEV, CONC_ELEV
               L = L + 1
               DO R = 1, NROWS
                  DO C = 1, NCOLS
                     SGRID(C, R, L, V) = CGRID(C, R, K, S)
                  END DO
               END DO
            END DO
         END DO

         CALL WR_CONC(JDATE, JTIME, TSTEP(1))

         CALL WR_ACONC(AGRID, JDATE, JTIME, TSTEP(1))

      END DO

      CALL WR_CGRID(CGRID, JDATE, JTIME, TSTEP(1))

C Write the final values of XFIRST to the file CTM_XFIRST_OUT
      IF ( MYPE == 0 ) THEN
         IUNIT = GETEFILE(CTM_XFIRST_OUT, .FALSE., .TRUE., PNAME)
         WRITE(IUNIT, '(L1)') XFIRST(1:NLAYS)
         CLOSE(IUNIT)
      END IF

C Close checkpoint files
      CALL CHK_FILE_CLOSE

C Clean up communications
      CALL PAR_TERM(PNAME, JDATE, JTIME, CLOCK)

      END PROGRAM DRIVER_FWD
