
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_bwd/driver_bwd.F,v 1.2 2013/05/15 20:14:11 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)driver.F  1.2 /project/mod3/CMAQ/src/driver/ctm/SCCS/s.driver.F 07 Jul 1997 09:20:18

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      PROGRAM DRIVER_BWD

C-----------------------------------------------------------------------
C Function:
C   Do backward run to compute adjoint sensitivities

C Revision History:
C   Apr 2007 Kumaresh Singh and Adrian Sandu at Virginia Tech: created
C   Aug 2010 by Peter Percell at UH-IMAQS: rewrote for adjoint/4dvar
C   May 2013 by Peter Percell at UH-IMAQS:
C   - support both constant and variable synchronization steps, depending on
C     the choice of checkpointing method
C   - changed the calling argument NADVSS back to the traditional ASTEP
C   - changed XFIRST to an array indexed by layer
C   Aug 2013 by Peter Percell at UH-IMAQS:
C   - made the LGRID output frequency selectable by user
C   Nov 2013 S.Capps (US EPA): shifted emissions-related sensitivity
C            arrays to module
C-----------------------------------------------------------------------

      USE PCGRID_DEFN           ! inherits HGRD_DEFN and CGRID_SPCS
      USE VGRD_DEFN             ! vertical layer specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE ADJOINT_FILES         ! adjoint file data and functions
      USE ADJ_ARRAY_MOD         ! adjoint sensitivity arrays
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines

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

      CHARACTER(16) :: PNAME = 'DRIVER_BWD'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(2) :: COLROW = 'CR' ! col/row arg list order

      INTEGER, ALLOCATABLE, SAVE :: ASTEP(:) ! advection time steps
      LOGICAL, ALLOCATABLE, SAVE :: XFIRST(:)
      REAL, POINTER :: LGRID(:, :, :, :)
      
      INTEGER :: ALLOCSTAT

      REAL :: CLOCK             ! Wall-clock time (sec) at initialization.
      INTEGER :: PAR_ERR        ! Error code from parallel initialization.

      INTEGER :: JDATE          ! current model date, coded YYYYDDD
      INTEGER :: JTIME          ! current model time, coded HHMMSS

      INTEGER :: SDATE          ! start date, format YYYYDDD
      INTEGER :: STIME          ! start time, format HHMMSS
      INTEGER :: EDATE          ! end date, format YYYYDDD
      INTEGER :: ETIME          ! end time, format HHMMSS
      INTEGER :: TSTEP(2)       ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
      INTEGER :: NSTEPS         ! run duration: number of output time steps

      INTEGER :: NREPS          ! number of sync steps per output step

      INTEGER :: IUNIT

      INTEGER :: STEP_NUM       ! current output time step number
      INTEGER :: IREP           ! model step number within this output step

      INTERFACE
         SUBROUTINE INIT_LGRID(LGRID, SDATE, STIME, EDATE, ETIME, TSTEP, NSTEPS)
            IMPLICIT NONE
            REAL, POINTER :: LGRID(:, :, :, :)
            INTEGER, INTENT(OUT) :: SDATE
            INTEGER, INTENT(OUT) :: STIME
            INTEGER, INTENT(OUT) :: EDATE
            INTEGER, INTENT(OUT) :: ETIME
            INTEGER, INTENT(OUT) :: TSTEP(2)
            INTEGER, INTENT(OUT) :: NSTEPS
         END SUBROUTINE INIT_LGRID
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
         SUBROUTINE SCIPROC_BWD(LGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)
            IMPLICIT NONE
            REAL, POINTER          :: LGRID(:, :, :, :)
            INTEGER, INTENT(INOUT) :: JDATE, JTIME
            INTEGER, INTENT(IN)    :: TSTEP(2)
            INTEGER, INTENT(IN)    :: ASTEP(:)
            LOGICAL, INTENT(INOUT) :: XFIRST(:)
         END SUBROUTINE SCIPROC_BWD
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
     &        ' to be written to the LGRID File'
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

C Allocate memory
      ALLOCATE(
     &     LGRID(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD),
     &     ASTEP(1:NLAYS),
     &     XFIRST(1:NLAYS),
     &     STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating LGRID, ASTEP, XFIRST'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

C Determine the simulation time period and time stepping constants;
C initialize the LGRID array.
      LGRID(:, :, :, :) = 0.0
      LGRID_FRC(:, :, :, :) = 0.0
      LGRID_FRC_TOT(:, :, :, :) = 0.0
      CALL INIT_LGRID(LGRID, SDATE, STIME, EDATE, ETIME, TSTEP, NSTEPS)

C Verify input file header consistency with COORD.EXT and run duration
      CALL FLCHECK(SDATE, STIME, TSTEP(1))

C If the sync step must be constant for this run, get a synchronization step,
C TSTEP(2), and number of repetitions per output step, NREPS, for the whole run
C (Note: For IOAPI compatibility, TSTEP(2) and NREPS must be constant for a run
C because checkpoints are written every sync step and the IOAPI tools require a
C constant output time step for files.)
      IF ( CONST_SYNC_STEP )
     &     CALL SYNCSTEP(SDATE, STIME, TSTEP, NSTEPS, NREPS)

C Create new adjoint output files and open for update
      CALL OP_ADJ_FILE(EDATE, ETIME, TSTEP, ADJ_LGRID)
      ! Also for emissions related sensitivity output (mdt|slc.11.2013)
      CALL OP_ADJ_FILE(EDATE, ETIME, TSTEP, ADJ_LGRID_EM)
      CALL OP_ADJ_FILE(EDATE, ETIME, TSTEP, ADJ_LGRID_EM_SF)
      CALL OP_ADJ_FILE(EDATE, ETIME, TSTEP, ADJ_LGRID_EM_NRM)

      CALL WR_ADJ_FILE(LGRID, EDATE, ETIME, ADJ_LGRID)

      IF ( TRIM(LGRID_OUTPUT_FREQ) /= 'ONCE' )
     &     CALL WR_ADJ_FILE(LGRID, EDATE, ETIME, ADJ_LGRID)

C Read XFIRST from the file CTM_XFIRST_IN
      IUNIT = GETEFILE(CTM_XFIRST_IN, .TRUE., .TRUE., PNAME)
      IF ( IUNIT > 0 ) THEN
         READ(IUNIT, '(L1)') XFIRST(1:NLAYS)
         CLOSE(IUNIT)
      ELSE
         XMSG = 'The XFIRST definition file does not exist'
         CALL M3EXIT(PNAME, EDATE, ETIME, XMSG, XSTAT1)
      END IF

C Set the current date and time to the end date and time
      JDATE = EDATE
      JTIME = ETIME

C Do backward (adjoint) solution and write sensitivities (step backwards from
C final time)
      DO STEP_NUM = NSTEPS, 1, -1 ! output time step loop
         SDATE = JDATE
         STIME = JTIME
         CALL NEXTIME(SDATE, STIME, -TSTEP(1))

C If the sync step does not need to be constant for this run, get a
C synchronization step, TSTEP(2), and number of repetitions per output
C step, NREPS, just for this output step
         IF ( .NOT. CONST_SYNC_STEP )
     &        CALL SYNCSTEP(SDATE, STIME, TSTEP, 1, NREPS)

C Get the advection time steps, ASTEP(L)
         CALL ADVSTEP(SDATE, STIME, TSTEP, ASTEP)

C Read adjoint forcing from file
         IF ( .NOT. OPEN3(ADJ_FORCE, FSREAD3, PNAME) ) THEN
            XMSG = 'Could not open ' // TRIM(ADJ_FORCE) // 
     &             ' file; exiting'
            CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT1)
         ELSE
            CALL RD_ADJ_FILE(LGRID_FRC_TOT, SDATE, STIME, ADJ_FORCE)
         END IF

C Since offline forcing is based on concentrations averaged over the output time step,
C   the forcing can simply be allocated evenly to each synchronization time step. 11.2013.slc
C        LGRID_FRC is in ADJ_ARRAY_MOD
         LGRID_FRC = LGRID_FRC_TOT / NSTEPS
         
         DO IREP = NREPS, 1, -1
            CALL SCIPROC_BWD(LGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)                     ! science processes
            IF ( TRIM(LGRID_OUTPUT_FREQ) == 'SYNC_STEP' ) THEN
                 CALL WR_ADJ_FILE(LGRID, JDATE, JTIME, ADJ_LGRID)
                 CALL WR_ADJ_FILE(LGRID_EM, JDATE, JTIME, 
     &                            ADJ_LGRID_EM)
                 CALL WR_ADJ_FILE(LGRID_EM_SF, JDATE, JTIME, 
     &                            ADJ_LGRID_EM_SF)
                 CALL WR_ADJ_FILE(LGRID_EM_NRM, JDATE, JTIME, 
     &                            ADJ_LGRID_EM_NRM)
            END IF
         END DO
         
         IF ( TRIM(LGRID_OUTPUT_FREQ) == 'OUTPUT_STEP' ) THEN
            CALL WR_ADJ_FILE(LGRID, JDATE, JTIME, ADJ_LGRID)
            CALL WR_ADJ_FILE(LGRID_EM, JDATE, JTIME, ADJ_LGRID_EM)
            CALL WR_ADJ_FILE(LGRID_EM_SF, JDATE, JTIME, ADJ_LGRID_EM_SF)
            CALL WR_ADJ_FILE(LGRID_EM_NRM, JDATE, JTIME, 
     &                       ADJ_LGRID_EM_NRM)
         END IF
     
      END DO
      IF ( TRIM(LGRID_OUTPUT_FREQ) == 'ONCE' ) THEN
           CALL WR_ADJ_FILE(LGRID, JDATE, JTIME, ADJ_LGRID)
           CALL WR_ADJ_FILE(LGRID_EM, JDATE, JTIME, ADJ_LGRID_EM)
           CALL WR_ADJ_FILE(LGRID_EM_SF, JDATE, JTIME, ADJ_LGRID_EM_SF)
           CALL WR_ADJ_FILE(LGRID_EM_NRM, JDATE, JTIME,
     &                      ADJ_LGRID_EM_NRM)
      END IF

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

      END PROGRAM DRIVER_BWD
