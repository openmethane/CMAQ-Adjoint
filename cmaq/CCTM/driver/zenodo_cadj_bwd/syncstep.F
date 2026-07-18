
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_bwd/syncstep.F,v 1.2 2013/05/15 20:14:11 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SYNCSTEP(JDATE, JTIME, TSTEP, NSTEPS, NREPS)

C-----------------------------------------------------------------------
C Function:
C    Use Courant-number conditions to determine the (time-splitting) model
C    base synchronization time step TSTEP(2) and number NREPS of sync time
C    steps per output time step TSTEP(1).

C Revision history:
C   Aug 2010 by Peter Percell at UH-IMAQS:
C   - adapted from advstep.F for adjoint/4dvar
C   May 2013 by Peter Percell at UH-IMAQS:
C   - changed default HDIV_LIM to 0.9 (as in CMAQ_5.0.1)
C   - support multiple calls so that the sync step can change with the
C     output step if called with NSTEPS=1
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE SUBST_MODULES         ! stenex

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Parameters:
      REAL, PARAMETER :: CC = 0.75 ! maximum Courant number allowed
c$$$      REAL, PARAMETER :: CC = 0.90 ! maximum Courant number allowed

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! current model simulation date (YYYYDDD)
      INTEGER, INTENT(IN) :: JTIME ! current model simulation time (HHMMSS)
      INTEGER, INTENT(INOUT) :: TSTEP(2) ! time step vector (HHMMSS)
                                !   TSTEP(1) = local output step
                                !   TSTEP(2) = sciproc sync. step (chem)
      INTEGER, INTENT(IN)  :: NSTEPS ! run duration: number of output time steps
      INTEGER, INTENT(OUT) :: NREPS ! number of sync steps per output step

C External Functions (not already declared by IODECL3.EXT):
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC, SECSDIFF, ENVINT
      LOGICAL, EXTERNAL :: CURRSTEP
      REAL, EXTERNAL    :: ENVREAL

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV   ! unit number for log device

      CHARACTER(16) :: PNAME = 'SYNCSTEP'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(16) :: CTM_MAXSYNC = 'CTM_MAXSYNC' ! env. var. max sync step
      CHARACTER(16) :: CTM_MINSYNC = 'CTM_MINSYNC' ! env. var. min sync step
      CHARACTER(16) :: SIGMA_SYNC_TOP = 'SIGMA_SYNC_TOP' ! env. var. adv layer
      CHARACTER(16) :: ADV_HDIV_LIM = 'ADV_HDIV_LIM' ! env. var. max horiz. div.

      INTEGER, SAVE :: MAXSYNC  ! force max TSTEP(2) (sec)
      INTEGER, SAVE :: MINSYNC  ! force min TSTEP(2) (sec)
      REAL, SAVE    :: SIGST    ! sigma_sync_top value
      INTEGER, SAVE :: ADVLAYR  ! top layer used to determine sync step

      INTEGER, SAVE :: OSTEP    ! output time step (sec)
      INTEGER, SAVE :: WSTEP = 0 ! wind file interpolation time step
      INTEGER, SAVE :: WSDATE, WSTIME ! wind file start date and time
      INTEGER, SAVE :: WEDATE, WETIME ! wind file end date and time
      INTEGER, SAVE :: FDATE, FTIME ! wind file current date and time
      INTEGER, SAVE :: FSTEP    ! wind file time step

      REAL, SAVE :: IDX1        ! 1/dx1
      REAL, SAVE :: IDX2        ! 1/dx2

      REAL, SAVE :: HDIV_LIM    ! cutoff for max horizontal divergence step adj

      INTEGER :: STATUS         ! ENVINT or ENVREAL status
      CHARACTER(80) :: VARDESC  ! environment variable description

      INTEGER :: MINSECS        ! effective MINSYNC (sec); must divide TSTEP(1)

      REAL :: UWIND(1:NCOLS+1, 1:NROWS+1, 1:NLAYS)
      REAL :: VWIND(1:NCOLS+1, 1:NROWS+1, 1:NLAYS)
      REAL :: HDIV(1:NCOLS, 1:NROWS, 1:NLAYS)

      INTEGER :: C, R, L, K     ! loop counters
      INTEGER :: EDATE, ETIME   ! end date & time
      INTEGER :: SYNC           ! synchronization step (sec)
      LOGICAL :: ADJFLG

      REAL :: MXADV             ! max acceptable advection step (sec)

      REAL :: WIND_DSI(1:NLAYS) ! max { component velocity / dXi }
      REAL :: MX_HDIV(1:NLAYS)  ! max(HDIV) for whole output step

      INTEGER :: MY_TEMP
      INTEGER, SAVE :: STARTCOL, ENDCOL
      INTEGER, SAVE :: STARTROW, ENDROW

      INTERFACE
        SUBROUTINE HVELOC(FDATE, FTIME, FSTEP, UWIND, VWIND, HDIV)
            IMPLICIT NONE
            INTEGER, INTENT(IN)  :: FDATE, FTIME, FSTEP
            REAL,    INTENT(OUT) :: UWIND(:, :, :), VWIND(:, :, :)
            REAL,    INTENT(OUT) :: HDIV(:, :, :)
         END SUBROUTINE HVELOC
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

         MAXSYNC = 720          ! default
         VARDESC = 'Maximum Synchronization Time Step (sec)'
         MAXSYNC = ENVINT(CTM_MAXSYNC, VARDESC, MAXSYNC, STATUS)
         IF ( STATUS /= 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
         IF ( STATUS == 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

         MINSYNC = 60           ! default
         VARDESC = 'Minimum Synchronization Time Step (sec)'
         MINSYNC = ENVINT(CTM_MINSYNC, VARDESC, MINSYNC, STATUS)
         IF ( STATUS /= 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
         IF ( STATUS == 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

C Make sure MINSYNC <= MAXSYNC
         IF ( MAXSYNC < MINSYNC ) THEN
            WRITE(LOGDEV, 92010) MAXSYNC, MINSYNC
            XMSG = ' '
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
         END IF

C Output time step (sec)
         OSTEP = TIME2SEC(TSTEP(1))

C Make sure MINSYNC <= output time step
         IF ( OSTEP < MINSYNC ) THEN
            WRITE(LOGDEV, 92009) OSTEP, MINSYNC
            XMSG = ' '
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
         END IF

C Make sure that the effective MINSYNC divides TSTEP(1): determine MINSECS,
C the smallest time step (sec) between MINSYNC and MAXSYNC that divides the
C output step
         IF ( MINSYNC <= 1 ) THEN
            MINSECS = 1
         ELSE
            DO MINSECS = MINSYNC, MAXSYNC
               IF ( MOD(OSTEP, MINSECS) == 0 ) THEN
                  IF ( MINSECS /= MINSYNC )
     &                 WRITE(LOGDEV, 92011) JDATE, JTIME, MINSYNC, MINSECS
                  EXIT
               ELSE IF ( MINSECS < MAXSYNC ) THEN
                  CYCLE
               END IF

C A satisfactory MINSECS could not be found
               WRITE(XMSG, 92012) MINSYNC, TSTEP(1)
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
            END DO
         END IF

C Open wind field file and get header data
         IF ( .NOT. OPEN3(MET_DOT_3D, FSREAD3, PNAME) ) THEN
            XMSG = 'Could not open ' // MET_DOT_3D // ' file'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

         IF ( .NOT. DESC3(MET_DOT_3D) ) THEN
            XMSG = 'Could not get ' // MET_DOT_3D // ' file description'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

C Get the wind file start date and time and the wind file time step
         WSDATE = SDATE3D
         WSTIME = STIME3D
         FSTEP = TSTEP3D

         WEDATE = SDATE3D
         WETIME = STIME3D
         CALL NEXTIME(WEDATE, WETIME, SEC2TIME((MXREC3D - 1) * TIME2SEC(FSTEP)))

         IF ( VGTYP_GD /= VGSGPN3 .AND.
     &        VGTYP_GD /= VGSGPH3 .AND.
     &        VGTYP_GD /= VGWRFEM ) THEN
            XMSG = 'Wrong vertical coordinate type'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

         SIGST = 0.7            ! default
         VARDESC = 'Minimum layer limit for which adv = sync'
         SIGST = ENVREAL(SIGMA_SYNC_TOP, VARDESC, SIGST, STATUS)
         IF ( STATUS /= 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS == 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

         IF ( SIGST > VGLVS_GD(2) .OR.
     &        SIGST < VGLVS_GD(NLAYS + 1) ) THEN
            WRITE(XMSG, '(A, 1PE12.3)') 'SIGMA_SYNC_TOP incorrect', SIGST
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

         IF ( SIGST == VGLVS_GD(NLAYS + 1) ) THEN
            L = NLAYS
         ELSE
            DO K = 1, NLAYS
               IF ( SIGST >= VGLVS_GD(K+1) ) THEN
                  IF ( SIGST - VGLVS_GD(K+1) <=
     &                 VGLVS_GD(K) - SIGST )  THEN
                     L = K
                  ELSE
                     L = K - 1
                  END IF
                  EXIT
               END IF
            END DO
         END IF

         ADVLAYR = L
         WRITE(LOGDEV, 92005) ADVLAYR

C Check file data
         IF ( XCELL3D /= XCELL_GD .OR. YCELL3D /= YCELL_GD ) THEN
            XMSG = 'File grid sizes do not match CTM domain definition'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
         END IF

C Get cell size in meters
         IF ( GDTYP_GD == LATGRD3 ) THEN
            IDX2 = 1.0 / (DG2M * YCELL_GD)
            IDX1 = 1.0 / (DG2M * XCELL_GD
     &           * COS(PI180 * (YORIG_GD + YCELL_GD * FLOAT(NROWS / 2))))
         ELSE
            IDX1 = 1.0 / XCELL_GD
            IDX2 = 1.0 / YCELL_GD
         END IF

         CALL SUBST_LOOP_INDEX('C', 1, NCOLS, 1, MY_TEMP,
     &        STARTCOL, ENDCOL)

         CALL SUBST_LOOP_INDEX('R', 1, NROWS, 1, MY_TEMP,
     &        STARTROW, ENDROW)

         HDIV_LIM = 0.9
         VARDESC = 'Maximum horiz. div. limit for adv step adjustment'
         HDIV_LIM = ENVREAL(ADV_HDIV_LIM, VARDESC, HDIV_LIM, STATUS)
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV,'(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

      END IF                    ! FIRSTIME

      WRITE(LOGDEV, 92013) JDATE, JTIME

C Establish end time for Courant number calculation
      EDATE = JDATE
      ETIME = JTIME
      CALL NEXTIME(EDATE, ETIME, SEC2TIME(NSTEPS * OSTEP))

C Make sure that there is sufficient wind data for Courant number calculation
      IF ( SECSDIFF(JDATE, JTIME, WSDATE, WSTIME) > 0 .OR.
     &     SECSDIFF(WEDATE, WETIME, EDATE, ETIME) > 0 ) THEN
         XMSG = 'Insufficient wind data for Courant number calculation'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
      END IF

C Get the initial wind file date and time needed for determining the
C synchronization time step
      IF ( .NOT. CURRSTEP(JDATE, JTIME, WSDATE, WSTIME, FSTEP,
     &     FDATE, FTIME) ) THEN
         XMSG = 'Cannot get initial wind file date and time'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
      END IF

C Determine the synchronization time step, TSTEP(2), and NREPS, the
C corresponding number of sync steps per output step
      WIND_DSI(:) = 0.0
      MX_HDIV(:) = 0.0
      DO WHILE ( SECSDIFF(FDATE, FTIME, WEDATE, WETIME) >= 0 )
         CALL HVELOC(FDATE, FTIME, WSTEP, UWIND, VWIND, HDIV)

         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = STARTCOL, ENDCOL
                  WIND_DSI(L) = MAX(WIND_DSI(L), ABS(UWIND(C, R, L)) * IDX1)
               END DO
            END DO

            DO R = STARTROW, ENDROW
               DO C = 1, NCOLS
                  WIND_DSI(L) = MAX(WIND_DSI(L), ABS(VWIND(C, R, L)) * IDX2)
               END DO
            END DO

            DO R = 1, NROWS
               DO C = 1, NCOLS
                  MX_HDIV(L) = MAX(MX_HDIV(L), HDIV(C, R, L))
               END DO
            END DO
         END DO

         IF ( SECSDIFF(EDATE, ETIME, FDATE, FTIME) >= 0 ) THEN
            EXIT
         ELSE
            CALL NEXTIME(FDATE, FTIME, FSTEP) ! next wind file time
         END IF
      END DO

      DO L = 1, NLAYS
         WIND_DSI(L) = SUBST_GLOBAL_MAX(WIND_DSI(L))
         MX_HDIV(L) = SUBST_GLOBAL_MAX(MX_HDIV(L))
      END DO

C Find the largest time step that satisfies the Courant condition in all of
C the lower ADVLAYR layers
      MXADV = MIN(CC / MAXVAL(WIND_DSI(1:ADVLAYR)),
     &     HDIV_LIM / MAXVAL(MX_HDIV(1:ADVLAYR)))

C If this is smaller than MINSECS, then set the sync step to MINSECS
      IF ( MXADV < MINSECS ) THEN
         ADJFLG = .TRUE.
         NREPS = OSTEP / MINSECS
         SYNC = MINSECS

C Otherwise, set the sync step to the largest time step not greater than
C MAXSYNC that both satisfies the Courant condition and is an integer divisor
C of the output step (existence of such a sync step is guaranteed because, at
C this point, MINSECS satisfies both conditions)
      ELSE
         ADJFLG = .FALSE.
         MXADV = MIN(MXADV, REAL(MAXSYNC)) ! respect MAXSYNC
C NREPS must be at least CEILING(OSTEP/MXADV) for the Courant condition to be
C satisfied, and does not need to be larger than OSTEP/MINSECS because MINSECS
C is an acceptable sync step
         DO NREPS = CEILING(OSTEP / MXADV), OSTEP / MINSECS
            IF ( MOD(OSTEP, NREPS) == 0 ) EXIT ! make TSTEP(2) divide TSTEP(1)
         END DO
         SYNC = OSTEP / NREPS
      END IF
      TSTEP(2) = SEC2TIME(SYNC)

      IF ( ADJFLG ) THEN
         WRITE(LOGDEV, 92019) MINSECS, TSTEP(2), NREPS
      ELSE
         WRITE(LOGDEV, 92021) TSTEP(2), NREPS
      END IF

      RETURN

C------------------  FORMAT  STATEMENTS   ------------------------------

92005 FORMAT(5X, 'Top layer thru which sync step determined:', I3)

92009 FORMAT(/ 5X, 'Output time step:', I8,
     &     1X, 'less than minimum synchronization step:', I8)

92010 FORMAT(/ 5X, 'Maximum synchronization step:', I8,
     &     1X, 'less than minimum synchronization step:', I8)

92011 FORMAT(/ 5X, 'From SYNCSTEP - date/time: ', I8, '/', I6.6
     &     / 5X, 'Minimum Synchronization Step adjusted from:', I8,
     &     1X, 'to:', I8)

92012 FORMAT(/ 5X, 'Starting from:', I6, ',',
     &     1X, 'could not determine minimum step that divides TSTEP',
     &     1X, 'for model step:', I7.6, ' HHMMSS')

92013 FORMAT(/ 5X, 'From SYNCSTEP - date/time: ', I8, '/', I6.6)

92019 FORMAT(/ 5X, 'Synchronization step adjusted up to mimimum (SEC):', I7
     &     /46X, '(HHMMSS): ', I6.6
     &     / 5X, 'Number of Synchronization steps:', I3 /)

92021 FORMAT(/ 5X, 'Computed synchronization step (HHMMSS): ', I6.6
     &     / 5X, 'Number of Synchronization steps:', I3 /)

      END SUBROUTINE SYNCSTEP
