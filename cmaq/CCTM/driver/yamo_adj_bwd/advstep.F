
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_bwd/advstep.F,v 1.3 2013/05/15 20:14:11 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ADVSTEP(JDATE, JTIME, TSTEP, ASTEP)

C-----------------------------------------------------------------------
C Function:
C    Use Courant-number conditions to determine the number of advection
C    sub-steps of the synchronization step needed in each layer.

C Preconditions:
C    The synchronization step, TSTEP(2), has already been determined.

C Revision history:
C   Aug 2010 by Peter Percell at UH-IMAQS: rewrote for adjoint/4dvar
C   May 2013 by Peter Percell at UH-IMAQS:
C   - changed default HDIV_LIM to 0.9 (as in CMAQ_5.0.1)
C   - changed the calling argument NADVSS back to the traditional ASTEP
C     which is once again required to be an integer divisor of the sync step
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
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                !   TSTEP(1) = local output step
                                !   TSTEP(2) = sciproc sync. step (chem)
      INTEGER, INTENT(OUT) :: ASTEP(:) ! advection time steps

C External Functions (not already declared by IODECL3.EXT):
      INTEGER, EXTERNAL :: SECSDIFF, TIME2SEC, SEC2TIME
      LOGICAL, EXTERNAL :: CURRSTEP
      REAL, EXTERNAL    :: ENVREAL

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV   ! unit number for log device

      CHARACTER(16) :: PNAME = 'ADVSTEP'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(16) :: ADV_HDIV_LIM = 'ADV_HDIV_LIM' ! env. var. max horiz. div.

      INTEGER, SAVE :: WSTEP = 0 ! wind file interpolation time step
      INTEGER, SAVE :: WSDATE, WSTIME ! wind file start date and time
      INTEGER, SAVE :: FDATE, FTIME ! wind file current date and time
      INTEGER, SAVE :: FSTEP    ! wind file time step

      REAL, SAVE :: IDX1        ! 1/dx1
      REAL, SAVE :: IDX2        ! 1/dx2

      REAL, SAVE :: HDIV_LIM    ! cutoff for max horizontal divergence step adj

      INTEGER :: STATUS         ! ENVREAL status
      CHARACTER(80) :: VARDESC  ! environment variable description

      REAL :: UWIND(1:NCOLS+1, 1:NROWS+1, 1:NLAYS)
      REAL :: VWIND(1:NCOLS+1, 1:NROWS+1, 1:NLAYS)
      REAL :: HDIV(1:NCOLS, 1:NROWS, 1:NLAYS)

      INTEGER :: C, R, L, K     ! loop counters
      INTEGER :: EDATE, ETIME   ! end date & time
      INTEGER :: SYNC           ! synchronization step (sec)

      REAL :: MXADV             ! max acceptable advection step (sec)
      INTEGER :: DT_INT         ! advection step (sec)

      REAL :: WIND_DSI(1:NLAYS) ! max { component velocity / dXi }
      REAL :: MX_HDIV(1:NLAYS)  ! max(HDIV) for whole output step

      LOGICAL :: WIND_DSI_DET   ! true if advection step is determined by WIND_DSI

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

C Get the wind file start date and time for this output step; do not interpolate
      IF ( .NOT. CURRSTEP(JDATE, JTIME, WSDATE, WSTIME, FSTEP,
     &     FDATE, FTIME) ) THEN
         XMSG = 'Cannot get starting date and time'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
      END IF

C Establish end time for Courant number calculation
      EDATE = JDATE
      ETIME = JTIME
      CALL NEXTIME(EDATE, ETIME, TSTEP(1))

      WIND_DSI(:) = 0.0
      MX_HDIV(:) = 0.0
      DO
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
            IF ( SECSDIFF(EDATE, ETIME, FDATE, FTIME) > 0 )
     &           CALL NEXTIME(FDATE, FTIME, -FSTEP) 
            EXIT
         ELSE
            CALL NEXTIME(FDATE, FTIME, FSTEP) ! next wind file time
         END IF
      END DO

      DO L = 1, NLAYS
         WIND_DSI(L) = SUBST_GLOBAL_MAX(WIND_DSI(L))
         MX_HDIV(L) = SUBST_GLOBAL_MAX(MX_HDIV(L))
      END DO

C For each layer, determine the largest integer advection time step (sec),
C that satisfies the Courant condition and is an integer divisor of the
C synchronization time step (sec)
      WRITE(LOGDEV, 92025)
      SYNC = TIME2SEC(TSTEP(2))
      DO L = NLAYS, 1, -1
         IF ( CC / WIND_DSI(L) <= HDIV_LIM / MX_HDIV(L) ) THEN
            MXADV = CC / WIND_DSI(L)
            WIND_DSI_DET = .TRUE.

            IF ( MXADV < 1 ) THEN
               WRITE(LOGDEV, 94013) TSTEP(1), TSTEP(2), L, WIND_DSI(L)
               XMSG = ' '
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
            END IF
         ELSE
            MXADV = HDIV_LIM / MX_HDIV(L)
            WIND_DSI_DET = .FALSE.

            IF ( MXADV < 1 ) THEN
               WRITE(LOGDEV, 94014) TSTEP(1), TSTEP(2), L, MX_HDIV(L)
               XMSG = ' '
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
            END IF
         END IF

C Note that this DO WHILE loop must end by the time DT_INT reaches 1 since
C MXADV & SYNC >= 1
         DT_INT = MIN(FLOOR(MXADV), SYNC)
         DO WHILE ( MOD(SYNC, DT_INT) /= 0 )
            DT_INT = DT_INT - 1
         END DO
         ASTEP(L) = SEC2TIME(DT_INT)

         IF ( WIND_DSI_DET ) THEN
            WRITE(LOGDEV, 92027) L, ASTEP(L), SYNC/DT_INT
         ELSE
            WRITE(LOGDEV, 92029) L, ASTEP(L), SYNC/DT_INT
         END IF
      END DO

      RETURN

C------------------  FORMAT  STATEMENTS   ------------------------------

92013 FORMAT(/ 5X, 'From ADVSTEP - date/time: ', I8, '/', I6.6)

92025 FORMAT(/ 5X, 'Layer', 3X, 'Advection', 3X, 'per Sync'
     &       /11X, 'Step (HHMMSS)', 2X, 'Step')

92027 FORMAT(5X, I4, 6X, I6.6, 6X, I2)
92029 FORMAT(5X, I4, 6X, I6.6, 6X, I2, 2X, '<- determined by HDIV_LIM')

94013 FORMAT(/ 5X, 'Could not determine satisfactory advection time step:',
     &       / 8X, 'output time step = ', I6.6, ' (HHMMSS)',
     &       / 8X, 'synchronization time step = ', I6.6, ' (HHMMSS)',
     &       / 8X, 'layer = ', I3
     &       / 8X, '(Max vel)/(dX) = ', 1PE10.3)
94014 FORMAT(/ 5X, 'Could not determine satisfactory advection time step:',
     &       / 8X, 'output time step = ', I6.6, ' (HHMMSS)',
     &       / 8X, 'synchronization time step = ', I6.6, ' (HHMMSS)',
     &       / 8X, 'layer = ', I3
     &       / 8X, 'Max horizontal divergence = ', 1PE10.3)

      END
