
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE VADV_FLUX_DEFN

C Function:
C   Measure what vertical advection moves across each layer interface, and in
C   particular across the model top.

C   CMAQ gives methane no top boundary condition. ZADV diagnoses VEL(NLAYS+1)
C   from the column air-mass budget residual rather than setting it to zero,
C   and VPPM returns the top layer's own concentration for air re-entering
C   through the lid, so the interface works as a valve: methane leaves at one
C   time's top-layer value and returns at a later one's. Net flux alone hides
C   that, because the valve can run hard with the two directions nearly
C   cancelling, so the up and down transports are kept apart here.

C   Two outputs, with different jobs:

C     - a log line per sync step, always on, carrying the domain-integrated
C       flux at every interface. Summing these closes the vertical mass budget
C       layer by layer.
C     - an optional gridded file holding the model-top interface alone, at the
C       output time step, which says where the lid flux happens. Written only
C       when CTM_VADV_TOPFLX names a file.

C   Both are instrumentation. Nothing here is read back into the solution.

C   Includes are kept inside the contained procedures: ZADV includes the same
C   species tables and I/O API headers itself, and use association would
C   collide with them.

C Revision History:
C   Sep 2026 by Open Methane: created for openmethane/openmethane#256
C-----------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER( 16 ), PARAMETER :: CTM_VADV_TOPFLX = 'CTM_VADV_TOPFLX'

C The advected variables, in ZADV's own order: the trace species followed by
C air in the last slot.

      INTEGER, SAVE :: N_FLX_VAR = 0

C Accumulators, in VPPM's units of advected quantity times sigma thickness.
C Scaling to mass is deferred to report time so the inner loops stay cheap.

C Every interface, summed over this processor's columns and over the sync
C step's sub-steps. Double precision because this is a domain sum of many
C small same-signed terms.

      REAL( 8 ), ALLOCATABLE, SAVE :: FLX_SYNC_UP( :,: ) ! ( NLAYS+1, N_FLX_VAR )
      REAL( 8 ), ALLOCATABLE, SAVE :: FLX_SYNC_DN( :,: )

C The model top interface alone, per column, summed over the output time step.

      REAL, ALLOCATABLE, SAVE :: FLX_TOP_UP( :,:,: ) ! ( MY_NCOLS, MY_NROWS, N_FLX_VAR )
      REAL, ALLOCATABLE, SAVE :: FLX_TOP_DN( :,:,: )

      LOGICAL, SAVE :: FLX_FILE = .FALSE. ! write the gridded model-top file

C kg per unit of accumulator. Zero marks a variable whose coupled units this
C module cannot convert, which is excluded from both outputs.

      REAL, ALLOCATABLE, SAVE :: FLX_CONV( : )
      CHARACTER( 16 ), ALLOCATABLE, SAVE :: FLX_NAME( : )

      INTEGER, SAVE :: OUT_SDATE, OUT_STIME ! start of the open output interval
      INTEGER, SAVE :: FLX_LOGDEV

      CONTAINS

C-----------------------------------------------------------------------
         SUBROUTINE VADV_FLUX_INIT ( JDATE, JTIME )

C Allocate the accumulators and work out each variable's mass conversion.

         USE GRID_CONF          ! horizontal & vertical domain specifications

         IMPLICIT NONE

         INCLUDE SUBST_GC_SPC   ! gas chemistry species table
         INCLUDE SUBST_AE_SPC   ! aerosol species table
         INCLUDE SUBST_NR_SPC   ! non-reactive species table
         INCLUDE SUBST_TR_SPC   ! tracer species table
         INCLUDE SUBST_GC_ADV   ! gas chem advection species and map table
         INCLUDE SUBST_AE_ADV   ! aerosol advection species and map table
         INCLUDE SUBST_NR_ADV   ! non-react advection species and map table
         INCLUDE SUBST_TR_ADV   ! tracer advection species and map table
         INCLUDE SUBST_CONST    ! physical constants
         INCLUDE SUBST_IOPARMS  ! I/O parameters definitions
         INCLUDE SUBST_IODECL   ! I/O definitions and declarations

         INTEGER, INTENT( IN ) :: JDATE ! current model date, coded YYYYDDD
         INTEGER, INTENT( IN ) :: JTIME ! current model time, coded HHMMSS

         CHARACTER( 16 ) :: PNAME = 'VADV_FLUX_INIT'
         CHARACTER( 96 ) :: XMSG = ' '
         CHARACTER( 80 ) :: VARDESC = ' '
         CHARACTER( 512 ) :: FLXFILE = ' ' ! only the presence of the value matters

         INTEGER ALLOCSTAT, STATUS
         INTEGER SPC, VAR
         REAL PPM_FAC           ! species conversion, less the molecular weight
         REAL CELLAREA          ! nominal grid cell area, m**2

C-----------------------------------------------------------------------

         FLX_LOGDEV = INIT3()

         N_FLX_VAR = N_GC_ADV + N_AE_ADV + N_NR_ADV + N_TR_ADV + 1

         ALLOCATE ( FLX_SYNC_UP( NLAYS+1, N_FLX_VAR ),
     &        FLX_SYNC_DN( NLAYS+1, N_FLX_VAR ),
     &        FLX_CONV( N_FLX_VAR ),
     &        FLX_NAME( N_FLX_VAR ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating vertical advection flux diagnostic'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         FLX_SYNC_UP = 0.0D0
         FLX_SYNC_DN = 0.0D0

C The advected quantity is a coupled concentration, so accumulator times cell
C area is already extensive. Gas, non-reactive and tracer species are coupled
C as ppmV x RHOJ, and RHOJ already carries the map scale factor
C (SqRDMT = Jacobian / m**2), which is why the nominal cell area is the right
C multiplier and no further map correction belongs here.

         CELLAREA = REAL( XCELL_GD * YCELL_GD )
         PPM_FAC = 1.0E-6 * CELLAREA / MWAIR

         SPC = 0

         DO VAR = 1, N_GC_ADV
            SPC = SPC + 1
            FLX_NAME( SPC ) = GC_ADV( VAR )
            FLX_CONV( SPC ) = PPM_FAC * GC_MOLWT( GC_ADV_MAP( VAR ) )
         END DO

C Aerosols are coupled through JACOBM in mass or number per volume, not through
C RHOJ, so the ppmV conversion above does not apply to them. The CH4-only
C mechanism this fork builds has none; anything else is reported as excluded
C rather than converted wrongly.

         DO VAR = 1, N_AE_ADV
            SPC = SPC + 1
            FLX_NAME( SPC ) = AE_ADV( VAR )
            FLX_CONV( SPC ) = 0.0
         END DO

         DO VAR = 1, N_NR_ADV
            SPC = SPC + 1
            FLX_NAME( SPC ) = NR_ADV( VAR )
            FLX_CONV( SPC ) = PPM_FAC * NR_MOLWT( NR_ADV_MAP( VAR ) )
         END DO

         DO VAR = 1, N_TR_ADV
            SPC = SPC + 1
            FLX_NAME( SPC ) = TR_ADV( VAR )
            FLX_CONV( SPC ) = PPM_FAC * TR_MOLWT( TR_ADV_MAP( VAR ) )
         END DO

C Air rides in the last slot, advected as RHOJ itself, so its accumulator is
C already a mass per unit area.

         FLX_NAME( N_FLX_VAR ) = 'AIR'
         FLX_CONV( N_FLX_VAR ) = CELLAREA

         IF ( N_AE_ADV .GT. 0 ) THEN
            XMSG = 'Aerosol species excluded from the vertical flux diagnostic'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
         END IF

         VARDESC = 'Gridded model-top advective flux file'
         CALL ENVSTR ( CTM_VADV_TOPFLX, VARDESC, ' ', FLXFILE, STATUS )
         FLX_FILE = ( STATUS .EQ. 0 )

         IF ( FLX_FILE ) THEN
            ALLOCATE ( FLX_TOP_UP( MY_NCOLS, MY_NROWS, N_FLX_VAR ),
     &           FLX_TOP_DN( MY_NCOLS, MY_NROWS, N_FLX_VAR ),
     &           STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating model-top flux grid'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            FLX_TOP_UP = 0.0
            FLX_TOP_DN = 0.0
            WRITE( FLX_LOGDEV,* ) ' '
            WRITE( FLX_LOGDEV,* ) '    Model-top advective flux written to '
     &           // CTM_VADV_TOPFLX
         ELSE
            WRITE( FLX_LOGDEV,* ) ' '
            WRITE( FLX_LOGDEV,* ) '    ' // CTM_VADV_TOPFLX
     &           // ' not set - gridded model-top flux file not written'
         END IF

         OUT_SDATE = JDATE
         OUT_STIME = JTIME

         RETURN
         END SUBROUTINE VADV_FLUX_INIT

C-----------------------------------------------------------------------
         CHARACTER( 16 ) FUNCTION FLX_VNAME ( SPC, SUFFIX )

C Name a variable in the gridded file. Species names are blank padded to 16,
C so they have to be trimmed before the direction is appended or the suffix
C lands past the end of the name.

         IMPLICIT NONE

         INTEGER, INTENT( IN ) :: SPC          ! advected variable
         CHARACTER( * ), INTENT( IN ) :: SUFFIX ! '_UP' or '_DN'

         INTEGER N

C-----------------------------------------------------------------------

         N = MIN( LEN_TRIM( FLX_NAME( SPC ) ), 16 - LEN( SUFFIX ) )
         FLX_VNAME = FLX_NAME( SPC )( 1:N ) // SUFFIX

         RETURN
         END FUNCTION FLX_VNAME

C-----------------------------------------------------------------------
         SUBROUTINE VADV_FLUX_ADD ( SPC, COL, ROW, FUP, FDN )

C Take up one VPPM call's interface fluxes for one column and one advected
C variable. Called from inside ZADV's column loop, once per advection sub-step.

         USE GRID_CONF          ! horizontal & vertical domain specifications

         IMPLICIT NONE

         INTEGER, INTENT( IN ) :: SPC   ! advected variable, in ZADV's ordering
         INTEGER, INTENT( IN ) :: COL   ! column, local to this processor
         INTEGER, INTENT( IN ) :: ROW   ! row, local to this processor
         REAL, INTENT( IN ) :: FUP( NLAYS+1 ) ! upward transport at each interface
         REAL, INTENT( IN ) :: FDN( NLAYS+1 ) ! downward transport at each interface

         INTEGER LVL

C-----------------------------------------------------------------------

         DO LVL = 1, NLAYS+1
            FLX_SYNC_UP( LVL,SPC ) = FLX_SYNC_UP( LVL,SPC ) + FUP( LVL )
            FLX_SYNC_DN( LVL,SPC ) = FLX_SYNC_DN( LVL,SPC ) + FDN( LVL )
         END DO

         IF ( FLX_FILE ) THEN
            FLX_TOP_UP( COL,ROW,SPC ) = FLX_TOP_UP( COL,ROW,SPC )
     &           + FUP( NLAYS+1 )
            FLX_TOP_DN( COL,ROW,SPC ) = FLX_TOP_DN( COL,ROW,SPC )
     &           + FDN( NLAYS+1 )
         END IF

         RETURN
         END SUBROUTINE VADV_FLUX_ADD

C-----------------------------------------------------------------------
         SUBROUTINE VADV_FLUX_REPORT ( JDATE, JTIME, TSTEP )

C Log this sync step's domain-integrated interface fluxes, and write the
C gridded model-top file when an output interval closes. Called once per sync
C step, after ZADV has advected every column.

         USE GRID_CONF          ! horizontal & vertical domain specifications

         USE SUBST_MODULES      ! stenex

         IMPLICIT NONE

         INTEGER, INTENT( IN ) :: JDATE      ! date at the start of the sync step
         INTEGER, INTENT( IN ) :: JTIME      ! time at the start of the sync step
         INTEGER, INTENT( IN ) :: TSTEP( 2 ) ! (1) output step, (2) sync step, HHMMSS

         INTEGER, EXTERNAL :: TIME2SEC
         INTEGER, EXTERNAL :: SECSDIFF

         INTEGER EDATE, ETIME   ! date/time at the end of the sync step
         INTEGER LVL, SPC
         REAL    DTSEC          ! sync step, seconds
         REAL    UP, DN
         REAL    PROFILE( NLAYS+1 )

C-----------------------------------------------------------------------

         DTSEC = FLOAT( TIME2SEC( TSTEP( 2 ) ) )

C One scalar reduction per interface per direction. The stenex no-op linked
C into the serial build has no array form, so this cannot be a single call.

         DO SPC = 1, N_FLX_VAR
            IF ( FLX_CONV( SPC ) .LE. 0.0 ) CYCLE
            DO LVL = 1, NLAYS+1
               FLX_SYNC_UP( LVL,SPC ) = SUBST_GLOBAL_SUM(
     &              REAL( FLX_SYNC_UP( LVL,SPC ) ) )
               FLX_SYNC_DN( LVL,SPC ) = SUBST_GLOBAL_SUM(
     &              REAL( FLX_SYNC_DN( LVL,SPC ) ) )
            END DO
         END DO

         IF ( MYPE .EQ. 0 ) THEN
            DO SPC = 1, N_FLX_VAR
               IF ( FLX_CONV( SPC ) .LE. 0.0 ) CYCLE

               UP = FLX_CONV( SPC ) * REAL( FLX_SYNC_UP( NLAYS+1,SPC ) ) / DTSEC
               DN = FLX_CONV( SPC ) * REAL( FLX_SYNC_DN( NLAYS+1,SPC ) ) / DTSEC
               WRITE( FLX_LOGDEV,2001 ) JDATE, JTIME, FLX_NAME( SPC )( 1:8 ),
     &              UP, DN, UP - DN

               DO LVL = 1, NLAYS+1
                  PROFILE( LVL ) =
     &                 FLX_CONV( SPC ) * REAL( FLX_SYNC_UP( LVL,SPC ) ) / DTSEC
               END DO
               WRITE( FLX_LOGDEV,2003 ) JDATE, JTIME, FLX_NAME( SPC )( 1:8 ),
     &              'UP', ( PROFILE( LVL ), LVL = 1, NLAYS+1 )

               DO LVL = 1, NLAYS+1
                  PROFILE( LVL ) =
     &                 FLX_CONV( SPC ) * REAL( FLX_SYNC_DN( LVL,SPC ) ) / DTSEC
               END DO
               WRITE( FLX_LOGDEV,2003 ) JDATE, JTIME, FLX_NAME( SPC )( 1:8 ),
     &              'DN', ( PROFILE( LVL ), LVL = 1, NLAYS+1 )
            END DO
         END IF

2001     FORMAT( ' VADVTOP ', I7.7, 1X, I6.6, 1X, A8,
     &        ' up', 1PE12.4, ' dn', 1PE12.4, ' net', 1PE12.4, ' kg/s' )
2003     FORMAT( ' VADVFLX ', I7.7, 1X, I6.6, 1X, A8, 1X, A2, ' kg/s',
     &        999( 1PE12.4 ) )

         FLX_SYNC_UP = 0.0D0
         FLX_SYNC_DN = 0.0D0

         IF ( .NOT. FLX_FILE ) RETURN

         EDATE = JDATE
         ETIME = JTIME
         CALL NEXTIME ( EDATE, ETIME, TSTEP( 2 ) )

         IF ( SECSDIFF( OUT_SDATE, OUT_STIME, EDATE, ETIME )
     &        .GE. TIME2SEC( TSTEP( 1 ) ) ) THEN
            CALL WR_VADV_TOPFLX ( EDATE, ETIME, TSTEP )
            OUT_SDATE = EDATE
            OUT_STIME = ETIME
         END IF

         RETURN
         END SUBROUTINE VADV_FLUX_REPORT

C-----------------------------------------------------------------------
         SUBROUTINE WR_VADV_TOPFLX ( EDATE, ETIME, TSTEP )

C Write the model-top interface flux for the output interval that has just
C closed, stamped with the interval's start as CMAQ's averaged files are.

         USE GRID_CONF          ! horizontal & vertical domain specifications

         IMPLICIT NONE

         INCLUDE SUBST_IOPARMS  ! I/O parameters definitions
         INCLUDE SUBST_IOFDESC  ! file header data structure
#include         SUBST_IODECL   ! I/O definitions and declarations

         INTEGER, INTENT( IN ) :: EDATE      ! date at the end of the interval
         INTEGER, INTENT( IN ) :: ETIME      ! time at the end of the interval
         INTEGER, INTENT( IN ) :: TSTEP( 2 ) ! (1) output step, (2) sync step, HHMMSS

         INTEGER, EXTERNAL :: SECSDIFF

         CHARACTER( 16 ) :: PNAME = 'WR_VADV_TOPFLX'
         CHARACTER( 96 ) :: XMSG = ' '
         CHARACTER( 16 ) :: VNAME

         LOGICAL, SAVE :: FIRSTIME = .TRUE.

         INTEGER COL, ROW, SPC
         REAL    DTSEC          ! length of the closed interval, seconds
         REAL    GRID( MY_NCOLS, MY_NROWS )

C-----------------------------------------------------------------------

         DTSEC = FLOAT( SECSDIFF( OUT_SDATE, OUT_STIME, EDATE, ETIME ) )

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            CALL OPEN_VADV_TOPFLX ( TSTEP )
         END IF

         DO SPC = 1, N_FLX_VAR
            IF ( FLX_CONV( SPC ) .LE. 0.0 ) CYCLE

            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  GRID( COL,ROW ) =
     &                 FLX_CONV( SPC ) * FLX_TOP_UP( COL,ROW,SPC ) / DTSEC
               END DO
            END DO
            VNAME = FLX_VNAME( SPC, '_UP' )
            IF ( .NOT. WRITE3( CTM_VADV_TOPFLX, VNAME,
     &           OUT_SDATE, OUT_STIME, GRID ) ) THEN
               XMSG = 'Could not write ' // VNAME // ' to ' // CTM_VADV_TOPFLX
               CALL M3EXIT ( PNAME, OUT_SDATE, OUT_STIME, XMSG, XSTAT1 )
            END IF

            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  GRID( COL,ROW ) =
     &                 FLX_CONV( SPC ) * FLX_TOP_DN( COL,ROW,SPC ) / DTSEC
               END DO
            END DO
            VNAME = FLX_VNAME( SPC, '_DN' )
            IF ( .NOT. WRITE3( CTM_VADV_TOPFLX, VNAME,
     &           OUT_SDATE, OUT_STIME, GRID ) ) THEN
               XMSG = 'Could not write ' // VNAME // ' to ' // CTM_VADV_TOPFLX
               CALL M3EXIT ( PNAME, OUT_SDATE, OUT_STIME, XMSG, XSTAT1 )
            END IF
         END DO

         FLX_TOP_UP = 0.0
         FLX_TOP_DN = 0.0

         RETURN
         END SUBROUTINE WR_VADV_TOPFLX

C-----------------------------------------------------------------------
         SUBROUTINE OPEN_VADV_TOPFLX ( TSTEP )

C Create the model-top flux file. The header is built from the model's own grid
C definition rather than copied from another open file, because MET_CRO_3D may
C cover a larger domain than the one being run.

         USE GRID_CONF          ! horizontal & vertical domain specifications

         USE SUBST_MODULES      ! stenex

         IMPLICIT NONE

         INCLUDE SUBST_IOPARMS  ! I/O parameters definitions
         INCLUDE SUBST_IOFDESC  ! file header data structure
#include         SUBST_IODECL   ! I/O definitions and declarations

         INTEGER, INTENT( IN ) :: TSTEP( 2 ) ! (1) output step, (2) sync step, HHMMSS

         CHARACTER( 16 ) :: PNAME = 'OPEN_VADV_TOPFLX'
         CHARACTER( 96 ) :: XMSG = ' '

         INTEGER SPC, VAR, K
         LOGICAL OK

C-----------------------------------------------------------------------

         CALL SUBST_BARRIER

         OK = OPEN3( CTM_VADV_TOPFLX, FSRDWR3, PNAME )
         CALL SUBST_GLOBAL_LOGICAL( OK, 'AND' )

         IF ( .NOT. OK .AND. MYPE .EQ. 0 ) THEN

            FTYPE3D = GRDDED3
            SDATE3D = OUT_SDATE
            STIME3D = OUT_STIME
            TSTEP3D = TSTEP( 1 )

            NCOLS3D = GL_NCOLS
            NROWS3D = GL_NROWS
            NLAYS3D = 1
            NTHIK3D = NTHIK
            GDTYP3D = GDTYP_GD
            P_ALP3D = P_ALP_GD
            P_BET3D = P_BET_GD
            P_GAM3D = P_GAM_GD
            XCENT3D = XCENT_GD
            YCENT3D = YCENT_GD
            XORIG3D = XORIG_GD
            YORIG3D = YORIG_GD
            XCELL3D = XCELL_GD
            YCELL3D = YCELL_GD
            VGTYP3D = VGTYP_GD
            VGTOP3D = VGTOP_GD
            VGLVS3D( 1 ) = VGLVS_GD( NLAYS )
            VGLVS3D( 2 ) = VGLVS_GD( NLAYS+1 )
            GDNAM3D = GRID_NAME

            FDESC3D( 1 ) = 'Advective flux across the model top interface'
            FDESC3D( 2 ) = 'Mean over the output time step'
            FDESC3D( 3 ) = 'Timestamp represents beginning computed date/time'
            FDESC3D( 4 ) = 'Upward and downward transports are reported'
            FDESC3D( 5 ) = 'separately; net upward flux is _UP minus _DN'
            DO K = 6, MXDESC3
               FDESC3D( K ) = ' '
            END DO

            VAR = 0
            DO SPC = 1, N_FLX_VAR
               IF ( FLX_CONV( SPC ) .LE. 0.0 ) CYCLE

               VAR = VAR + 1
               VTYPE3D( VAR ) = M3REAL
               VNAME3D( VAR ) = FLX_VNAME( SPC, '_UP' )
               UNITS3D( VAR ) = 'kg/s'
               VDESC3D( VAR ) = TRIM( FLX_NAME( SPC ) )
     &              // ' carried up through the model top'

               VAR = VAR + 1
               VTYPE3D( VAR ) = M3REAL
               VNAME3D( VAR ) = FLX_VNAME( SPC, '_DN' )
               UNITS3D( VAR ) = 'kg/s'
               VDESC3D( VAR ) = TRIM( FLX_NAME( SPC ) )
     &              // ' carried down through the model top'
            END DO
            NVARS3D = VAR

            IF ( .NOT. OPEN3( CTM_VADV_TOPFLX, FSNEW3, PNAME ) ) THEN
               XMSG = 'Could not open ' // CTM_VADV_TOPFLX // ' file'
               CALL M3EXIT ( PNAME, OUT_SDATE, OUT_STIME, XMSG, XSTAT1 )
            END IF

         END IF

         CALL SUBST_BARRIER

         RETURN
         END SUBROUTINE OPEN_VADV_TOPFLX

      END MODULE VADV_FLUX_DEFN
