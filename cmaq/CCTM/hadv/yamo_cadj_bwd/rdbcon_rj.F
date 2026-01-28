
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_bwd/rdbcon_rj.F,v 1.1 2013/05/14 16:19:14 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE RDBCON_RJ ( JDATE, JTIME, TSTEP, LVL, BCON_RJ )

C-----------------------------------------------------------------------
C Function:
C   Read boundary values for RHOJ for advection
C   Load Air Density X SqRDMT = RHOJ into BCON_RJ for advection

C Revision History:
C   Apr 2013 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE

C Includes:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP ! timestep
      INTEGER, INTENT(IN) :: LVL ! layer
      REAL, INTENT(OUT)   :: BCON_RJ( : ) ! boundary values for RHOJ

C File variables:

      REAL :: RHOJ_BND ( NBNDY,NLAYS ) ! mid-layer boundary RhoJ

      REAL, ALLOCATABLE, SAVE :: RHOJ_BUF ( :,: ) ! mid-layer grid RhoJ

C External Functions not declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC
#ifdef parallel
      LOGICAL, EXTERNAL :: INTERPB
#endif

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'RDBCON_RJ'
      CHARACTER( 16 ) :: VNAME

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER :: MDATE          ! mid-advection date
      INTEGER :: MTIME          ! mid-advection time
      INTEGER :: STEP           ! advection time step in seconds

      INTEGER :: COL, ROW       ! loop counters
      INTEGER :: ALLOCSTAT

      INTEGER :: COUNT          ! Counter for constructing boundary arrays

      INTEGER :: GXOFF, GYOFF   ! global origin offset from file
      LOGICAL, SAVE :: WINDOW = .FALSE. ! posit same file and global
                                        ! processing domain
C for INTERPX
      INTEGER, SAVE :: STRTCOL, ENDCOL, STRTROW, ENDROW

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &        STRTCOL, ENDCOL, STRTROW, ENDROW )

C currently not implemented: case where only one origin component matches file's
         IF ( GXOFF .NE. 0 .AND. GYOFF .NE. 0 ) THEN
            WINDOW = .TRUE.     ! windowing from file
            STRTCOL = STRTCOL - 1
            ENDCOL  = ENDCOL  + 1
            STRTROW = STRTROW - 1
            ENDROW  = ENDROW  + 1
         END IF

         IF ( WINDOW ) THEN
            ALLOCATE ( RHOJ_BUF( 0:MY_NCOLS+1,0:MY_NROWS+1 ),
     &           STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating RHOJ_BUF'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF
      END IF                    ! FIRSTIME

      MDATE  = JDATE
      MTIME  = JTIME
      STEP   = TIME2SEC( TSTEP )
      CALL NEXTIME( MDATE, MTIME, SEC2TIME( STEP / 2 ) )

C Read & interpolate boundary RhoJ

      BCON_RJ = 0.0

      IF ( WINDOW ) THEN

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPX ( MET_CRO_3D, VNAME, PNAME,
     &        STRTCOL,ENDCOL, STRTROW,ENDROW, LVL,LVL,
     &        MDATE, MTIME, RHOJ_BUF ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

C Fill in DENSJ array for boundaries

         COUNT = 0
         DO ROW = 0, 0          ! South
            DO COL = 1, MY_NCOLS+1
               COUNT = COUNT + 1
               BCON_RJ ( COUNT ) = RHOJ_BUF ( COL,ROW )
            END DO
         END DO
         DO ROW = 1, MY_NROWS+1 ! East
            DO COL = MY_NCOLS+1, MY_NCOLS+1
               COUNT = COUNT + 1
               BCON_RJ ( COUNT ) = RHOJ_BUF ( COL,ROW )
            END DO
         END DO
         DO ROW = MY_NROWS+1, MY_NROWS+1 ! North
            DO COL = 0, MY_NCOLS
               COUNT = COUNT + 1
               BCON_RJ ( COUNT ) = RHOJ_BUF ( COL,ROW )
            END DO
         END DO
         DO ROW = 0, MY_NROWS   ! West
            DO COL = 0, 0
               COUNT = COUNT + 1
               BCON_RJ ( COUNT ) = RHOJ_BUF ( COL,ROW )
            END DO
         END DO

      ELSE

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPB ( MET_BDY_3D, VNAME, PNAME,
     &        MDATE, MTIME, NBNDY*NLAYS,
     &        RHOJ_BND ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_BDY_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

         BCON_RJ ( 1:NBNDY ) = RHOJ_BND( 1:NBNDY,LVL )

      END IF                    ! WINDOW

      RETURN
      END
