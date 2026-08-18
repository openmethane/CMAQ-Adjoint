
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_fwd/hadv.F,v 1.1 2013/05/14 15:53:42 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HADV_BWD(LGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)

C-----------------------------------------------------------------------
C Function:
C   Advection in the horizontal plane

C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.

C Revision History:
C    Feb 2013 by Peter Percell at UH-IMAQS:
C    - Modified hadvyppm.F for use with adjoint/4dvar code
C
C    May 2019 by Shunliu Zhao at Carleton:
C    - Modified for the adjoint
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines
      USE SUBST_MODULES         ! stenex
      USE ADV_CADJ

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

C Parameters:
      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV ! number of advected species
     &     + N_AE_ADV
     &     + N_NR_ADV
     &     + N_TR_ADV
     &     + 1                  ! for advecting RHO*SqRDMT

C Arguments:
      REAL, POINTER :: LGRID(:,:,:,:) ! concentrations
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
      INTEGER, INTENT(IN) :: ASTEP(:) ! advection time steps for layers (HHMMSS)
      LOGICAL, INTENT(IN OUT) :: XFIRST(:)

C External Functions:
      INTEGER, EXTERNAL :: TIME2SEC

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER      GXOFF, GYOFF ! global origin offset from file
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER      EDATE, ETIME

      INTEGER, SAVE :: ASPC     ! RHOJ (air) index in CGRID

      LOGICAL, SAVE :: W_BND_PE, E_BND_PE, S_BND_PE, N_BND_PE

      CHARACTER(16) :: PNAME = 'HADV_BWD'
      CHARACTER(96) :: XMSG = ' '

      REAL, ALLOCATABLE, SAVE :: ADJFAC(:, :) ! Yamo's adjustment
      REAL, ALLOCATABLE, SAVE :: XVEL(:, :) ! adjusted velocities along rows for substeps
      REAL, ALLOCATABLE, SAVE :: YVEL(:, :) ! adjusted velocities along columns for substeps

      INTEGER :: ALLOCSTAT

      LOGICAL :: XFIRST_L

      INTEGER :: SYNC           ! synchronization step (sec)
      INTEGER :: ADATE          ! current advection date, coded YYYYDDD
      INTEGER :: ATIME          ! current advection time, coded HHMMSS
      INTEGER :: LAYER          ! layer loop counter
      INTEGER :: ASTEP_L        ! advection step for layer and sync step (HHMMSS)
      INTEGER :: ASTEP_SEC      ! advection step for layer and sync step (sec)
      INTEGER :: NADVSS         ! number of advection sub-steps per sync step
      INTEGER :: SS             ! loop counter

      REAL :: BCON_RJ(1:NBNDY)     ! boundary conditions for RHOJ
      REAL :: BCON_ADJ(1:NBNDY, 1:N_SPC_ADV) ! boundary conditions for lgrid

      REAL :: RHOJ_BSS(1:NCOLS, 1:NROWS) ! RHOJ at beginning of advection sub-step
      REAL :: RHOJM( NCOLS,NROWS,NLAYS ) ! RhoJ (Kg/m**2) from Met file

      INTERFACE
         SUBROUTINE RDBCON_RJ(ADATE, ATIME, ASTEP_L, LAYER, BCON_RJ)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: ADATE, ATIME
            INTEGER, INTENT(IN) :: ASTEP_L
            INTEGER, INTENT(IN) :: LAYER
            REAL, INTENT(OUT)   :: BCON_RJ(:)
         END SUBROUTINE RDBCON_RJ
         SUBROUTINE RDBCON_ADJ(ADATE, ATIME, ASTEP_L, LAYER, BCON_ADJ)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: ADATE, ATIME
            INTEGER, INTENT(IN) :: ASTEP_L
            INTEGER, INTENT(IN) :: LAYER
            REAL, INTENT(OUT)   :: BCON_ADJ(:,:)
         END SUBROUTINE RDBCON_ADJ
         SUBROUTINE XADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &     BCON_RJ, XVEL, ADJFAC)
            IMPLICIT NONE
            REAL, INTENT(IN OUT) :: ADV_RHOJ(:, :, :, :)
            INTEGER, INTENT(IN)  :: ADATE, ATIME
            INTEGER, INTENT(IN)  :: ASTEP_L
            INTEGER, INTENT(IN)  :: LAYER
            REAL, INTENT(IN)     :: BCON_RJ(:)
            REAL, INTENT(OUT)    :: XVEL(:, :)
            REAL, INTENT(IN), OPTIONAL :: ADJFAC(0:, 0:)
         END SUBROUTINE XADV_RJ_FWD
         SUBROUTINE YADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &     BCON_RJ, YVEL, ADJFAC)
            IMPLICIT NONE
            REAL, INTENT(IN OUT) :: ADV_RHOJ(:, :, :, :)
            INTEGER, INTENT(IN)  :: ADATE, ATIME
            INTEGER, INTENT(IN)  :: ASTEP_L
            INTEGER, INTENT(IN)  :: LAYER
            REAL, INTENT(IN)     :: BCON_RJ(:)
            REAL, INTENT(OUT)    :: YVEL(:, :)
            REAL, INTENT(IN), OPTIONAL :: ADJFAC(0:, 0:)
         END SUBROUTINE YADV_RJ_FWD
         SUBROUTINE XADV_BWD(LGRID, ASTEP_L, LAYER, XVEL, BCON_ADJ)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: ASTEP_L
            INTEGER, INTENT(IN) :: LAYER
            REAL, INTENT(IN)    :: XVEL(:, :)
            REAL, INTENT(IN)    :: BCON_ADJ(:, :)
         END SUBROUTINE XADV_BWD
         SUBROUTINE YADV_BWD(LGRID, ASTEP_L, LAYER, YVEL, BCON_ADJ)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: ASTEP_L
            INTEGER, INTENT(IN) :: LAYER
            REAL, INTENT(IN)    :: YVEL(:, :)
            REAL, INTENT(IN)    :: BCON_ADJ(:, :)
         END SUBROUTINE YADV_BWD
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         W_BND_PE = ( COLSX_PE( 1,MYPE+1 ) .EQ. 1 )
         E_BND_PE = ( COLSX_PE( 2,MYPE+1 ) .EQ. GL_NCOLS )
         S_BND_PE = ( ROWSX_PE( 1,MYPE+1 ) .EQ. 1 )
         N_BND_PE = ( ROWSX_PE( 2,MYPE+1 ) .EQ. GL_NROWS )

         ALLOCATE(
     &        ADJFAC(0:NCOLS+1, 0:NROWS+1),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating ADJFAC'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

         ALLOCATE(
     &        XVEL(1:NCOLS+1, 1:NROWS),
     &        YVEL(1:NCOLS, 1:NROWS+1),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating XVEL, YVEL'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

         IF ( W_BND_PE ) ADJFAC(0, 1:NROWS) = 1.0
         IF ( E_BND_PE ) ADJFAC(NCOLS+1, 1:NROWS) = 1.0
         IF ( S_BND_PE ) ADJFAC(1:NCOLS, 0) = 1.0
         IF ( N_BND_PE ) ADJFAC(1:NCOLS, NROWS+1) = 1.0

! get file local domain offsets
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &        STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         CALL ADV_CADJ_INIT( JDATE, JTIME)
      END IF                    ! FIRSTIME

! Get rho*J at end of sync step
      EDATE = JDATE
      ETIME = JTIME
      CALL NEXTIME( EDATE, ETIME, TSTEP( 2 ) )

      IF ( .NOT. INTERPX( MET_CRO_3D, 'DENSA_J', PNAME,
     &     STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &     EDATE, ETIME, RHOJM ) ) THEN
         XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, EDATE, ETIME, XMSG, XSTAT1 )
      END IF

      ADV_RHOJ(:,:,:,1)=RHOJM(:,:,:)

      SYNC = TIME2SEC(TSTEP(2))
      DO LAYER = 1, NLAYS       ! LAYERS
         ADATE = EDATE
         ATIME = ETIME
         ASTEP_L = ASTEP(LAYER)
         CALL NEXTIME( ADATE, ATIME, - ASTEP_L)
         ASTEP_SEC = TIME2SEC(ASTEP_L)

C Note that advstep.F ensures that ASTEP_SEC is an integer divisor of SYNC
         NADVSS = SYNC / ASTEP_SEC
         DO SS = 1, NADVSS      ! SUBSTEPS
            CALL RDBCON_RJ(ADATE, ATIME, ASTEP_L, LAYER, BCON_RJ)
            CALL RDBCON_ADJ(ADATE, ATIME, ASTEP_L, LAYER, BCON_ADJ) !slz

            RHOJ_BSS(1:NCOLS, 1:NROWS) = ADV_RHOJ(1:NCOLS, 1:NROWS, LAYER, 1)

            IF ( XFIRST_L ) THEN
               CALL XADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &              BCON_RJ, XVEL(1:NCOLS+1, 1:NROWS))
               CALL XADV_BWD(LGRID, ASTEP_L, LAYER,
     &              XVEL(1:NCOLS+1, 1:NROWS), BCON_ADJ)
            ELSE
               CALL YADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &              BCON_RJ, YVEL(1:NCOLS, 1:NROWS+1))
               CALL YADV_BWD(LGRID, ASTEP_L, LAYER,
     &              YVEL(1:NCOLS, 1:NROWS+1), BCON_ADJ)
            END IF

            ADJFAC(1:NCOLS, 1:NROWS) =
     &           RHOJ_BSS(1:NCOLS, 1:NROWS)
     &           / ADV_RHOJ(1:NCOLS, 1:NROWS, LAYER, 1)

            CALL SUBST_COMM(ADJFAC, DSPL_N1_E1_S1_W1, DRCN_N_E_S_W, '1 0 2 0')

            IF ( XFIRST_L ) THEN
               CALL YADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &              BCON_RJ, YVEL(1:NCOLS, 1:NROWS+1), ADJFAC)
               CALL YADV_BWD(LGRID, ASTEP_L, LAYER,
     &              YVEL(1:NCOLS, 1:NROWS+1), BCON_ADJ)
            ELSE
               CALL XADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &              BCON_RJ, XVEL(1:NCOLS+1, 1:NROWS), ADJFAC)
               CALL XADV_BWD(LGRID, ASTEP_L, LAYER,
     &              XVEL(1:NCOLS+1, 1:NROWS), BCON_ADJ)
            END IF

            XFIRST_L = ( .NOT. XFIRST_L )

            CALL NEXTIME(ADATE, ATIME, - ASTEP_L)
         END DO                 ! SUBSTEPS
      END DO                    ! LAYERS

      RETURN
      END
