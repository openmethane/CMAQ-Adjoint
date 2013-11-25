
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
      SUBROUTINE HADV(CGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)

C-----------------------------------------------------------------------
C Function:
C   Advection in the horizontal plane

C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.

C Revision History:
C    Feb 2013 by Peter Percell at UH-IMAQS:
C    - Modified hadvyppm.F for use with adjoint/4dvar code
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines
      USE SUBST_MODULES         ! stenex

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction

C Parameters:
      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV ! number of advected species
     &     + N_AE_ADV
     &     + N_NR_ADV
     &     + N_TR_ADV
     &     + 1                  ! for advecting RHO*SqRDMT

C Arguments:
      REAL, POINTER :: CGRID(:,:,:,:) ! concentrations
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

      INTEGER, SAVE :: ASPC     ! RHOJ (air) index in CGRID

      LOGICAL, SAVE :: W_BND_PE, E_BND_PE, S_BND_PE, N_BND_PE

      CHARACTER(16) :: PNAME = 'HADV'
      CHARACTER(96) :: XMSG = ' '

      REAL, ALLOCATABLE, SAVE :: ADJFAC(:, :) ! Yamo's adjustment

      INTEGER :: ALLOCSTAT

      INTEGER :: SYNC           ! synchronization step (sec)
      INTEGER :: ADATE          ! current advection date, coded YYYYDDD
      INTEGER :: ATIME          ! current advection time, coded HHMMSS
      INTEGER :: LAYER          ! layer loop counter
      INTEGER :: ASTEP_L        ! advection step for layer and sync step (HHMMSS)
      INTEGER :: ASTEP_SEC      ! advection step for layer and sync step (sec)
      INTEGER :: NADVSS         ! number of advection sub-steps per sync step
      INTEGER :: SS             ! loop counter

      REAL :: BCON(1:NBNDY, 1:N_SPC_ADV) ! boundary concentrations (in transport units)

      REAL :: RHOJ_BSS(1:NCOLS, 1:NROWS) ! RHOJ at beginning of advection sub-step

      INTERFACE
         SUBROUTINE RDBCON(ADATE, ATIME, ASTEP_L, LAYER, BCON)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: ADATE, ATIME
            INTEGER, INTENT(IN) :: ASTEP_L
            INTEGER, INTENT(IN) :: LAYER
            REAL, INTENT(OUT)   :: BCON(:, :)
         END SUBROUTINE RDBCON
         SUBROUTINE XADV(CGRID, ADATE, ATIME, ASTEP_L, LAYER, BCON, ADJFAC)
            IMPLICIT NONE
            REAL, POINTER       :: CGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: ADATE, ATIME
            INTEGER, INTENT(IN) :: ASTEP_L
            INTEGER, INTENT(IN) :: LAYER
            REAL, INTENT(IN)    :: BCON(:, :)
            REAL, INTENT(IN), OPTIONAL :: ADJFAC(0:, 0:)
         END SUBROUTINE XADV
         SUBROUTINE YADV(CGRID, ADATE, ATIME, ASTEP_L, LAYER, BCON, ADJFAC)
            IMPLICIT NONE
            REAL, POINTER       :: CGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: ADATE, ATIME
            INTEGER, INTENT(IN) :: ASTEP_L
            INTEGER, INTENT(IN) :: LAYER
            REAL, INTENT(IN)    :: BCON(:, :)
            REAL, INTENT(IN), OPTIONAL :: ADJFAC(0:, 0:)
         END SUBROUTINE YADV
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

C Get CGRID offsets
         CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

         ASPC = GC_STRT - 1 + N_GC_SPCD

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

         IF ( W_BND_PE ) ADJFAC(0, 1:NROWS) = 1.0
         IF ( E_BND_PE ) ADJFAC(NCOLS+1, 1:NROWS) = 1.0
         IF ( S_BND_PE ) ADJFAC(1:NCOLS, 0) = 1.0
         IF ( N_BND_PE ) ADJFAC(1:NCOLS, NROWS+1) = 1.0

C Create new adjoint file and open for update
         CALL CHK_FILE_OPEN(ADJ_HA_RHOJ_CHK, JDATE, JTIME, .TRUE., TSTEP)
      END IF                    ! FIRSTIME

C Write horizontal advection checkpoint at sync step
      CALL CHK_FILE_WRITE(ADJ_HA_RHOJ_CHK, JDATE, JTIME,
     &     CGRID(1:NCOLS, 1:NROWS, 1:NLAYS, ASPC:ASPC))

      SYNC = TIME2SEC(TSTEP(2))
      DO LAYER = 1, NLAYS       ! LAYERS
         ADATE = JDATE
         ATIME = JTIME
         ASTEP_L = ASTEP(LAYER)
         ASTEP_SEC = TIME2SEC(ASTEP_L)

C Note that advstep.F ensures that ASTEP_SEC is an integer divisor of SYNC
         NADVSS = SYNC / ASTEP_SEC
         DO SS = 1, NADVSS      ! SUBSTEPS
            CALL RDBCON(ADATE, ATIME, ASTEP_L, LAYER, BCON)

            RHOJ_BSS(1:NCOLS, 1:NROWS) = CGRID(1:NCOLS, 1:NROWS, LAYER, ASPC)

            IF ( XFIRST(LAYER) ) THEN
               CALL XADV(CGRID, ADATE, ATIME, ASTEP_L, LAYER, BCON)
            ELSE
               CALL YADV(CGRID, ADATE, ATIME, ASTEP_L, LAYER, BCON)
            END IF

            ADJFAC(1:NCOLS, 1:NROWS) =
     &           RHOJ_BSS(1:NCOLS, 1:NROWS) / CGRID(1:NCOLS, 1:NROWS, LAYER, ASPC)

            CALL SUBST_COMM(ADJFAC, DSPL_N1_E1_S1_W1, DRCN_N_E_S_W, '1 0 2 0')

            IF ( XFIRST(LAYER) ) THEN
               CALL YADV(CGRID, ADATE, ATIME, ASTEP_L, LAYER, BCON, ADJFAC)
            ELSE
               CALL XADV(CGRID, ADATE, ATIME, ASTEP_L, LAYER, BCON, ADJFAC)
            END IF

            XFIRST(LAYER) = ( .NOT. XFIRST(LAYER) )

            CALL NEXTIME(ADATE, ATIME, ASTEP_L)
         END DO                 ! SUBSTEPS
      END DO                    ! LAYERS

      RETURN
      END
