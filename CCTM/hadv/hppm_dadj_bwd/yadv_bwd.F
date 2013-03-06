
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/hppm_dadj_bwd/yadv_bwd.F,v 1.1 2011/02/21 23:54:17 ppercell Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE YADV_BWD(LGRID, JDATE, JTIME, TSTEP, NADVSS)

C-----------------------------------------------------------------------
C Function:
C   Discrete adjoint for advection in the horizontal plane; x2-direction.

C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.  All boundary input variables have the
C   same boundary perimeter structure with a thickness of 1

C Revision history:
C   Sep 2010 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines
      USE SUBST_MODULES         ! stenex

#ifdef parallel
      USE SWAP_SANDIA
#endif

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

C Parameters:
      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV ! advected species dimension
     &     + N_AE_ADV
     &     + N_NR_ADV
     &     + N_TR_ADV
     &     + 1                  ! for advecting RHO*SqRDMT

#ifdef parallel
      INTEGER, PARAMETER :: SWP = 3
#else
      INTEGER, PARAMETER :: SWP = 1
#endif

C Arguments:
      REAL, POINTER :: LGRID(:, :, :, :)
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
      INTEGER, INTENT(IN) :: NADVSS(:) ! number of advection sub-steps per sync step

C External Functions not declared in IODECL3.EXT:
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC
      REAL,    EXTERNAL :: ZFDBC

C File Variables:
      REAL :: ADV_CHK_SYNC(1:NCOLS, 1:NROWS, 1:NLAYS, 1:N_SPC_ADV) ! adv checkpt at sync step
      REAL :: VHAT(1:NCOLS+1, 1:NROWS+1) ! x2-component CX-velocity
      REAL :: BCON(1:NBNDY, 1:N_SPC_ADV) ! boundary concentrations

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: PNAME = 'YADV_BWD'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(16) :: UORV = 'X2VEL' ! flag for HCONTVEL

      INTEGER, SAVE :: SFX      ! fixed parameter for southern boundary
      INTEGER, SAVE :: NFX      ! fixed parameter for northern boundary

      INTEGER, SAVE :: ADV_MAP(1:N_SPC_ADV) ! global adv map to LGRID

      LOGICAL, SAVE :: BNDY_LO_PE, BNDY_HI_PE

      INTEGER, SAVE :: MX_NSS = 0 ! max number of advection sub-steps per sync step

      INTEGER :: ALLOCSTAT

      INTEGER :: SYNC           ! sync step (sec)
      INTEGER :: ASTEP          ! advection time step for layer (HHMMSS)
      INTEGER :: SS             ! loop counter

      REAL :: DT                ! advection time step for layer (sec)
      INTEGER :: DT_INT         ! advection time step for layer (sec)

      INTEGER :: ADATE          ! interpolation date
      INTEGER :: ATIME          ! interpolation time

      INTEGER :: C, R, L, S, V  ! loop counters
      INTEGER :: A2C

      REAL, ALLOCATABLE, SAVE, TARGET :: ADV_CHK(:, :, :, :)
#ifdef parallel
      REAL, POINTER :: SEND(:, :, :)
      REAL, POINTER :: RECV(:, :, :)
#endif
      REAL, ALLOCATABLE, SAVE :: VEL(:, :, :)

      REAL :: CON(1:NCOLS, 1-SWP:NROWS+SWP, 1:N_SPC_ADV)
      REAL, TARGET :: LAMBDA(1:NCOLS, 1-SWP:NROWS+SWP, 1:N_SPC_ADV)

C Required interfaces
      INTERFACE
         SUBROUTINE RDBCON(JDATE, JTIME, ASTEP, LVL, BCON)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: JDATE, JTIME, ASTEP, LVL
            REAL, INTENT(OUT)   :: BCON(:, :)
         END SUBROUTINE RDBCON
         SUBROUTINE HCONTVEL(JDATE, JTIME, ASTEP, LVL, UORV, UHAT)
            IMPLICIT NONE
            INTEGER, INTENT(IN)       :: JDATE, JTIME, ASTEP, LVL
            CHARACTER(16), INTENT(IN) :: UORV
            REAL, INTENT(OUT)         :: UHAT(:, :)
         END SUBROUTINE HCONTVEL
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         SFX = 0
         NFX = NCOLS + NROWS + 3

C Get CGRID offsets
         CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

C Create global map to LGRID
         S = 0
         DO V = 1, N_GC_ADV
            S = S + 1
            ADV_MAP(S) = GC_STRT - 1 + GC_ADV_MAP(V)
         END DO
         DO V = 1, N_AE_ADV
            S = S + 1
            ADV_MAP(S) = AE_STRT - 1 + AE_ADV_MAP(V)
         END DO
         DO V = 1, N_NR_ADV
            S = S + 1
            ADV_MAP(S) = NR_STRT - 1 + NR_ADV_MAP(V)
         END DO
         DO V = 1, N_TR_ADV
            S = S + 1
            ADV_MAP(S) = TR_STRT - 1 + TR_ADV_MAP(V)
         END DO

         ADV_MAP(N_SPC_ADV) = N_GC_SPCD

         CALL SUBST_HI_LO_BND_PE('R', BNDY_LO_PE, BNDY_HI_PE)

C Open checkpoint file as read-only
         CALL CHK_FILE_OPEN ( ADJ_YADV_CHK, JDATE, JTIME, .FALSE. )

      END IF                    ! FIRSTIME

      IF ( MAXVAL(NADVSS(1:NLAYS)) > MX_NSS ) THEN
         IF ( ALLOCATED(ADV_CHK) ) DEALLOCATE(ADV_CHK)
         IF ( ALLOCATED(VEL) ) DEALLOCATE(VEL)
         MX_NSS = MAXVAL(NADVSS(1:NLAYS))
         ALLOCATE(
     &        ADV_CHK(1:NCOLS, 1-SWP:NROWS+SWP, 1:N_SPC_ADV, 1:MX_NSS),
     &        VEL(1:NCOLS, 1:NROWS+1, 1:MX_NSS),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating ADV_CHK, VEL'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF
      END IF

C Read horizontal advection checkpoint at sync step
      CALL CHK_FILE_READ ( ADJ_YADV_CHK, JDATE, JTIME, ADV_CHK_SYNC )

      SYNC = TIME2SEC(TSTEP(2))
      DO L = 1, NLAYS

C Load the interior concentrations into ADV_CHK for the first sub-step
         ADV_CHK(1:NCOLS, 1:NROWS, 1:N_SPC_ADV, 1) =
     &        ADV_CHK_SYNC(1:NCOLS, 1:NROWS, L, 1:N_SPC_ADV)

C Get the data needed for the sub-steps of the reverse solution:
C   1. wind velocities at the middle of each sub-step
C   2. concentrations at the beginning of each sub-step
C If NADVSS(L) > 1, then it is necessary to redo sub-steps of the forward
C solution, because global checkpoints were only saved at the beginning of
C sync steps.
         ADATE = JDATE
         ATIME = JTIME
         DO SS = 1, NADVSS(L)

C ASTEP must be an integer in order to make the IO/API happy. Also, the sum of
C the NADVSS(L) steps must equal SYNC. We would like to use equal steps, but
C this is not possible unless SYNC is an integer multiple of NADVSS(L). In
C general, with these restrictions, the best we can do is to make all the steps
C differ from DT = REAL(SYNC) / NADVSS(L) by at most 1 second, i.e., replace
C DT by either either FLOOR(DT) = SYNC / NADVSS(L) or CEILING(DT) = FLOOR(DT) + 1.
            DT_INT = SYNC / NADVSS(L)
            IF ( MOD(SYNC, NADVSS(L)) /= 0 .AND.
     &           SS <= MOD(SYNC, NADVSS(L)) ) THEN
               DT_INT = DT_INT + 1
            END IF
            DT = REAL(DT_INT)
            ASTEP = SEC2TIME(DT_INT)

C Get the contravariant x2 velocity component
            CALL HCONTVEL(ADATE, ATIME, ASTEP, L, UORV, VHAT)

#ifdef parallel
C swap2d interface ...
C     swap2d( 1st value of array to be sent,
C             1st value of array to be received,
C             number of values to send in 1st dimension of array,
C             number of values to send in 2nd dimension of array,
C             size of 1st dimension in arrays,
C             direction to receive from: NSEW )

            CALL SWAP2D(VHAT(1, 1), VHAT(1, NROWS+1), NCOLS, 1, NCOLS+1, NORTH)
#endif

            VEL(1:NCOLS, 1:NROWS+1, SS) = VHAT(1:NCOLS, 1:NROWS+1)

C Get the boundary concentration values (converted to transport units)
            CALL RDBCON(ADATE, ATIME, ASTEP, L, BCON)

C Insert the boundary conditions into ADV_CHK
            IF ( BNDY_LO_PE ) THEN ! south boundary
               DO S = 1, N_SPC_ADV
                  DO C = 1, NCOLS
                     IF ( VHAT(C, 1) < 0.0 ) THEN ! outflow
                        ADV_CHK(C, 1-SWP:0, S, SS) = ZFDBC(
     &                       ADV_CHK(C, 1, S, SS), ADV_CHK(C, 2, S, SS),
     &                       VHAT(C, 1), VHAT(C, 2) )
                     ELSE       ! inflow
                        ADV_CHK(C, 1-SWP:0, S, SS) = BCON(SFX+C, S)
                     END IF
                  END DO
               END DO
            END IF

            IF ( BNDY_HI_PE ) THEN ! north boundary
               DO S = 1, N_SPC_ADV
                  DO C = 1, NCOLS
                     IF ( VHAT(C, NROWS+1) > 0.0 ) THEN ! outflow
                        ADV_CHK(C, NROWS+1:NROWS+SWP, S, SS) = ZFDBC(
     &                       ADV_CHK(C, NROWS, S, SS), ADV_CHK(C, NROWS-1, S, SS),
     &                       VHAT(C, NROWS+1), VHAT(C, NROWS) )
                     ELSE       ! inflow
                        ADV_CHK(C, NROWS+1:NROWS+SWP, S, SS) = BCON(NFX+C, S)
                     END IF
                  END DO
               END DO
            END IF

#ifdef parallel
C swap3d interface ...
C     swap3d( 1st value of array to be sent,
C             1st value of array to be received,
C             number of values to send in 1st dimension of array,
C             number of values to send in 2nd dimension of array,
C             number of values to send in 3rd dimension of array,
C             size of 1st dimension in arrays,
C             size of 2nd dimension in arrays,
C             direction to receive from: NSEW )

            SEND => ADV_CHK(1:NCOLS, 1:SWP, 1:N_SPC_ADV, SS)
            RECV => ADV_CHK(1:NCOLS, NROWS+1:NROWS+SWP, 1:N_SPC_ADV, SS)
            CALL SWAP3D(SEND, RECV, NCOLS, SWP, N_SPC_ADV,
     &           NCOLS, SWP, NORTH)

            SEND => ADV_CHK(1:NCOLS, NROWS-SWP+1:NROWS, 1:N_SPC_ADV, SS)
            RECV => ADV_CHK(1:NCOLS, 1-SWP:0, 1:N_SPC_ADV, SS)
            CALL SWAP3D(SEND, RECV, NCOLS, SWP, N_SPC_ADV,
     &           NCOLS, SWP, SOUTH)
#endif

            IF ( SS == NADVSS(L) ) EXIT

C Do the gridded computation for horizontal advection using the PPM scheme
            CON(:, :, :) = ADV_CHK(:, :, :, SS)
            DO C = 1, NCOLS     ! x2-direction
               CALL HPPM(NROWS, N_SPC_ADV, CON(C, :, :), VHAT(C, :),
     &              DT, REAL(YCELL_GD), 'R')
            END DO
            ADV_CHK(:, :, :, SS+1) = CON(:, :, :)

            CALL NEXTIME(ADATE, ATIME, ASTEP)

         END DO                 ! SS

C Reverse solution

C Load the interior sensitivities into LAMBDA
         DO S = 1, N_SPC_ADV
            A2C = ADV_MAP(S)
            LAMBDA(1:NCOLS, 1:NROWS, S) = LGRID(1:NCOLS, 1:NROWS, L, A2C)
         END DO

C Insert the boundary conditions into LAMBDA
         IF ( BNDY_LO_PE ) THEN ! south boundary
            LAMBDA(:, 1-SWP:0, :) = 0.0
         END IF

         IF ( BNDY_HI_PE ) THEN ! north boundary
            LAMBDA(:, NROWS+1:NROWS+SWP, :) = 0.0
         END IF

         DO SS = NADVSS(L), 1, -1
            DT_INT = SYNC / NADVSS(L)
            IF ( MOD(SYNC, NADVSS(L)) /= 0 .AND.
     &           SS <= MOD(SYNC, NADVSS(L)) ) THEN
               DT_INT = DT_INT + 1
            END IF
            DT = REAL(DT_INT)

#ifdef parallel
C swap3d interface ...
C     swap3d( 1st value of array to be sent,
C             1st value of array to be received,
C             number of values to send in 1st dimension of array,
C             number of values to send in 2nd dimension of array,
C             number of values to send in 3rd dimension of array,
C             size of 1st dimension in arrays,
C             size of 2nd dimension in arrays,
C             direction to receive from: NSEW )

            SEND => LAMBDA(1:NCOLS, 1:SWP, 1:N_SPC_ADV)
            RECV => LAMBDA(1:NCOLS, NROWS+1:NROWS+SWP, 1:N_SPC_ADV)
            CALL SWAP3D(SEND, RECV, NCOLS, SWP, N_SPC_ADV,
     &           NCOLS, SWP, NORTH)

            SEND => LAMBDA(1:NCOLS, NROWS-SWP+1:NROWS, 1:N_SPC_ADV)
            RECV => LAMBDA(1:NCOLS, 1-SWP:0, 1:N_SPC_ADV)
            CALL SWAP3D(SEND, RECV, NCOLS, SWP, N_SPC_ADV,
     &           NCOLS, SWP, SOUTH)
#endif

C Do the gridded computation for horizontal advection using the PPM scheme
            DO C = 1, NCOLS     ! x2-direction
               CALL HPPM_BWD(NROWS, N_SPC_ADV,
     &              ADV_CHK(C, :, :, SS), VEL(C, :, SS),
     &              DT, REAL(YCELL_GD), 'R', LAMBDA(C, :, :))
            END DO
         END DO                 ! SS

C Transfer updated interior sensitivities from LAMBDA into LGRID
         DO S = 1, N_SPC_ADV
            A2C = ADV_MAP(S)
            LGRID(1:NCOLS, 1:NROWS, L, A2C) = LAMBDA(1:NCOLS, 1:NROWS, S)
         END DO

      END DO                    ! NLAYS

      RETURN
      END
