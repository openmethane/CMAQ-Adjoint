
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_bwd/yadv_bwd.F,v 1.1 2013/05/14 16:19:14 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE YADV_BWD(LGRID, ASTEP_L, LAYER, YVEL, BCON_ADJ)

C-----------------------------------------------------------------------
C Function:
C   Continuous adjoint for advection in the horizontal plane; x2-direction.

C Revision history:
C   Apr 2013 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets

#ifdef parallel
      USE SWAP_SANDIA
#endif

      IMPLICIT NONE

C Include Files:
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
      INTEGER, INTENT(IN) :: ASTEP_L ! advection time step for layer (HHMMSS)
      INTEGER, INTENT(IN) :: LAYER ! layer
      REAL, INTENT(IN) :: YVEL(:, :) ! adjusted velocities along columns
      REAL, INTENT(IN) :: BCON_ADJ(:, :) ! boundary values

C External Functions not declared in IODECL3.EXT:
      INTEGER, EXTERNAL :: TIME2SEC
      REAL,    EXTERNAL :: ZFDBC

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: PNAME = 'YADV_BWD'
      CHARACTER(96) :: XMSG = ' '

      REAL, ALLOCATABLE, SAVE, TARGET :: LAMBDA(:, :, :)
      INTEGER :: ALLOCSTAT

      INTEGER, SAVE :: SFX      ! fixed parameter for southern boundary
      INTEGER, SAVE :: NFX      ! fixed parameter for northern boundary

      INTEGER, SAVE :: ADV_MAP(1:N_SPC_ADV) ! global adv map to CGRID

      REAL, SAVE :: DS          ! distance step

      REAL :: DT                ! advection time step for layer (sec)

      INTEGER :: C, R, S, V     ! loop counters
      INTEGER :: A2C

      REAL :: NEG_VEL(1:NCOLS, 1:NROWS+1) ! -YVEL

#ifdef parallel
      REAL, POINTER :: SEND(:, :, :)
      REAL, POINTER :: RECV(:, :, :)
#endif

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         SFX = 0
         NFX = NCOLS + NROWS + 3

         ALLOCATE(
     &        LAMBDA(1:NCOLS, 1-SWP:NROWS+SWP, 1:N_SPC_ADV),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating LAMBDA'
            CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT1)
         END IF

C Get CGRID offsets
         CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

C Create global map to CGRID
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

         DS = REAL(YCELL_GD)
      END IF                    ! FIRSTIME

      DT = TIME2SEC(ASTEP_L)

C Load the interior concentrations into LAMBDA
      DO S = 1, N_SPC_ADV
         A2C = ADV_MAP(S)
         LAMBDA(1:NCOLS, 1:NROWS, S) = LGRID(1:NCOLS, 1:NROWS, LAYER, A2C)
      END DO

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
     &     NCOLS, SWP, NORTH)

      SEND => LAMBDA(1:NCOLS, NROWS-SWP+1:NROWS, 1:N_SPC_ADV)
      RECV => LAMBDA(1:NCOLS, 1-SWP:0, 1:N_SPC_ADV)
      CALL SWAP3D(SEND, RECV, NCOLS, SWP, N_SPC_ADV,
     &     NCOLS, SWP, SOUTH)
#endif

C Do the gridded computation for horizontal advection using the PPM scheme
      !slzx NEG_VEL(1:NCOLS, 1:NROWS+1) = -YVEL(1:NCOLS, 1:NROWS+1)
      NEG_VEL(1:NCOLS, 1:NROWS+1) = YVEL(1:NCOLS, 1:NROWS+1)
      !slzx velocity direction revsersed in yadv_rj...
      DO C = 1, NCOLS
C Insert the boundary conditions into LAMBDA
         IF ( ROWSX_PE(1, MYPE+1) == 1 ) THEN ! south boundary
!slzx            LAMBDA(:, 1-SWP:0, :) = 0.0
            DO S = 1, N_SPC_ADV
               IF ( NEG_VEL(C, 1) < 0.0 ) THEN ! outflow
                  LAMBDA(C, 1-SWP:0, S) = ZFDBC(
     &                 LAMBDA(C, 1, S), LAMBDA(C, 2, S),
     &                 NEG_VEL(C, 1), NEG_VEL(C, 2) )
               ELSE             ! inflow
                  LAMBDA(C, 1-SWP:0, S) = BCON_ADJ(SFX+C, S)
               END IF
            END DO
         END IF

         IF ( ROWSX_PE(2, MYPE+1) == GL_NROWS ) THEN ! north boundary
!slzx            LAMBDA(:, NROWS+1:NROWS+SWP, :) = 0.0
            DO S = 1, N_SPC_ADV
               IF ( NEG_VEL(C, NROWS+1) > 0.0 ) THEN ! outflow
                  LAMBDA(C, NROWS+1:NROWS+SWP, S) = ZFDBC(
     &                 LAMBDA(C, NROWS, S), LAMBDA(C, NROWS-1, S),
     &                 NEG_VEL(C, NROWS+1), NEG_VEL(C, NROWS) )
               ELSE             ! inflow
                  LAMBDA(C, NROWS+1:NROWS+SWP, S) = BCON_ADJ(NFX+C, S)
               END IF
            END DO
         END IF

         CALL HPPM(NROWS, N_SPC_ADV, LAMBDA(C, :, :), NEG_VEL(C, :),
     &        DT, DS, 'R')
      END DO

C Transfer updated interior concentrations from LAMBDA into LGRID
      DO S = 1, N_SPC_ADV
         A2C = ADV_MAP(S)
         LGRID(1:NCOLS, 1:NROWS, LAYER, A2C) = LAMBDA(1:NCOLS, 1:NROWS, S)
      END DO

      RETURN
      END
