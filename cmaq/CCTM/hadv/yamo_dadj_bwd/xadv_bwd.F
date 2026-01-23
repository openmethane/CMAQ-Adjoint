
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_bwd/xadv_bwd.F,v 1.1 2013/05/14 16:19:14 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE XADV_BWD(LGRID, ASTEP_L, LAYER, CON, XVEL)

C-----------------------------------------------------------------------
C Function:
C   Continuous adjoint for advection in the horizontal plane; x1-direction.

C Revision history:
C   June 2013 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
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
      REAL, INTENT(IN), TARGET :: CON(1-SWP:, 1:, 1:) ! advected concentrations
      REAL, INTENT(IN)    :: XVEL(:, :) ! adjusted velocities along rows

C External Functions not declared in IODECL3.EXT:
      INTEGER, EXTERNAL :: TIME2SEC

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: PNAME = 'XADV_BWD'
      CHARACTER(96) :: XMSG = ' '

      REAL, ALLOCATABLE, SAVE, TARGET :: LAMBDA(:, :, :)
      INTEGER :: ALLOCSTAT

      INTEGER, SAVE :: ADV_MAP(1:N_SPC_ADV) ! global adv map to CGRID

      REAL, SAVE :: DS          ! distance step

      REAL :: DT                ! advection time step for layer (sec)

      INTEGER :: C, R, S, V     ! loop counters
      INTEGER :: A2C

#ifdef parallel
      REAL, POINTER :: SEND(:, :, :)
      REAL, POINTER :: RECV(:, :, :)
#endif

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         ALLOCATE(
     &        LAMBDA(1-SWP:NCOLS+SWP, 1:NROWS, 1:N_SPC_ADV),
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

         DS = REAL(XCELL_GD)
      END IF                    ! FIRSTIME

      DT = TIME2SEC(ASTEP_L)

C Load the interior values of LGRID into LAMBDA
      DO S = 1, N_SPC_ADV
         A2C = ADV_MAP(S)
         LAMBDA(1:NCOLS, 1:NROWS, S) = LGRID(1:NCOLS, 1:NROWS, LAYER, A2C)
      END DO

C Set the exterior values of LAMBDA to zero
      LAMBDA(1-SWP:0, 1:NROWS, 1:N_SPC_ADV) = 0.0
      LAMBDA(NCOLS+1:NCOLS+SWP, 1:NROWS, 1:N_SPC_ADV) = 0.0

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

      SEND => CON(1:SWP, 1:NROWS, 1:N_SPC_ADV)
      RECV => CON(NCOLS+1:NCOLS+SWP, 1:NROWS, 1:N_SPC_ADV)
      CALL SWAP3D(SEND, RECV, SWP, NROWS, N_SPC_ADV,
     &     SWP, NROWS, EAST)

      SEND => CON(NCOLS-SWP+1:NCOLS, 1:NROWS, 1:N_SPC_ADV)
      RECV => CON(1-SWP:0, 1:NROWS, 1:N_SPC_ADV)
      CALL SWAP3D(SEND, RECV, SWP, NROWS, N_SPC_ADV,
     &     SWP, NROWS, WEST)
#endif

C Do the gridded computation for backward horizontal advection using
C the discrete adjoint PPM scheme
      DO R = 1, NROWS
         CALL HPPM_BWD(NCOLS, N_SPC_ADV, CON(:, R, :), LAMBDA(:, R, :),
     &        XVEL(:, R), DT, DS, 'C')
      END DO

C Transfer updated interior concentrations from LAMBDA into LGRID
      DO S = 1, N_SPC_ADV
         A2C = ADV_MAP(S)
         LGRID(1:NCOLS, 1:NROWS, LAYER, A2C) = LAMBDA(1:NCOLS, 1:NROWS, S)
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

C Get values of LAMBDA inside this processor's domain which were updated by
C neighboring processors
      SEND => LAMBDA(1-SWP:0, 1:NROWS, 1:N_SPC_ADV)
      RECV => LAMBDA(NCOLS-SWP+1:NCOLS, 1:NROWS, 1:N_SPC_ADV)
      CALL SWAP3D(SEND, RECV, SWP, NROWS, N_SPC_ADV,
     &     SWP, NROWS, EAST)

      SEND => LAMBDA(NCOLS+1:NCOLS+SWP, 1:NROWS, 1:N_SPC_ADV)
      RECV => LAMBDA(1:SWP, 1:NROWS, 1:N_SPC_ADV)
      CALL SWAP3D(SEND, RECV, SWP, NROWS, N_SPC_ADV,
     &     SWP, NROWS, WEST)

C Update LGRID with values of LAMBDA from neighboring processors
      IF ( COLSX_PE(1, MYPE+1) /= 1 ) THEN ! not west boundary
         DO S = 1, N_SPC_ADV
            A2C = ADV_MAP(S)
            LGRID(1:SWP, 1:NROWS, LAYER, A2C) =
     &           LGRID(1:SWP, 1:NROWS, LAYER, A2C) +
     &           LAMBDA(1:SWP, 1:NROWS, S)
         END DO
      END IF

      IF ( COLSX_PE(2, MYPE+1) /= GL_NCOLS ) THEN ! not east boundary
         DO S = 1, N_SPC_ADV
            A2C = ADV_MAP(S)
            LGRID(NCOLS-SWP+1:NCOLS, 1:NROWS, LAYER, A2C) =
     &           LGRID(NCOLS-SWP+1:NCOLS, 1:NROWS, LAYER, A2C) +
     &           LAMBDA(NCOLS-SWP+1:NCOLS, 1:NROWS, S)
         END DO
      END IF
#endif

      RETURN
      END
