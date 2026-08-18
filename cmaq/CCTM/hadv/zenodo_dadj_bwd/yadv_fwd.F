
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_fwd/yadv.F,v 1.1 2013/05/14 15:53:42 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE YADV_FWD(CON, ADATE, ATIME, ASTEP_L, LAYER, BCON, YVEL, ADJFAC)

C-----------------------------------------------------------------------
C Function:
C   Forward advection in the x2-direction of the horizontal plane; output
C   (adjusted) velocities along columns

C Revision history:
C   June 2013 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C-----------------------------------------------------------------------

      USE HGRD_DEFN             ! horizontal domain specifications

#ifdef parallel
      USE SWAP_SANDIA
#endif

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table

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
      REAL, INTENT(IN OUT), TARGET :: CON(1:, 1-SWP:, 1:) ! advected concentrations
      INTEGER, INTENT(IN) :: ADATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: ATIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: ASTEP_L ! advection time step for layer (HHMMSS)
      INTEGER, INTENT(IN) :: LAYER ! layer
      REAL, INTENT(IN)    :: BCON(:, :) ! boundary concentrations
      REAL, INTENT(OUT)   :: YVEL(:, :) ! adjusted velocities along columns
      REAL, INTENT(IN), OPTIONAL :: ADJFAC(0:, 0:) ! Yamo's velocity adjustment factor

C External Functions not declared in IODECL3.EXT:
      INTEGER, EXTERNAL :: TIME2SEC
      REAL,    EXTERNAL :: ZFDBC

C File Variables:
      REAL :: VHAT(1:NCOLS+1, 1:NROWS+1) ! x2-component CX-velocity

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: UORV = 'X2VEL' ! flag for HCONTVEL

      INTEGER, SAVE :: SFX      ! fixed parameter for southern boundary
      INTEGER, SAVE :: NFX      ! fixed parameter for northern boundary

      LOGICAL, SAVE :: BNDY_LO_PE, BNDY_HI_PE

      REAL, SAVE :: DS          ! distance step

      REAL :: DT                ! advection time step for layer (sec)

      INTEGER :: C, R, S        ! loop counters

#ifdef parallel
      REAL, POINTER :: SEND(:, :, :)
      REAL, POINTER :: RECV(:, :, :)
#endif

C Required interfaces
      INTERFACE
         SUBROUTINE HCONTVEL(ADATE, ATIME, ASTEP_L, LAYER, UORV, VHAT)
            IMPLICIT NONE
            INTEGER, INTENT(IN)       :: ADATE, ATIME, ASTEP_L, LAYER
            CHARACTER(16), INTENT(IN) :: UORV
            REAL, INTENT(OUT)         :: VHAT(:, :)
         END SUBROUTINE HCONTVEL
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         SFX = 0
         NFX = NCOLS + NROWS + 3

         BNDY_LO_PE = ( ROWSX_PE(1, MYPE+1) == 1 )
         BNDY_HI_PE = ( ROWSX_PE(2, MYPE+1) == GL_NROWS )

         DS = REAL(YCELL_GD)
      END IF                    ! FIRSTIME

      DT = TIME2SEC(ASTEP_L)

C Get the contravariant x2 velocity component
      CALL HCONTVEL(ADATE, ATIME, ASTEP_L, LAYER, UORV, VHAT)

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

      YVEL(1:NCOLS, 1:NROWS+1) = VHAT(1:NCOLS, 1:NROWS+1)

      IF ( PRESENT(ADJFAC) ) THEN
C Get Yamo's adjustment to velocity from the upwind cell
         DO R = 1, NROWS + 1
            DO C = 1, NCOLS
               IF ( YVEL(C, R) < 0.0 ) THEN
                  YVEL(C, R) = YVEL(C, R) * ADJFAC(C, R)
               ELSE IF ( YVEL(C, R) > 0.0 ) THEN
                  YVEL(C, R) = YVEL(C, R) * ADJFAC(C, R-1)
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

      SEND => CON(1:NCOLS, 1:SWP, 1:N_SPC_ADV)
      RECV => CON(1:NCOLS, NROWS+1:NROWS+SWP, 1:N_SPC_ADV)
      CALL SWAP3D(SEND, RECV, NCOLS, SWP, N_SPC_ADV,
     &     NCOLS, SWP, NORTH)

      SEND => CON(1:NCOLS, NROWS-SWP+1:NROWS, 1:N_SPC_ADV)
      RECV => CON(1:NCOLS, 1-SWP:0, 1:N_SPC_ADV)
      CALL SWAP3D(SEND, RECV, NCOLS, SWP, N_SPC_ADV,
     &     NCOLS, SWP, SOUTH)
#endif

C Insert the boundary conditions into CON
      DO C = 1, NCOLS
         IF ( BNDY_LO_PE ) THEN ! south boundary
            DO S = 1, N_SPC_ADV
               IF ( YVEL(C, 1) < 0.0 ) THEN ! outflow
                  CON(C, 1-SWP:0, S) = ZFDBC(
     &                 CON(C, 1, S), CON(C, 2, S),
     &                 YVEL(C, 1), YVEL(C, 2) )
               ELSE             ! inflow
                  CON(C, 1-SWP:0, S) = BCON(SFX+C, S)
               END IF
            END DO
         END IF

         IF ( BNDY_HI_PE ) THEN ! north boundary
            DO S = 1, N_SPC_ADV
               IF ( YVEL(C, NROWS+1) > 0.0 ) THEN ! outflow
                  CON(C, NROWS+1:NROWS+SWP, S) = ZFDBC(
     &                 CON(C, NROWS, S), CON(C, NROWS-1, S),
     &                 YVEL(C, NROWS+1), YVEL(C, NROWS) )
               ELSE             ! inflow
                  CON(C, NROWS+1:NROWS+SWP, S) = BCON(NFX+C, S)
               END IF
            END DO
         END IF

C Do the gridded computation for horizontal advection using the PPM scheme
         CALL HPPM(NROWS, N_SPC_ADV, CON(C, :, :), YVEL(C, :), DT, DS, 'R')
      END DO

      RETURN
      END
