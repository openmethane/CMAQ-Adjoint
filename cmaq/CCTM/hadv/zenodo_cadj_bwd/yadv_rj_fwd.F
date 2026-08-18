
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_bwd/yadv_rj_fwd.F,v 1.1 2013/05/14 16:19:14 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE YADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &     BCON_RJ, YVEL, ADJFAC)

C-----------------------------------------------------------------------
C Function:
C   Advect RHOJ along columns in the horizontal plane to get YVEL, adjusted
C   velocities along columns

C Revision history:
C   Apr 2013 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C-----------------------------------------------------------------------

      USE HGRD_DEFN             ! horizontal domain specifications

#ifdef parallel
      USE SWAP_SANDIA
#endif

      IMPLICIT NONE

C Parameters:
#ifdef parallel
      INTEGER, PARAMETER :: SWP = 3
#else
      INTEGER, PARAMETER :: SWP = 1
#endif

C Arguments:
      REAL, INTENT(IN OUT) :: ADV_RHOJ(:, :, :, :) ! advected RHOJ
      INTEGER, INTENT(IN)  :: ADATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN)  :: ATIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN)  :: ASTEP_L ! advection time step for layer (HHMMSS)
      INTEGER, INTENT(IN)  :: LAYER ! layer
      REAL, INTENT(IN)     :: BCON_RJ(:) ! boundary values for RHOJ
      REAL, INTENT(OUT)    :: YVEL(:, :) ! adjusted velocities along columns
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

      INTEGER :: C, R           ! loop counters

      REAL, TARGET :: RHOJ(1:NCOLS, 1-SWP:NROWS+SWP, 1:1)
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

      !slzx YVEL(1:NCOLS, 1:NROWS+1) = VHAT(1:NCOLS, 1:NROWS+1)
      YVEL(1:NCOLS, 1:NROWS+1) = - VHAT(1:NCOLS, 1:NROWS+1)

C Load the interior concentrations into RHOJ
      RHOJ(1:NCOLS, 1:NROWS, 1) = ADV_RHOJ(1:NCOLS, 1:NROWS, LAYER, 1)

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

      SEND => RHOJ(1:NCOLS, 1:SWP, 1:1)
      RECV => RHOJ(1:NCOLS, NROWS+1:NROWS+SWP, 1:1)
      CALL SWAP3D(SEND, RECV, NCOLS, SWP, 1, NCOLS, SWP, NORTH)

      SEND => RHOJ(1:NCOLS, NROWS-SWP+1:NROWS, 1:1)
      RECV => RHOJ(1:NCOLS, 1-SWP:0, 1:1)
      CALL SWAP3D(SEND, RECV, NCOLS, SWP, 1, NCOLS, SWP, SOUTH)
#endif

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

C Insert the boundary conditions into RHOJ
      DO C = 1, NCOLS
         IF ( BNDY_LO_PE ) THEN ! south boundary
            IF ( YVEL(C, 1) < 0.0 ) THEN ! outflow
               RHOJ(C, 1-SWP:0, 1) = ZFDBC(
     &              RHOJ(C, 1, 1), RHOJ(C, 2, 1),
     &              YVEL(C, 1), YVEL(C, 2) )
            ELSE                ! inflow
               RHOJ(C, 1-SWP:0, 1) = BCON_RJ(SFX+C)
            END IF
         END IF

         IF ( BNDY_HI_PE ) THEN ! north boundary
            IF ( YVEL(C, NROWS+1) > 0.0 ) THEN ! outflow
               RHOJ(C, NROWS+1:NROWS+SWP, 1) = ZFDBC(
     &              RHOJ(C, NROWS, 1), RHOJ(C, NROWS-1, 1),
     &              YVEL(C, NROWS+1), YVEL(C, NROWS) )
            ELSE                ! inflow
               RHOJ(C, NROWS+1:NROWS+SWP, 1) = BCON_RJ(NFX+C)
            END IF
         END IF

C Do the gridded computation for horizontal advection using the PPM scheme
         CALL HPPM(NROWS, 1, RHOJ(C, :, :), YVEL(C, :), DT, DS, 'R')
      END DO

C Transfer updated interior concentrations from RHOJ into ADV_RHOJ
      ADV_RHOJ(1:NCOLS, 1:NROWS, LAYER, 1) = RHOJ(1:NCOLS, 1:NROWS, 1)

      RETURN
      END
