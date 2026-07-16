
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_bwd/xadv_rj_fwd.F,v 1.1 2013/05/14 16:19:14 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE XADV_RJ_FWD(ADV_RHOJ, ADATE, ATIME, ASTEP_L, LAYER,
     &     BCON_RJ, XVEL, ADJFAC)

C-----------------------------------------------------------------------
C Function:
C   Advect RHOJ along rows in the horizontal plane to get XVEL, adjusted
C   velocities along rows

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
      REAL, INTENT(OUT)    :: XVEL(:, :) ! adjusted velocities along rows
      REAL, INTENT(IN), OPTIONAL :: ADJFAC(0:, 0:) ! Yamo's velocity adjustment factor

C External Functions not declared in IODECL3.EXT:
      INTEGER, EXTERNAL :: TIME2SEC
      REAL,    EXTERNAL :: ZFDBC

C File Variables:
      REAL :: UHAT(1:NCOLS+1, 1:NROWS+1) ! x1-component CX-velocity

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: UORV = 'X1VEL' ! flag for HCONTVEL

      INTEGER, SAVE :: EFX      ! fixed parameter for eastern boundary
      INTEGER, SAVE :: WFX      ! fixed parameter for western boundary

      LOGICAL, SAVE :: BNDY_LO_PE, BNDY_HI_PE

      REAL, SAVE :: DS          ! distance step

      REAL :: DT                ! advection time step for layer (sec)

      INTEGER :: C, R           ! loop counters

      REAL, TARGET :: RHOJ(1-SWP:NCOLS+SWP, 1:NROWS, 1:1)
#ifdef parallel
      REAL, POINTER :: SEND(:, :, :)
      REAL, POINTER :: RECV(:, :, :)
#endif

C Required interfaces
      INTERFACE
         SUBROUTINE HCONTVEL(ADATE, ATIME, ASTEP_L, LAYER, UORV, UHAT)
            IMPLICIT NONE
            INTEGER, INTENT(IN)       :: ADATE, ATIME, ASTEP_L, LAYER
            CHARACTER(16), INTENT(IN) :: UORV
            REAL, INTENT(OUT)         :: UHAT(:, :)
         END SUBROUTINE HCONTVEL
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         EFX = NCOLS + 1
         WFX = 2 * NCOLS + NROWS + 4

         BNDY_LO_PE = ( COLSX_PE(1, MYPE+1) == 1 )
         BNDY_HI_PE = ( COLSX_PE(2, MYPE+1) == GL_NCOLS )

         DS = REAL(XCELL_GD)
      END IF                    ! FIRSTIME

      DT = TIME2SEC(ASTEP_L)

C Get the contravariant x1 velocity component
      CALL HCONTVEL(ADATE, ATIME, ASTEP_L, LAYER, UORV, UHAT)

#ifdef parallel
C swap2d interface ...
C     swap2d( 1st value of array to be sent,
C             1st value of array to be received,
C             number of values to send in 1st dimension of array,
C             number of values to send in 2nd dimension of array,
C             size of 1st dimension in arrays,
C             direction to receive from: NSEW )

      CALL SWAP2D(UHAT(1, 1), UHAT(NCOLS+1, 1), 1, NROWS, NCOLS+1, EAST)
#endif

      !slzx XVEL(1:NCOLS+1, 1:NROWS) = UHAT(1:NCOLS+1, 1:NROWS)
      XVEL(1:NCOLS+1, 1:NROWS) = - UHAT(1:NCOLS+1, 1:NROWS)
	!slx negative velocities used for velocity adjustment

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

      SEND => RHOJ(1:SWP, 1:NROWS, 1:1)
      RECV => RHOJ(NCOLS+1:NCOLS+SWP, 1:NROWS, 1:1)
      CALL SWAP3D(SEND, RECV, SWP, NROWS, 1, SWP, NROWS, EAST)

      SEND => RHOJ(NCOLS-SWP+1:NCOLS, 1:NROWS, 1:1)
      RECV => RHOJ(1-SWP:0, 1:NROWS, 1:1)
      CALL SWAP3D(SEND, RECV, SWP, NROWS, 1, SWP, NROWS, WEST)
#endif

      IF ( PRESENT(ADJFAC) ) THEN
C Get Yamo's adjustment to velocity from the upwind cell
         DO R = 1, NROWS
            DO C = 1, NCOLS + 1
               IF ( XVEL(C, R) < 0.0 ) THEN
                  XVEL(C, R) = XVEL(C, R) * ADJFAC(C, R)
               ELSE IF ( XVEL(C, R) > 0.0 ) THEN
                  XVEL(C, R) = XVEL(C, R) * ADJFAC(C-1, R)
               END IF
            END DO
         END DO
      END IF

C Insert the boundary conditions into RHOJ
      DO R = 1, NROWS
         IF ( BNDY_LO_PE ) THEN ! west boundary
            IF ( XVEL(1, R) < 0.0 ) THEN ! outflow
               RHOJ(1-SWP:0, R, 1) = ZFDBC(
     &              RHOJ(1, R, 1), RHOJ(2, R, 1),
     &              XVEL(1, R), XVEL(2, R) )
            ELSE                ! inflow
               RHOJ(1-SWP:0, R, 1) = BCON_RJ(WFX+R)
            END IF
         END IF

         IF ( BNDY_HI_PE ) THEN ! east boundary
            IF ( XVEL(NCOLS+1, R) > 0.0 ) THEN ! outflow
               RHOJ(NCOLS+1:NCOLS+SWP, R, 1) = ZFDBC(
     &              RHOJ(NCOLS, R, 1), RHOJ(NCOLS-1, R, 1),
     &              XVEL(NCOLS+1, R), XVEL(NCOLS, R) )
            ELSE                ! inflow
               RHOJ(NCOLS+1:NCOLS+SWP, R, 1) = BCON_RJ(EFX+R)
            END IF
         END IF

C Do the gridded computation for horizontal advection using the PPM scheme
         CALL HPPM(NCOLS, 1, RHOJ(:, R, :), XVEL(:, R), DT, DS, 'C')
      END DO

C Transfer updated interior concentrations from RHOJ into ADV_RHOJ
      ADV_RHOJ(1:NCOLS, 1:NROWS, LAYER, 1) = RHOJ(1:NCOLS, 1:NROWS, 1)

      RETURN
      END
