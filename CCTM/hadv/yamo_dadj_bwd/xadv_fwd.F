
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_fwd/xadv.F,v 1.1 2013/05/14 15:53:42 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE XADV_FWD(CON, ADATE, ATIME, ASTEP_L, LAYER, BCON, XVEL, ADJFAC)

C-----------------------------------------------------------------------
C Function:
C   Forward advection in the x1-direction of the horizontal plane; output
C   (adjusted) velocities along rows

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
      REAL, INTENT(IN OUT), TARGET :: CON(1-SWP:, 1:, 1:) ! advected concentrations
      INTEGER, INTENT(IN) :: ADATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: ATIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: ASTEP_L ! advection time step for layer (HHMMSS)
      INTEGER, INTENT(IN) :: LAYER ! layer
      REAL, INTENT(IN)    :: BCON(:, :) ! boundary concentrations
      REAL, INTENT(OUT)   :: XVEL(:, :) ! adjusted velocities along rows
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

      INTEGER :: C, R, S        ! loop counters

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

      XVEL(1:NCOLS+1, 1:NROWS) = UHAT(1:NCOLS+1, 1:NROWS)

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

C Insert the boundary conditions into CON
      DO R = 1, NROWS
         IF ( BNDY_LO_PE ) THEN ! west boundary
            DO S = 1, N_SPC_ADV
               IF ( XVEL(1, R) < 0.0 ) THEN ! outflow
                  CON(1-SWP:0, R, S) = ZFDBC(
     &                 CON(1, R, S), CON(2, R, S),
     &                 XVEL(1, R), XVEL(2, R) )
               ELSE             ! inflow
                  CON(1-SWP:0, R, S) = BCON(WFX+R, S)
               END IF
            END DO
         END IF

         IF ( BNDY_HI_PE ) THEN ! east boundary
            DO S = 1, N_SPC_ADV
               IF ( XVEL(NCOLS+1, R) > 0.0 ) THEN ! outflow
                  CON(NCOLS+1:NCOLS+SWP, R, S) = ZFDBC(
     &                 CON(NCOLS, R, S), CON(NCOLS-1, R, S),
     &                 XVEL(NCOLS+1, R), XVEL(NCOLS, R) )
               ELSE             ! inflow
                  CON(NCOLS+1:NCOLS+SWP, R, S) = BCON(EFX+R, S)
               END IF
            END DO
         END IF

C Do the gridded computation for horizontal advection using the PPM scheme
         CALL HPPM(NCOLS, N_SPC_ADV, CON(:, R, :), XVEL(:, R), DT, DS, 'C')
      END DO

      RETURN
      END
