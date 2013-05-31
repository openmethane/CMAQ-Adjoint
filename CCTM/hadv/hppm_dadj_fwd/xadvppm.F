
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/hppm_dadj_fwd/xadvppm.F,v 1.4 2013/05/14 21:06:52 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE XADV(CGRID, JDATE, JTIME, TSTEP, NADVSS)

C-----------------------------------------------------------------------
C Function:
C   Advection in the horizontal plane; x1-direction:
C   The process time step is set equal to TSTEP(2). Boundary concentrations
C   are coupled in RDBCON with SqRDMT = Sq. Root [det ( metric tensor )]
C   = Jacobian / (map scale factor)**2
C   where Air Density X SqRDMT is loaded into last BCON slot for advection.

C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times.  All boundary input variables have the
C   same boundary perimeter structure with a thickness of 1
C   CGRID in transport units: SQRT{DET[metric tensor]}*concentration (Mass/Vol)

C Revision history:
C   01/20/09 by Peter Percell at UH-IMAQS: rewrote for adjoint/4dvar
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
      REAL, POINTER :: CGRID(:, :, :, :)
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
      REAL :: UHAT(1:NCOLS+1, 1:NROWS+1) ! x1-component CX-velocity
      REAL :: BCON(1:NBNDY, 1:N_SPC_ADV) ! boundary concentrations

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: UORV = 'X1VEL' ! flag for HCONTVEL

      INTEGER, SAVE :: EFX      ! fixed parameter for eastern boundary
      INTEGER, SAVE :: WFX      ! fixed parameter for western boundary

      INTEGER, SAVE :: ADV_MAP(1:N_SPC_ADV) ! global adv map to CGRID

      LOGICAL, SAVE :: BNDY_LO_PE, BNDY_HI_PE

      INTEGER :: SYNC           ! sync step (sec)
      INTEGER :: ASTEP          ! advection time step for layer (HHMMSS)
      INTEGER :: SS             ! loop counter

      REAL :: DT                ! advection time step for layer (sec)
      INTEGER :: DT_INT         ! advection time step for layer (sec)

      INTEGER :: ADATE          ! interpolation date
      INTEGER :: ATIME          ! interpolation time

      INTEGER :: C, R, L, S, V  ! loop counters
      INTEGER :: A2C

      REAL, TARGET :: CON(1-SWP:NCOLS+SWP, 1:NROWS, 1:N_SPC_ADV)
#ifdef parallel
      REAL, POINTER :: SEND(:, :, :)
      REAL, POINTER :: RECV(:, :, :)
#endif

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

         EFX = NCOLS + 1
         WFX = 2 * NCOLS + NROWS + 4

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

         CALL SUBST_HI_LO_BND_PE('C', BNDY_LO_PE, BNDY_HI_PE)

C Create new adjoint file and open for update
         CALL CHK_FILE_OPEN ( ADJ_XADV_CHK, JDATE, JTIME, .TRUE., TSTEP )

      END IF                    ! FIRSTIME

C Write horizontal advection checkpoint at sync step
      CALL CHK_FILE_WRITE ( ADJ_XADV_CHK, JDATE, JTIME,
     &     CGRID(1:NCOLS, 1:NROWS, 1:NLAYS, ADV_MAP(1:N_SPC_ADV)) )

      SYNC = TIME2SEC(TSTEP(2))
      DO L = 1, NLAYS

C Load the interior concentrations into CON
         DO S = 1, N_SPC_ADV
            A2C = ADV_MAP(S)
            CON(1:NCOLS, 1:NROWS, S) = CGRID(1:NCOLS, 1:NROWS, L, A2C)
         END DO

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

C Get the contravariant x1 velocity component
            CALL HCONTVEL(ADATE, ATIME, ASTEP, L, UORV, UHAT)

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

C Get the boundary concentration values (converted to transport units)
            CALL RDBCON(ADATE, ATIME, ASTEP, L, BCON)

C Insert the boundary conditions into CON
            IF ( BNDY_LO_PE ) THEN ! west boundary
               DO S = 1, N_SPC_ADV
                  DO R = 1, NROWS
                     IF ( UHAT(1, R) < 0.0 ) THEN ! outflow
                        CON(1-SWP:0, R, S) = ZFDBC(
     &                       CON(1, R, S), CON(2, R, S),
     &                       UHAT(1, R), UHAT(2, R) )
                     ELSE       ! inflow
                        CON(1-SWP:0, R, S) = BCON(WFX+R, S)
                     END IF
                  END DO
               END DO
            END IF

            IF ( BNDY_HI_PE ) THEN ! east boundary
               DO S = 1, N_SPC_ADV
                  DO R = 1, NROWS
                     IF ( UHAT(NCOLS+1, R) > 0.0 ) THEN ! outflow
                        CON(NCOLS+1:NCOLS+SWP, R, S) = ZFDBC(
     &                       CON(NCOLS, R, S), CON(NCOLS-1, R, S),
     &                       UHAT(NCOLS+1, R), UHAT(NCOLS, R) )
                     ELSE       ! inflow
                        CON(NCOLS+1:NCOLS+SWP, R, S) = BCON(EFX+R, S)
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
     &           SWP, NROWS, EAST)

            SEND => CON(NCOLS-SWP+1:NCOLS, 1:NROWS, 1:N_SPC_ADV)
            RECV => CON(1-SWP:0, 1:NROWS, 1:N_SPC_ADV)
            CALL SWAP3D(SEND, RECV, SWP, NROWS, N_SPC_ADV,
     &           SWP, NROWS, WEST)
#endif

C Do the gridded computation for horizontal advection using the PPM scheme
            DO R = 1, NROWS     ! x1-direction
               CALL HPPM(NCOLS, N_SPC_ADV, CON(:, R, :), UHAT(:, R),
     &              DT, REAL(XCELL_GD), 'C')
            END DO

            CALL NEXTIME(ADATE, ATIME, ASTEP)

         END DO                 ! SS

C Transfer updated interior concentrations from CON into CGRID
         DO S = 1, N_SPC_ADV
            A2C = ADV_MAP(S)
            CGRID(1:NCOLS, 1:NROWS, L, A2C) = CON(1:NCOLS, 1:NROWS, S)
         END DO

      END DO                    ! NLAYS

      RETURN
      END
