C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      MODULE ADJ_ARRAY_MOD

C-----------------------------------------------------------------------
C Function:
C   Define adjoint forcing and sensitivity arrays

C Revision History:
C   Nov 2013 by Shannon Capps (US EPA)
C            based on the FOURDVAR_MOD created by Peter Percell (UH)
C            left significant portion commented for potential use in 4D-Var
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

C-----------------------------------
C Data
C-----------------------------------

CC Number of observed species
C      INTEGER :: N_OBS_SPC
C
CC Observed species names and units
C      CHARACTER(16), ALLOCATABLE :: OBS_SPC_NAME(:)
C      CHARACTER(16), ALLOCATABLE :: OBS_SPC_UNITS(:)
C
CC Map from the observed species list to the CGRID species list
C      INTEGER, ALLOCATABLE :: OBS_TO_CGRID_MAP(:)
C
CC Number of control species
C      INTEGER :: N_CNTL_SPC
C
CC Control species names and units
C      CHARACTER(16), ALLOCATABLE :: CNTL_SPC_NAME(:)
C      CHARACTER(16), ALLOCATABLE :: CNTL_SPC_UNITS(:)
C
CC Maps between the control species list and the CGRID, advected and diffused
CC species lists
C      INTEGER, ALLOCATABLE :: CNTL_TO_CGRID_MAP(:)
C      INTEGER, ALLOCATABLE :: CNTL_TO_ADV_MAP(:)
C      INTEGER, ALLOCATABLE :: ADV_TO_CNTL_MAP(:)
C      INTEGER, ALLOCATABLE :: CNTL_TO_DIFF_MAP(:)
C      INTEGER, ALLOCATABLE :: DIFF_TO_CNTL_MAP(:)
C
CC COST FUNCTION
C      REAL              :: COST_FUNC
C      INTEGER, ALLOCATABLE :: RECPTR_DEF(:,:)
C      LOGICAL           :: output_save = .false.
C
C! Mortality Data and CONC_AVG
C      REAL, ALLOCATABLE :: MORTALITY(:,:)
C      REAL, ALLOCATABLE :: CONC_AVG(:,:,:,:)
C      REAL, PARAMETER    :: BETA = 0.005827 ! total CRF of PM2.5
C
C! Variables for receptor species
C      INTEGER :: SPC, SPC_CONC
C      INTEGER :: CF_MYPE = -1
C      INTEGER :: BCOL, ECOL, BROW, EROW, BLEV, ELEV
C
CC First timestep for CF calculation
C      CHARACTER(16)     :: CF_BEGIN_DATE = "CF_BEGIN_DATE"
C      CHARACTER(16)     :: CF_BEGIN_TIME = "CF_BEGIN_TIME"
C      INTEGER           :: CF_STDATE
C      CHARACTER(16)     :: CF_STDATE_STR
C      INTEGER           :: CF_STTIME
C
CC Final timestep for CF calculation
C      CHARACTER(16)     :: CF_END_DATE = "CF_END_DATE"
C      CHARACTER(16)     :: CF_END_TIME = "CF_END_TIME"
C      INTEGER           :: CF_EDATE
C      CHARACTER(16)     :: CF_EDATE_STR
C      INTEGER           :: CF_ETIME
C
CC Background values of control variables
C      REAL, ALLOCATABLE :: BG_GRID(:, :, :, :)
C
CC Current values of control variables
C      INTEGER    :: CURRENT_TIME, ROW_SAVE, COL_SAVE, LAY_SAVE   !debug, mdt
C      INTEGER    :: ROW_SAVE1, COL_SAVE1, LAY_SAVE1   !debug, mdt
C
CC Current values of control variables
C      INTEGER    :: CURRENT_DATE
C
CC Current values of control variables
C      REAL, ALLOCATABLE :: CNTL_GRID(:, :, :, :)
C
CC Emissions scaling factor
C      REAL, ALLOCATABLE :: EM_SF(:, :, :, :)
C
CC debugt, mdt :: use coagulation?
C      LOGICAL :: USE_COAG = .FALSE.
C
CC debug, mdt :: use mode merging?
C      LOGICAL :: USE_MM = .TRUE.
C
C *** For now, keep this declaration outside of module & pass b/t subroutines
CC Gradient of cost function (unitless)
C      REAL, ALLOCATABLE :: LGRID(:, :, :, :)  ! adjoint accumulation variable

C Adjoint forcing function - at synchronization time step (unitless)
      REAL, ALLOCATABLE :: LGRID_FRC(:, :, :, :)   ! adjoint forcing at sync step

C Adjoint forcing function - at output time step (unitless)
      REAL, ALLOCATABLE :: LGRID_FRC_TOT(:, :, :, :)   ! adjoint forcing at output step from file

C Gradient of cost function wrt emissions
      REAL, ALLOCATABLE :: LGRID_EM(:, :, :, :)

C Gradient of cost function wrt emissions scaling factor
      REAL, ALLOCATABLE :: LGRID_EM_SF(:, :, :, :)

C Gradient of cost function wrt emissions - fully normalized
      REAL, ALLOCATABLE :: LGRID_EM_NRM(:, :, :, :)

CC Boundary condition scaling factor
C      REAL, ALLOCATABLE :: BC_SF(:, :, :)
C
C Gradient of cost function wrt boundary condition scaling factor
C      REAL, ALLOCATABLE :: LGRID_BC_SF(:, :, :)
C
CC Background values of control variables in boundary condition time series
C      REAL, ALLOCATABLE :: BG_BC_TS(:, :, :, :)
C
CC Current values of control variables in boundary condition time series
C      REAL, ALLOCATABLE :: CNTL_BC_TS(:, :, :, :)
C
CC Gradient of cost function wrt control variables in boundary condition time
CC series
C      REAL, ALLOCATABLE :: LGRID_BC_TS(:, :, :, :)
C
CC Map from advected species to BC species names
C      CHARACTER(16), ALLOCATABLE :: BCNAME_ADV(:)
C
CC Map from control species to BC species names
C      CHARACTER(16), ALLOCATABLE :: BCNAME_CNTL(:)

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      CONTAINS

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Initialize data
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FOURDVAR_INIT(JDATE, JTIME)

      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_GC_DIFF     ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF     ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF     ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF     ! tracer diffusion species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! starting date (YYYYDDD)
      INTEGER, INTENT(IN) :: JTIME ! starting time (HHMMSS)

C External Functions (not already declared by IODECL3.EXT):
      INTEGER, EXTERNAL :: INDEX1

C Local Variables:
      CHARACTER(16) :: PNAME = 'FOURDVAR_INIT'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(16) :: SPC_LIST(N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC)
      LOGICAL, ALLOCATABLE :: SPC_KNOWN(:)
      INTEGER :: N_UNKNOWN

      INTEGER :: ALLOCSTAT

      INTEGER :: S              ! species index
      INTEGER :: STRT           ! loop limit
      INTEGER :: INDX

      INTERFACE
         SUBROUTINE GET_ENVLIST(ENV_VAR, NVARS, V_LIST)
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN) :: ENV_VAR
            INTEGER, INTENT(OUT) :: NVARS
            CHARACTER(16), INTENT(OUT) :: V_LIST(:)
         END SUBROUTINE GET_ENVLIST
      END INTERFACE

C-----------------------------------------------------------------------

C Get CGRID offsets
      CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

C Get the list of observed species
      CALL GET_ENVLIST('OBS_SPCS', N_OBS_SPC, SPC_LIST)

      IF ( N_OBS_SPC == 1 .AND. TRIM(SPC_LIST(1)) == 'ALL' ) THEN
         N_OBS_SPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC ! everything

         STRT = 0
         DO INDX = 1, N_GC_SPC
            SPC_LIST(STRT + INDX) = GC_SPC(INDX)
         END DO

         STRT = N_GC_SPC
         DO INDX = 1, N_AE_SPC
            SPC_LIST(STRT + INDX) = AE_SPC(INDX)
         END DO

         STRT = N_GC_SPC + N_AE_SPC
         DO INDX = 1, N_NR_SPC
            SPC_LIST(STRT + INDX) = NR_SPC(INDX)
         END DO

         STRT = N_GC_SPC + N_AE_SPC + N_NR_SPC
         DO INDX = 1, N_TR_SPC
            SPC_LIST(STRT + INDX) = TR_SPC(INDX)
         END DO
      END IF

      IF ( N_OBS_SPC > 0 ) THEN

C Allocate memory for the map from the observed species list to the CGRID
C species list and for the lists of observed species names and units
         ALLOCATE(
     &        OBS_TO_CGRID_MAP(1:N_OBS_SPC),
     &        OBS_SPC_NAME(1:N_OBS_SPC),
     &        OBS_SPC_UNITS(1:N_OBS_SPC),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT /= 0 ) THEN
            XMSG = 'Failure allocating memory for observed species data'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF
         OBS_TO_CGRID_MAP(:) = 0

C Create map from the observed species list to the CGRID species list and
C create the lists of observed species names and units
         ALLOCATE(
     &        SPC_KNOWN(1:N_OBS_SPC),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT /= 0 ) THEN
            XMSG = 'Failure allocating memory for SPC_KNOWN'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF
         SPC_KNOWN = .FALSE.

         DO S = 1, N_OBS_SPC
            INDX = INDEX1(SPC_LIST(S), N_GC_SPC, GC_SPC)
            IF ( INDX > 0 ) THEN
               OBS_TO_CGRID_MAP(S) = GC_STRT + INDX - 1
               OBS_SPC_NAME(S) = SPC_LIST(S)
               OBS_SPC_UNITS(S) = 'ppmV'
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         DO S = 1, N_OBS_SPC
            INDX = INDEX1(SPC_LIST(S), N_AE_SPC, AE_SPC)
            IF ( INDX > 0 ) THEN
               OBS_TO_CGRID_MAP(S) = AE_STRT + INDX - 1
               OBS_SPC_NAME(S) = SPC_LIST(S)
               IF ( OBS_SPC_NAME(S)(1:3) == 'NUM' ) THEN
                  OBS_SPC_UNITS(S) = 'number/m**3'
               ELSE IF ( OBS_SPC_NAME(S)(1:3) == 'SRF' ) THEN
                  OBS_SPC_UNITS(S) = 'm**2/m**3'
               ELSE
                  OBS_SPC_UNITS(S) = 'microg/m**3'
               END IF
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         DO S = 1, N_OBS_SPC
            INDX = INDEX1(SPC_LIST(S), N_NR_SPC, NR_SPC)
            IF ( INDX > 0 ) THEN
               OBS_TO_CGRID_MAP(S) = NR_STRT + INDX - 1
               OBS_SPC_NAME(S) = SPC_LIST(S)
               OBS_SPC_UNITS(S) = 'ppmV'
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         DO S = 1, N_OBS_SPC
            INDX = INDEX1(SPC_LIST(S), N_TR_SPC, TR_SPC)
            IF ( INDX > 0 ) THEN
               OBS_TO_CGRID_MAP(S) = TR_STRT + INDX - 1
               OBS_SPC_NAME(S) = SPC_LIST(S)
               OBS_SPC_UNITS(S) = 'ppmV'
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         N_UNKNOWN = 0
         DO S = 1, N_OBS_SPC
            IF ( .NOT. SPC_KNOWN(S) ) THEN
               XMSG = 'Species ' // TRIM(SPC_LIST(S))
     &              // ' is unknown for this model'
               CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
               N_UNKNOWN = N_UNKNOWN + 1
            END IF
         END DO
         IF ( N_UNKNOWN > 0 ) CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)

         IF ( ALLOCATED(SPC_KNOWN) ) DEALLOCATE(SPC_KNOWN)

      END IF                    ! N_OBS_SPC > 0

C Get the list of control species
      CALL GET_ENVLIST('CNTL_SPCS', N_CNTL_SPC, SPC_LIST)

      IF ( N_CNTL_SPC == 1 .AND. TRIM(SPC_LIST(1)) == 'ALL' ) THEN
         N_CNTL_SPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC ! everything

         STRT = 0
         DO INDX = 1, N_GC_SPC
            SPC_LIST(STRT + INDX) = GC_SPC(INDX)
         END DO

         STRT = N_GC_SPC
         DO INDX = 1, N_AE_SPC
            SPC_LIST(STRT + INDX) = AE_SPC(INDX)
         END DO

         STRT = N_GC_SPC + N_AE_SPC
         DO INDX = 1, N_NR_SPC
            SPC_LIST(STRT + INDX) = NR_SPC(INDX)
         END DO

         STRT = N_GC_SPC + N_AE_SPC + N_NR_SPC
         DO INDX = 1, N_TR_SPC
            SPC_LIST(STRT + INDX) = TR_SPC(INDX)
         END DO
      END IF

      IF ( N_CNTL_SPC > 0 ) THEN

C Allocate memory for maps between the control species list and the CGRID,
C advected and diffused species lists and for the lists of control species
C names and units
         ALLOCATE(
     &        CNTL_TO_CGRID_MAP(1:N_CNTL_SPC),
     &        CNTL_TO_ADV_MAP(1:N_CNTL_SPC),
     &        ADV_TO_CNTL_MAP(N_GC_ADV+N_AE_ADV+N_NR_ADV+N_TR_ADV+1),
     &        CNTL_TO_DIFF_MAP(1:N_CNTL_SPC),
     &        DIFF_TO_CNTL_MAP(N_GC_DIFF+N_AE_DIFF+N_NR_DIFF+N_TR_DIFF),
     &        CNTL_SPC_NAME(1:N_CNTL_SPC),
     &        CNTL_SPC_UNITS(1:N_CNTL_SPC),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT /= 0 ) THEN
            XMSG = 'Failure allocating memory for control species data'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF
         CNTL_TO_CGRID_MAP(:) = 0
         CNTL_TO_ADV_MAP(:) = 0
         ADV_TO_CNTL_MAP(:) = 0
         CNTL_TO_DIFF_MAP(:) = 0
         DIFF_TO_CNTL_MAP(:) = 0

C Create map from the control species list to the CGRID species list and
C create the lists of control species names and units
         ALLOCATE(
     &        SPC_KNOWN(1:N_CNTL_SPC),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT /= 0 ) THEN
            XMSG = 'Failure allocating memory for SPC_KNOWN'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF
         SPC_KNOWN = .FALSE.

         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(SPC_LIST(S), N_GC_SPC, GC_SPC)
            IF ( INDX > 0 ) THEN
               CNTL_TO_CGRID_MAP(S) = GC_STRT + INDX - 1
               CNTL_SPC_NAME(S) = SPC_LIST(S)
               CNTL_SPC_UNITS(S) = 'ppmV'
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(SPC_LIST(S), N_AE_SPC, AE_SPC)
            IF ( INDX > 0 ) THEN
               CNTL_TO_CGRID_MAP(S) = AE_STRT + INDX - 1
               CNTL_SPC_NAME(S) = SPC_LIST(S)
               IF ( CNTL_SPC_NAME(S)(1:3) == 'NUM' ) THEN
                  CNTL_SPC_UNITS(S) = 'number/m**3'
               ELSE IF ( CNTL_SPC_NAME(S)(1:3) == 'SRF' ) THEN
                  CNTL_SPC_UNITS(S) = 'm**2/m**3'
               ELSE
                  CNTL_SPC_UNITS(S) = 'microg/m**3'
               END IF
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(SPC_LIST(S), N_NR_SPC, NR_SPC)
            IF ( INDX > 0 ) THEN
               CNTL_TO_CGRID_MAP(S) = NR_STRT + INDX - 1
               CNTL_SPC_NAME(S) = SPC_LIST(S)
               CNTL_SPC_UNITS(S) = 'ppmV'
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(SPC_LIST(S), N_TR_SPC, TR_SPC)
            IF ( INDX > 0 ) THEN
               CNTL_TO_CGRID_MAP(S) = TR_STRT + INDX - 1
               CNTL_SPC_NAME(S) = SPC_LIST(S)
               CNTL_SPC_UNITS(S) = 'ppmV'
               SPC_KNOWN(S) = .TRUE.
            END IF
         END DO

         N_UNKNOWN = 0
         DO S = 1, N_CNTL_SPC
            IF ( .NOT. SPC_KNOWN(S) ) THEN
               XMSG = 'Species ' // TRIM(SPC_LIST(S))
     &              // ' is unknown for this model'
               CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
               N_UNKNOWN = N_UNKNOWN + 1
            END IF
         END DO
         IF ( N_UNKNOWN > 0 ) CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)

         IF ( ALLOCATED(SPC_KNOWN) ) DEALLOCATE(SPC_KNOWN)

C Create maps between the control species list and the advected species list
         STRT = 0
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_GC_ADV, GC_ADV)
            IF ( INDX > 0 ) THEN
               CNTL_TO_ADV_MAP(S) = STRT + INDX
               ADV_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

         STRT = N_GC_ADV
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_AE_ADV, AE_ADV)
            IF ( INDX > 0 ) THEN
               CNTL_TO_ADV_MAP(S) = STRT + INDX
               ADV_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

         STRT = N_GC_ADV + N_AE_ADV
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_NR_ADV, NR_ADV)
            IF ( INDX > 0 ) THEN
               CNTL_TO_ADV_MAP(S) = STRT + INDX
               ADV_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

         STRT = N_GC_ADV + N_AE_ADV + N_NR_ADV
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_TR_ADV, TR_ADV)
            IF ( INDX > 0 ) THEN
               CNTL_TO_ADV_MAP(S) = STRT + INDX
               ADV_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

C Create maps between the control species list and the diffused species list
         STRT = 0
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_GC_DIFF, GC_DIFF)
            IF ( INDX > 0 ) THEN
               CNTL_TO_DIFF_MAP(S) = STRT + INDX
               DIFF_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

         STRT = N_GC_DIFF
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_AE_DIFF, AE_DIFF)
            IF ( INDX > 0 ) THEN
               CNTL_TO_DIFF_MAP(S) = STRT + INDX
               DIFF_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

         STRT = N_GC_DIFF + N_AE_DIFF
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_NR_DIFF, NR_DIFF)
            IF ( INDX > 0 ) THEN
               CNTL_TO_DIFF_MAP(S) = STRT + INDX
               DIFF_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

         STRT = N_GC_DIFF + N_AE_DIFF + N_NR_DIFF
         DO S = 1, N_CNTL_SPC
            INDX = INDEX1(CNTL_SPC_NAME(S), N_TR_DIFF, TR_DIFF)
            IF ( INDX > 0 ) THEN
               CNTL_TO_DIFF_MAP(S) = STRT + INDX
               DIFF_TO_CNTL_MAP(STRT + INDX) = S
            END IF
         END DO

      END IF                    ! N_CNTL_SPC > 0

      END SUBROUTINE FOURDVAR_INIT

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      END MODULE FOURDVAR_MOD