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
C Emissions scaling factor
      REAL, ALLOCATABLE :: EM_SF(:, :, :, :)
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

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      END MODULE ADJ_ARRAY_MOD
