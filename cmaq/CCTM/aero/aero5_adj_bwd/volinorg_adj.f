
!***********************************************************************
!   Portions of Models-3/CMAQ software were developed or based on      *
!   information from various groups: Federal Government employees,     *
!   contractors working on a United States Government contract, and    *
!   non-Federal sources (including research institutions).  These      *
!   research institutions have given the Government permission to      *
!   use, prepare derivative works, and distribute copies of their      *
!   work in Models-3/CMAQ to the public and to permit others to do     *
!   so.  EPA therefore grants similar permissions for use of the       *
!   Models-3/CMAQ software, but users are requested to provide copies  *
!   of derivative works to the Government without restrictions as to   *
!   use by others.  Users are responsible for acquiring their own      *
!   copies of commercial software associated with Models-3/CMAQ and    *
!   for complying with vendor requirements.  Software copyrights by    *
!   the MCNC Environmental Modeling Center are used with their         *
!   permissions subject to the above restrictions.                     *
!***********************************************************************

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VOLINORG_ADJ( DT, COL, ROW, LAYER, DV_SO4, CGR, 
     &                         M3_WET_FLAG )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of subroutine that calculates the partitioning of inorganic
!   components (CL, NO3, NH4, SO4) between the aerosol and gas phase
!   over the operator synchronization timestep (DT).  

! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      USE SOA_DEFN_B
      USE MET_DATA
      USE AERO_DATA_B
      USE PRECURSOR_DATA_B

      IMPLICIT NONE

! *** Arguments:

      REAL    DT              ! time step [sec]
      INTEGER COL             ! grid column index
      INTEGER ROW             ! grid row index
      INTEGER LAYER           ! model layer index
      REAL    DV_SO4          ! molecular diffusivity of H2SO4 vapor 
                              ! after correction for ambient conditions
      REAL( 8 ) :: CGR( N_MODE-1 ) ! 3rd moment SO4 growth rate [m^3/m^3-s]
      LOGICAL, INTENT ( IN ) ::  M3_WET_FLAG  ! flag to include water in GETPAR update

! *** Parameters: 

      INTEGER, PARAMETER :: NINORG = 6          ! number of inorganic species
      INTEGER, PARAMETER :: NVOLINORG = 3       ! number of volatile inorganic species

      ! indices for inorganic species
      INTEGER, PARAMETER :: KNH4 = 1, KNO3 = 2, KCL  = 3, KSO4 = 4, KNA = 5, KHP = 6

      REAL( 8 ), PARAMETER :: MWH2SO4 = 98.0 ! molecular weight for H2SO4
      REAL( 8 ), PARAMETER :: FAERH2SO4 = 1.0E-6 / MWH2SO4
      REAL( 8 ), PARAMETER :: D_TWOTHIRDS = 2.0D0 / 3.0D0

      REAL, PARAMETER :: CUTOFF = 0.05   ! [ug/m**3]

      REAL, PARAMETER :: ALPHSULF = 0.1 ! Accommodation coefficient for sulfuric acid
                                        ! see Capaldo et al. (2000)

! *** Local Variables:

! *** Inputs to subroutine HCOND3

      REAL, SAVE :: COFCBAR_SO4  ! Temperature-independent coefficients
                                 ! for caculating molecular vel [m/s]
                                 ! = sqrt((8*Rgas)/(pi*MW)) 
      REAL         CBAR_SO4      ! molecular velocity of H2SO4                      

      REAL( 8 ) :: AM0( N_MODE ) ! zeroth moments
      REAL( 8 ) :: AM1( N_MODE ) ! first moments
      REAL( 8 ) :: AM2( N_MODE ) ! second moments

! *** Outputs from HCOND3: size-dependent term in the condensational-growth 
!     expressions defined in Equations A13-A14 of [Binkowski & Shankar,1995]
      REAL( 8 ) :: FCONC_SO4( N_MODE,2 )  ! All sizes 2nd and 3rd moments
      REAL( 8 ) :: FCONCM1_SO4       ! reciprocals of total SO4 cond rates

! *** Modal partition factors [ dimensionless ]
!     defined in Equations A17-A18 of [Binkowski & Shankar,1995]
      REAL( 8 ) :: OMEGA_AT_SO4  ! Aitken mode 2nd and 3rd moments
      REAL( 8 ) :: OMEGA_AC_SO4  ! Accumulation mode 2nd and 3rd moments
      REAL( 8 ) :: OMEGA( 2 )    ! partitioning coefficient for equilibrium PM mass

! *** Variables for new particle formation:
      REAL XH2SO4            ! steady state H2SO4 concentration
      REAL( 8 ) :: DMDT_SO4  ! particle mass production rate [ ug/m**3 s ]
      REAL( 8 ) :: DNDT      ! particle number production rate [ # / m**3 s ]
      REAL( 8 ) :: DM2DT     ! second moment production rate [ m**2 / m**3 s]
      REAL( 8 ) :: SCONDRATE ! SO4 condensation rate [ ug/m**3 s ]

! *** Mode-specific sulfate production rate [ ug/m**3 s ]
      REAL(8) :: CONDSO4( N_MODE )      ! sulfate condensation rate [ ug/m**3 s ]
      REAL( 8 ) :: RATE                 ! CONDSO4 or cond+nucl rate

! *** Size-dependent portion of mass-transfer rate equation
      REAL( 8 ) :: GRFAC1( N_MODE )     ! 2nd moment [ m**2/m**3-s ] 
      REAL( 8 ) :: GRFAC2( N_MODE )     ! 3rd moment [ m**3/m**3-s ] 

! *** ISORROPIA input variables
      REAL( 8 ) :: WI( NINORG - 1 )     ! species array
      REAL( 8 ) :: RHI                  ! relative humidity
      REAL( 8 ) :: TEMPI                ! temperature
      REAL( 8 ) :: CNTRL( 2 )           ! control parameters 

! *** ISORROPIA output variables
      REAL( 8 ) :: WT( NINORG - 1 )     ! species output array
      REAL( 8 ) :: GAS( 3 )             ! gas-phase   "     " 
      REAL( 8 ) :: AERLIQ( 12 )         ! liq aerosol "     " 
      REAL( 8 ) :: AERSLD( 9 )          ! solid "     "     " 
      REAL( 8 ) :: OTHER( 6 )           ! supplmentary output array
      CHARACTER( 15 ) :: SCASI          ! subcase number output

! *** Variables to account for mass conservation violations in ISRP3F
      LOGICAL TRUSTNH4                  ! false if ISOROPIA's partitioning
                                        !  of NH4/NH3 is to be ignored
      LOGICAL TRUSTCL                   ! false if ISOROPIA's partitioning       
                                        !  of Cl/HCl is to be ignored

! *** Initial (double-precision) concentrations [ug/m3]
      REAL( 8 ) :: GNH3R8               ! gas-phase ammonia
      REAL( 8 ) :: GNO3R8               ! gas-phase nitric acid
      REAL( 8 ) :: GCLR8                ! gas-phase hydrochloric acid
!      REAL( 8 ) :: H2OR8                ! aerosol LWC

! *** Variables for volatile species mass transfer between gas and aerosol and
!     mass partitioning between the modes 
      LOGICAL HYBRID ! mass transfer option flag (mass transfer if .TRUE.)
      REAL( 8 ) :: DELT                 ! time step DT [s]
      REAL( 8 ) :: HPLUS( N_MODE )      ! scratch var for H+ [umol/m**3]

      REAL(8), SAVE :: H2SO4RATM1       ! Mol. wt. ratio of SO4/H2SO4

!      REAL( 8 ) :: JH2SO4 ! flux of cond./evap. H2SO4 (mol/m3/s)
!      REAL( 8 ) :: CH2SO4 ! concentration of H2SO4 before cond/evap (mol/m3)
      REAL( 8 ) :: DVOLINORG( NVOLINORG ) ! vol inorg spcs mass to be xferred [mol/m3]
      REAL( 8 ) :: DVOLMAX                ! max value for DVOLINORG 
      REAL( 8 ) :: CINORG( NINORG,N_MODE ) ! scratch array for inorg spcs [ug/m**3]
      REAL( 8 ) :: INT_TIME               ! internal mass transfer time (s)
      REAL( 8 ) :: TSTEP                  ! mass transfer time step [s]
      REAL( 8 ) :: DRYM20, Y              ! scratch vars for 2nd moment [m**2/m**3]
      REAL( 8 ) :: M3OTHR, M3SOA          ! vars for 3rd moment calculation [m**3/m**3]
      REAL( 8 ), SAVE :: DF( NVOLINORG )  ! scratch array for mole -> ug conversion
      REAL( 8 ), SAVE :: DFH2OR8          ! mole -> ug conversion for H2O
      REAL( 8 ) :: CONVFAC           ! conversion factor for SO4 conc from DH2SO4/DT
      REAL( 8 ) :: J(NVOLINORG)      ! condensation/evaporation flux [mol/m**3-s]
      REAL( 8 ) :: CFINAL( NVOLINORG,N_MODE ) ! conc after mass xfer step [ug/m**3]
      REAL( 8 ) :: H2O                  ! Scratch LWC variable for output
      REAL( 8 ) :: H2O_NEW              ! Update of LWC after condensation 
      REAL( 8 ) :: SO4                  ! modal SO4 after condensation or cond + nucl
      REAL( 8 ) :: DMASS( NVOLINORG )   ! excess mass in fine mode partitioning
      REAL( 8 ) :: DDRYM3DT             ! rate of 3rd moment transfer - dry inorg spcs
      REAL( 8 ) :: DDRYM2DT             ! rate of 2nd moment transfer -  "     "    "
      REAL( 8 ) :: DRYM3, WETM3         ! scratch vars for 3rd moment calc [m**3/m**3]
      REAL( 8 ) :: DRYM2, WETM2         ! scratch vars for 2nd moment calc [m**2/m**3]
      REAL( 8 ) :: LOSS                 ! rate of loss of second moment [1/s] 

      REAL( 8 ) :: TMP
      INTEGER I, N, B                   ! loop and array indices
      INTEGER IMODE                     ! mode loop index  
      INTEGER ISTEP                     ! loop index, mass transfer time step loop
      INTEGER ISP                       ! loop index, species loop
!      INTEGER LOGDEV
      LOGICAL, PARAMETER :: LIMIT_Sg = .TRUE.   ! fix Sg at current value?
      LOGICAL TrustIso                  ! For negative vap. press., TrustIso = F 

! *** Local Saved Variables 
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL, SAVE :: FIRSTWRITE = .TRUE.
      REAL( 8 ), SAVE :: SO4FAC         !  F6DPIM9 / RHOSO4
      REAL( 8 ), SAVE :: SOILFAC        !  F6DPIM9 / RHOSOIL
      REAL( 8 ), SAVE :: ANTHFAC        !  F6DPIM9 / RHOANTH
      REAL( 8 ), SAVE :: NH3RAT         ! Mol. wt ratios NH3/NH4
      REAL( 8 ), SAVE :: HNO3RAT        ! Mol. wt ratios HNO3/NO3
      REAL( 8 ), SAVE :: HCLRAT         ! Mol. wt ratios HCL/CL

! *** adjoint variables
      REAL( 8 ) :: CGR_ADJ ( N_MODE - 1 )
      REAL( 8 ) :: AM0_ADJ_0
      REAL( 8 ) :: WI_0( 5 )
      REAL( 8 ) :: AM1_ADJ_0
      REAL( 8 ) :: AM2_ADJ_0
      REAL( 8 ) :: AM0_ADJ ( N_MODE )
      REAL( 8 ) :: AM1_ADJ ( N_MODE )
      REAL( 8 ) :: AM2_ADJ ( N_MODE )
      REAL( 8 ) :: FCONC_SO4_ADJ ( N_MODE, 2 )
      REAL( 8 ) :: FCONCM1_SO4_ADJ
      REAL( 8 ) :: OMEGA_AT_SO4_ADJ
      REAL( 8 ) :: OMEGA_AC_SO4_ADJ
      REAL( 8 ) :: OMEGA_ADJ ( 2 )
      REAL      :: XH2SO4_ADJ
      REAL( 8 ) :: DMDT_SO4_ADJ
      REAL( 8 ) :: DNDT_ADJ
      REAL( 8 ) :: DM2DT_ADJ
      REAL( 8 ) :: SCONDRATE_ADJ
      REAL( 8 ) :: CONDSO4_ADJ ( N_MODE )
      REAL( 8 ) :: RATE_ADJ
      REAL( 8 ) :: GRFAC1_ADJ ( N_MODE )
      REAL( 8 ) :: GRFAC2_ADJ ( N_MODE )
      REAL( 8 ) :: WI_ADJ ( NINORG - 1 )
      REAL( 8 ) :: WI_ADJ_SAVE ( NINORG - 1 )  
      REAL( 8 ) :: GAS_ADJ ( 3 )
      REAL( 8 ) :: AERLIQ_ADJ ( 12 )
      REAL( 8 ) :: GNH3R8_ADJ
      REAL( 8 ) :: GNO3R8_ADJ
      REAL( 8 ) :: GCLR8_ADJ
      REAL( 8 ) :: DVOLINORG_ADJ ( NVOLINORG )
      REAL( 8 ) :: DVOLMAX_ADJ
      REAL( 8 ) :: HPLUS_ADJ ( N_MODE )
      REAL( 8 ) :: CINORG_ADJ ( NINORG, N_MODE )
      REAL( 8 ) :: DRYM20_ADJ
      REAL( 8 ) :: Y_ADJ
      REAL( 8 ) :: M3OTHR_ADJ
      REAL( 8 ) :: M3SOA_ADJ
      REAL( 8 ) :: J_ADJ ( NVOLINORG )
      REAL( 8 ) :: CFINAL_ADJ ( NVOLINORG, N_MODE )
      REAL( 8 ) :: H2O_ADJ
      REAL( 8 ) :: H2O_NEW_ADJ
      REAL( 8 ) :: SO4_ADJ
      REAL( 8 ) :: DMASS_ADJ ( NVOLINORG )
      REAL( 8 ) :: DDRYM3DT_ADJ
      REAL( 8 ) :: DDRYM2DT_ADJ 
      REAL( 8 ) :: DRYM3_ADJ
      REAL( 8 ) :: WETM3_ADJ
      REAL( 8 ) :: DRYM2_ADJ
      REAL( 8 ) :: WETM2_ADJ
      REAL( 8 ) :: LOSS_ADJ
      REAL( 8 ) :: TMP_ADJ
      REAL      :: SO4RATE_ADJ
      INTEGER   :: CALC_H2O_STEP
      INTEGER   :: COUNTER
      REAL( 8 ) :: TEMP4, TEMP4B5, TEMP3, TEMP3B9, TEMP2, TEMP2B5
      REAL( 8 ) :: TEMP1, TEMP1_ADJ

! Checkpointing variables
      REAL( 8 ) :: WETM3_CHECK( 2 ), DRYM2_CHECK( 2 ), DRYM3_CHECK( 2 )
      REAL( 8 ) :: WETM3_CHECK1( 2 ), WETM2_CHECK1( 2 ), DRYM3_CHECK1( 2 )
      REAL( 8 ) :: WETM3_CHECK2( 10 ), DRYM2_CHECK2( 10 )
!      REAL( 8 ), ALLOCATABLE :: WETM3_CHECK2( : ), DRYM2_CHECK2( : )
!      REAL( 8 ), ALLOCATABLE :: DRYM3_CHECK2( : ), Y_CHECK1( : )
      REAL( 8 ) :: DRYM3_CHECK2( 10 ), Y_CHECK1( 10 )
      REAL( 8 ) :: Y_CHECK( 2 ), WI_CHECK( 5 ), GAS_CHECK( 3 )
      REAL( 8 ) :: DDRYM2DT_CHECK( 10 ), J_CHECK2( 10, 3 )
      REAL( 8 ) :: DRYM20_CHECK( 10 ), LOSS_CHECK( 10 )
      REAL( 8 ) :: DDRYM3DT_CHECK( 10 ), GRFAC1_CHECK( 10, 3 )
      REAL( 8 ) :: GRFAC2_CHECK( 10, 3 ), WI_CHECK1( 10, 5 )
      REAL( 8 ) :: H2O_NEW_CHECK( 10 ), CINORG_CHECK( 10, 6, 3 )
      REAL( 8 ) :: J_CHECK( 10, 3 ), GNH3R8_CHECK( 10 ), GNO3R8_CHECK( 10 )
      REAL( 8 ) :: GCLR8_CHECK( 10 ), GAS_CHECK1( 10, 3 )
      REAL( 8 ) :: GRFAC2_CHECK1( 10, 3 ), HPLUS_CHECK( 10, 3 )
      REAL( 8 ) :: RATE_CHECK( 10 ), J_CHECK1( 10, 3 )
      REAL( 8 ) :: WI_CHECK2( 10, 5 )
!      REAL( 8 ), ALLOCATABLE :: DDRYM2DT_CHECK( : ), J_CHECK2( :, : )
!      REAL( 8 ), ALLOCATABLE :: DRYM20_CHECK( : ), LOSS_CHECK( : )
!      REAL( 8 ), ALLOCATABLE :: DDRYM3DT_CHECK( : ), GRFAC1_CHECK( :, : )
!      REAL( 8 ), ALLOCATABLE :: GRFAC2_CHECK( :, : ), WI_CHECK1( :, : )
!      REAL( 8 ), ALLOCATABLE :: H2O_NEW_CHECK( : ), CINORG_CHECK( :, :, : )
!      REAL( 8 ), ALLOCATABLE :: J_CHECK( :, : ), GNH3R8_CHECK( : ), GNO3R8_CHECK( : )
!      REAL( 8 ), ALLOCATABLE :: GCLR8_CHECK( : ), GAS_CHECK1( :, : )
!      REAL( 8 ), ALLOCATABLE :: GRFAC2_CHECK1( :, : ), HPLUS_CHECK( :, : )
!      REAL( 8 ), ALLOCATABLE :: RATE_CHECK( : ), J_CHECK1( :, : )
!      REAL( 8 ), ALLOCATABLE :: WI_CHECK2( :, : ) 
      REAL( 8 ) :: WT_CHECK( 5 )
      REAL( 8 ) :: AERLIQ_CHECK( 12 ), AERSLD_CHECK( 9 )
      REAL( 8 ) :: OTHER_CHECK( 6 )
!      REAL( 8 ), ALLOCATABLE :: WT_CHECK1( :, : )
      REAL( 8 ) :: WT_CHECK1( 10, 5 )
      REAL( 8 ) :: CNTRL_CHECK1( 10, 2 ), GAS_CHECK2( 10, 3 )
      REAL( 8 ) :: AERLIQ_CHECK1( 10, 12 ), AERSLD_CHECK1( 10, 9 )
      REAL( 8 ) :: OTHER_CHECK1( 10, 6 )
      REAL( 8 ) :: WETM3_CHECK3( 10 ), DRYM3_CHECK3( 10 )
      REAL( 8 ) :: WETM2_CHECK3( 10 ), AM0_CHECK( 10 ), AM1_CHECK( 10 )
      REAL( 8 ) :: AM2_CHECK( 10 ), FCONC_SO4_CHECK( 10, 2 )
!      REAL( 8 ), ALLOCATABLE :: CNTRL_CHECK1( :, : ), GAS_CHECK2( :, : )
!      REAL( 8 ), ALLOCATABLE :: AERLIQ_CHECK1( :, : ), AERSLD_CHECK1( :, : )
!      REAL( 8 ), ALLOCATABLE :: OTHER_CHECK1( :, : )
!      REAL( 8 ), ALLOCATABLE :: WETM3_CHECK3( : ), DRYM3_CHECK3( : )
!      REAL( 8 ), ALLOCATABLE :: WETM2_CHECK3( : ), AM0_CHECK( : ), AM1_CHECK( : )
!      REAL( 8 ), ALLOCATABLE :: AM2_CHECK( : ), FCONC_SO4_CHECK( :, : )
      REAL( 8 ) :: AM0_CHECK1( 3 ), AM1_CHECK1( 3 ), AM2_CHECK1( 3 )
      REAL( 8 ) :: FCONC_SO4_CHECK1( 3, 2 ), LOSS_CHECK1( 2  )
      REAL( 8 ) :: DDRYM2DT_CHECK1( 2 ), DRYM20_CHECK1( 2 )
      REAL( 8 ) :: DDRYM3DT_CHECK1( 2 ), WI_CHECK4( 2, 5 )
      REAL( 8 ) :: CFINAL_CHECK( 3, 2 ), CFINAL_CHECK1( 3, 2 )
      REAL( 8 ) :: DVOLINORG_CHECK( NVOLINORG ) 
      REAL( 8 ) :: DNDT_CHECK, DMDT_SO4_CHECK, DM2DT_CHECK
      REAL      :: aeromode_sdev_check( 10, 3 ), aeromode_mass_check( 10, 3 )
      REAL      :: aeromode_diam_check( 10, 3 ), aerospc_conc_check( 10, N_AEROSPC, N_MODE )
      REAL      :: moment2_conc_check( 10, 3 ), moment3_conc_check( 10, 3 )
      REAL      :: moment0_conc_check( 10, 3 )
      REAL      :: AEROMODE_MASS_CHK( N_MODE ), AEROMODE_SDEV_CHK( N_MODE )
      REAL      :: AEROMODE_DIAM_CHK( N_MODE ), AEROSPC_CONC_CHK( N_AEROSPC, N_MODE )
      REAL      :: MOMENT3_CONC_CHK( N_MODE ), MOMENT2_CONC_CHK( N_MODE )
      REAL      :: MOMENT0_CONC_CHK( N_MODE )
!      REAL( 8 ), ALLOCATABLE :: aeromode_sdev_check( :, : ), aeromode_mass_check( :, : )
!      REAL( 8 ), ALLOCATABLE :: aeromode_diam_check( :, : )
!      REAL( 8 ), ALLOCATABLE :: moment2_conc_check( :, : ), moment3_conc_check( :, : )
      LOGICAL   :: WET_CHECK( 2, N_AEROSPC ), TRUSTCL_CHECK

! FD Test variables 
      INTEGER :: N_P, X, P_I, X_I
      REAL*8 :: IN_P( 5 ), IN_P_ADJ( 5 ), OUT_P, OUT_P_ADJ
      REAL*8 :: AM0_0, AM1_0, AM2_0, FCONC_SO4_P( 2 )
      REAL   :: XH2SO4_ADJ_P, SO4RATE_ADJ_P
      REAL*8 :: GNH3R8_0, DNDT_0
      REAL*8 :: GNH3R8_ADJ_P, GNO3R8_ADJ_P, GCLR8_ADJ_P
      REAL   :: XH2SO4_0, SO4RATE_0
      REAL*8 :: GNO3R8_0, DMDT_SO4_0, DM2DT_0
      REAL*8 :: GCLR8_0, AM0_P, AM1_P, AM2_P
      REAL*8 :: GAS_0( 3 ), WI_TMP( 5 )
      REAL*8 :: GNH3R8_TMP, GNO3R8_TMP, GCLR8_TMP
      REAL*8 :: GRFAC2_TMP, RATE_TMP, AERLIQ_TMP1
      REAL*8 :: GRFAC2_0, GAS_ADJ_P( 3 ), GRFAC2_ADJ_P, AERLIQ_ADJ_P
      REAL*8 :: AERLIQ_0, RATE_ADJ_P, J_ADJ_P( 3 )
      REAL*8 :: RATE_0, DNDT_ADJ_P, DMDT_SO4_ADJ_P, DM2DT_ADJ_P
      REAL*8 :: J_0( 3 )
      REAL*8 :: J_TMP( 3 )
      CHARACTER( 80 ) :: filename, jdate_str, jtime_str
      CHARACTER( 80 ) :: filename_vol(30,3)
      CHARACTER( 80 ) :: filename1, filename2, i_str, counter_str
      CHARACTER( 80 ) :: filename3, filename4, filename5, filename6
      CHARACTER( 80 ) :: filename7, filename8, filename9, filename10
      CHARACTER( 80 ) :: filename11, filename12, filename13, filename14
      CHARACTER( 80 ) :: filename15, filename16, filename17, filename18
      CHARACTER( 80 ) :: filename19, filename20, filename21
      CHARACTER( 80 ) :: var_str, mode_str
      LOGICAL :: EXISTS
      REAL ::   moment0_conc_0( 3 )
      REAL ::   moment3_conc_0( 3 )
      REAL ::   moment2_conc_0( 3 )
      REAL ::   aerospc_conc_0( 30,3 )
      REAL ::   aeromode_diam_0( 3 ) 
      REAL ::   aeromode_dens_0( 3 )
      REAL ::   aeromode_sdev_0( 3 )
      REAL ::   aeromode_mass_0( 3 )
      REAL ::   D01, D02, D03, D04, D05, D06, D07, D08, D09, D10
      REAL ::   D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21
      REAL ::   D_VOL(30,3)

      REAL*8 :: WI_TMP1( NINORG - 1 )
      REAL*8 :: RHI_TMP
      REAL*8 :: WI_ADJ_CHK(5)
      REAL*8 :: TEMPI_TMP, WT_TMP( NINORG - 1 ), GAS_TMP( 3 )
      REAL*8 :: AERLIQ_TMP( 12 ), AERSLD_TMP( 9 ), OTHER_TMP( 6 )
      REAL*8 :: WI_ADJ_TMP( NINORG - 1 ), GAS_ADJ_TMP( 3 )
      REAL*8 :: AERLIQ_ADJ_TMP( 12 )
      CHARACTER( 15 ) :: SCASI_TMP
      LOGICAL :: TRUSTISO_TMP

!-----------------------------------------------------------------------

! adjoint variable declaration:

      CGR_ADJ = 0.0
      AM0_ADJ = 0.0
      AM1_ADJ = 0.0
      AM2_ADJ = 0.0
      FCONC_SO4_ADJ  = 0.0
      FCONCM1_SO4_ADJ = 0.0
      OMEGA_AT_SO4_ADJ = 0.0
      OMEGA_AC_SO4_ADJ = 0.0
      OMEGA_ADJ = 0.0
      XH2SO4_ADJ = 0.0
      DMDT_SO4_ADJ = 0.0
      DNDT_ADJ = 0.0
      DM2DT_ADJ = 0.0
      SCONDRATE_ADJ  = 0.0
      CONDSO4_ADJ = 0.0
      RATE_ADJ = 0.0
      GRFAC1_ADJ = 0.0
      GRFAC2_ADJ = 0.0
      WI_ADJ = 0.0
      GAS_ADJ = 0.0
      AERLIQ_ADJ = 0.0
      GNH3R8_ADJ = 0.0
      GNO3R8_ADJ = 0.0
      GCLR8_ADJ = 0.0
      DVOLINORG_ADJ = 0.0
      DVOLMAX_ADJ  = 0.0
      HPLUS_ADJ = 0.0
      CINORG_ADJ = 0.0
      DRYM20_ADJ = 0.0
      Y_ADJ = 0.0
      M3OTHR_ADJ = 0.0
      M3SOA_ADJ = 0.0
      J_ADJ = 0.0
      CFINAL_ADJ = 0.0
      H2O_ADJ = 0.0
      H2O_NEW_ADJ = 0.0
      SO4_ADJ = 0.0
      DMASS_ADJ = 0.0
      DDRYM3DT_ADJ = 0.0
      DDRYM2DT_ADJ = 0.0
      DRYM3_ADJ = 0.0
      WETM3_ADJ = 0.0
      DRYM2_ADJ = 0.0
      WETM2_ADJ = 0.0
      LOSS_ADJ = 0.0
      TMP_ADJ = 0.0
      CALC_H2O_STEP = 0
      COUNTER = 0
      SO4RATE_ADJ = 0.0

!  Code:

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         COFCBAR_SO4 = SQRT( 8.0 * RGASUNIV / ( PI * MWH2SO4 * 1.0E-3 ) )
         H2SO4RATM1 = AEROSPC_MW( ASO4_IDX ) / MWH2SO4
         SO4FAC  = 1.0E-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY
         SOILFAC = 1.0E-9 * F6DPI / AEROSPC( ASOIL_IDX )%DENSITY
         ANTHFAC = 1.0E-9 * F6DPI / AEROSPC( ACORS_IDX )%DENSITY

         DF( KNH4 ) = 1.0D6 * AEROSPC_MW( ANH4_IDX )
         DF( KNO3 ) = 1.0D6 * AEROSPC_MW( ANO3_IDX )
         DF( KCL )  = 1.0D6 * AEROSPC_MW( ACL_IDX )
         DFH2OR8    = 1.0D6 * MWWAT      ! aerospc_mw(AH2O_IDX) (18.0 != 18.0153)

         ! Mol. wt ratios NH3/NH4, HNO3/NO3, HCL/CL
         NH3RAT  = PRECURSOR_MW( NH3_IDX )  / AEROSPC_MW( ANH4_IDX )
         HNO3RAT = PRECURSOR_MW( HNO3_IDX ) / AEROSPC_MW( ANO3_IDX )
         HCLRAT  = PRECURSOR_MW( HCL_IDX )  / AEROSPC_MW( ACL_IDX )

      END IF

C *** Determine if Hybrid
      TMP = 0.0
      DO I = 1, N_AEROSPC
         IF ( AEROSPC( I )%CHARGE .NE. 0 ) TMP = TMP + AEROSPC_CONC( I,N_MODE )
      END DO
      HYBRID = ( TMP .GE. CUTOFF ) .AND. ( AIRRH .GE. 0.18 )

      HYBRID = .FALSE.    !debug, mdt

      DELT  = DBLE( DT )
      CONVFAC = DELT * H2SO4RATM1

      TEMPI = AIRTEMP             ! assume const within synch timestep
      RHI   = MIN( 0.95, AIRRH )  ! "        "     "      "     "

C *** Calculate molecular velocities (temperature dependent) and
C     H+ concentration

      CBAR_SO4 = COFCBAR_SO4 * SQRT( AIRTEMP )

      HPLUS = 0.0
      DO I = 1, N_MODE
         DO N = 1, N_AEROSPC
            HPLUS( I ) = HPLUS( I )
     &                 + AEROSPC( N )%CHARGE * AEROSPC_CONC( N,I ) / AEROSPC_MW( N )
         END DO
      END DO

C *** Condensational Growth (Size-dependent terms)
C     Calculate intermediate variables needed to determine the 2nd and
C     3rd moment condensational-growth rates.  3rd moment terms are 
C     needed for the calculation of new particle production.  See 
C     Section 3.3 of Jiang & Roth (2003) for a detailed discussion.
C    
C *** Calculate first moments using Equation 4 of Binkowski & Shankar
C     (1995) or Equation 3 of Binkowski and Roselle (2003).
C     N.B: these are for a "wet" size distribution

      DO I = 1, N_MODE
         AM0( I ) = MOMENT0_CONC( I ) 
         AM1( I ) = MOMENT0_CONC( I ) * AEROMODE_DIAM( I )
     &            * EXP(0.5 * AEROMODE_SDEV( I ) * AEROMODE_SDEV( I ) )
         AM2( I ) = MOMENT2_CONC( I )
      END DO
      
C *** Calculate the size-dependent terms in the condensational-
C     growth factor expressions for sulfate using 
C     Equations A13-A14 of Binkowski & Shankar (1995). 

      DO I = 1, N_MODE

         ! Checkpointing
         AM0_CHECK1( I ) = AM0( I )
         AM1_CHECK1( I ) = AM1( I )
         AM2_CHECK1( I ) = AM2( I )

         CALL HCOND3( AM0( I ), AM1( I ), AM2( I ),
     &                DV_SO4, ALPHSULF, CBAR_SO4, FCONC_SO4( I,: ) )

         FCONC_SO4_CHECK1( I, : ) = FCONC_SO4( I, : )
      END DO 

      IF ( .NOT.HYBRID ) THEN
         FCONC_SO4( N_MODE,1 ) = 0.D0
         FCONC_SO4( N_MODE,2 ) = 0.D0
      END IF

      DO I = 1, N_MODE
         GRFAC1( I ) = FCONC_SO4( I,1 )
         GRFAC2( I ) = FCONC_SO4( I,2 )
      END DO

C *** New Particle Production
C     Calculate the new particle production rate due to binary
C     nucleation of H2O and H2SO4.  These calculations are performed 
C     only when the gas-phase production rate of H2SO4 (i.e., SO4RATE) 
C     is non-zero.  The condensation rate of H2SO4 is calculated as the
C     gas-phase production rate minus the new particle production rate.

C *** Initialize Variables
      DMDT_SO4  = 0.0D0
      DNDT      = 0.0D0
      DM2DT     = 0.0D0
      SCONDRATE = 0.0D0

C *** Produce new particles only during time steps when the gas-phase 
C     production rate of H2SO4 is non-zero

      IF ( SO4RATE .NE. 0.0D0 ) THEN

C *** Adjust sulfuric acid vapor concentration to a value in
C     equilibrium with the production of new particles and the
C     condensation of sulfuric acid vapor on existing particles, based 
C     on Equations A21 and A23 of Binkowski & Shankar (1995).
         TMP = 0.0
         DO I = 1, N_MODE
            TMP = TMP + FCONC_SO4( I,2 )
         END DO

         XH2SO4 = SO4RATE / TMP
         XH2SO4 = MAX( XH2SO4, CONMIN )
         PRECURSOR_CONC( SULF_IDX ) = XH2SO4

         DNDT_CHECK = DNDT
         DMDT_SO4_CHECK = DMDT_SO4
         DM2DT_CHECK = DM2DT

C *** Calculate new particle production rate for 0th, 2nd, & 3rd moments

         CALL NEWPART3 ( AIRRH, AIRTEMP, XH2SO4, SO4RATE,
     &                   DNDT, DMDT_SO4, DM2DT )
 
C *** Calculate sulfate condensation rate as the gas-phase production 
C     rate minus the new particle production rate, following Equation
C     3.23 of Jiang & Roth (2003).
         SCONDRATE = MAX( DBLE(SO4RATE) - DMDT_SO4, 0.0D0 )

      END IF   ! SO4RATE .NE. 0

C *** Sulfate Condensation (Size-resolved)
C     Calculate rate at which condensing sulfate should be added to each
C     mode.  The "omega" factors are defined in Equations 7a and 7b of
C     Binkowski & Shankar (1995). The i-mode and j-mode factors are 
C     calculated using Equation A17 of Binkowski & Shankar (1995). The 
C     condensation rate for accumulation mode (fine-equilibrium scheme) or 
C     coarse mode (hybrid and dynamic schemes) is computed by difference, 
C     to avoid mass conservation violations arising from numerical error.

      TMP = 0.0
      DO I = 1, N_MODE
         TMP = TMP + FCONC_SO4( I,2 )
      END DO

      FCONCM1_SO4  = 1.0D0 / TMP
      OMEGA_AT_SO4 = FCONCM1_SO4 * FCONC_SO4( 1,2 )
      OMEGA_AC_SO4 = FCONCM1_SO4 * FCONC_SO4( 2,2 )

C *** Growth values for mode merge condition
      CGR( 1 ) = SO4FAC * SCONDRATE * OMEGA_AT_SO4
      CGR( 2 ) = SO4FAC * SCONDRATE * OMEGA_AC_SO4

      CONDSO4( 1 ) = OMEGA_AT_SO4 * SCONDRATE      
      
      IF ( HYBRID ) THEN 
         CONDSO4( 2 ) = OMEGA_AC_SO4 * SCONDRATE      
         CONDSO4( 3 ) = SCONDRATE - ( CONDSO4( 1 ) + CONDSO4( 2 ) ) 
      ELSE                                  ! fine equilibrium
         CONDSO4( 2 ) = SCONDRATE - CONDSO4( 1 )
         CONDSO4( 3 ) = 0.0D0               ! no coarse mode chemistry
      END IF

C *** For Hybrid approach, calculate dynamic mass trasfer for
C     semi-volatile components of coarse mode (NO3, CL, NH4)  

      IF ( HYBRID ) THEN 

         INT_TIME = 0.0D0
         TSTEP    = 90.0D0
         ISTEP    = 1
         IMODE    = 3
         TrustIso = .TRUE.
 
         DO WHILE ( INT_TIME .LT. DELT ) 

            COUNTER = COUNTER + 1

            IF ( INT_TIME + TSTEP .GT. DELT ) TSTEP = DELT - INT_TIME 
            INT_TIME = INT_TIME + TSTEP 
            ISTEP = ISTEP + 1   

            IF ( ISTEP .GT. 1 ) THEN

C *** Calculate first moments using Equation 4 of Binkowski & Shankar
C     (1995) or Equation 3 of Binkowski and Roselle (2003).
C     N.B: these are for a "wet" size distribution

               AM0( 3 ) = MOMENT0_CONC( 3 )
               AM1( 3 ) = MOMENT0_CONC( 3 ) * AEROMODE_DIAM( 3 )
     &                  * EXP( 0.5 * AEROMODE_SDEV( 3 ) * AEROMODE_SDEV( 3 ) )
               AM2( 3 ) = MOMENT2_CONC( 3 )

C *** Calculate the size-dependent terms in the condensational-
C     growth factor expressions for sulfate using 
C     Equations A13-A14 of Binkowski & Shankar (1995). 
               
               ! Checkpointing
               AM0_CHECK( COUNTER ) = AM0( 3 )
               AM1_CHECK( COUNTER ) = AM1( 3 )
               AM2_CHECK( COUNTER ) = AM2( 3 )

               CALL HCOND3( AM0( 3 ), AM1( 3 ), AM2( 3 ), DV_SO4, ALPHSULF, 
     &                      CBAR_SO4, FCONC_SO4( 3,: ) )  ! adapted from Eq A14

               FCONC_SO4_CHECK( COUNTER, : ) = FCONC_SO4( 3, : )

               GRFAC1( 3 ) = FCONC_SO4( 3,1 ) 
               GRFAC2( 3 ) = FCONC_SO4( 3,2 ) 

            END IF               ! if ISTEP .GT. 1

C *** Set conc array to aerosol concentrations prior to mass transfer

            CINORG( KNH4,IMODE ) = AEROSPC_CONC( ANH4_IDX,IMODE )
            CINORG( KNO3,IMODE ) = AEROSPC_CONC( ANO3_IDX,IMODE )
            CINORG( KCL, IMODE ) = AEROSPC_CONC( ACL_IDX,IMODE )
            CINORG( KSO4,IMODE ) = AEROSPC_CONC( ASO4_IDX,IMODE )
            CINORG( KNA, IMODE ) = AEROSPC_CONC( ANA_IDX,IMODE )
            CINORG( KHP, IMODE ) = HPLUS( IMODE )

            M3OTHR = SOILFAC * AEROSPC_CONC( ASOIL_IDX,IMODE )
     &             + ANTHFAC * AEROSPC_CONC( ACORS_IDX,IMODE )
            M3SOA  = 0.0D0
            WETM3  = MOMENT3_CONC( IMODE )
            WETM2  = MOMENT2_CONC( IMODE )
            DRYM3  = WETM3 - DBLE( H2OFAC * AEROSPC_CONC( AH2O_IDX, IMODE ) )
            DRYM20 = WETM2 * ( DRYM3 / WETM3 ) ** D_TWOTHIRDS

            ! Checkpointing
            WETM3_CHECK3( COUNTER ) = WETM3
            DRYM3_CHECK3( COUNTER ) = DRYM3
            WETM2_CHECK3( COUNTER ) = WETM3

C *** Initial vapor-phase concentrations [ug/m3]
            GNO3R8 = PRECURSOR_CONC( HNO3_IDX )
            GNH3R8 = PRECURSOR_CONC( NH3_IDX )
            GCLR8  = PRECURSOR_CONC( HCL_IDX )

C *** This section of code calculates the distribution of ammonia/
C     ammonium, nitric acid/nitrate, and water between the gas and 
C     aerosol phases as a function of total sulfate, total ammonia, 
C     total nitrate, relative humidity, and temperature.  It is assumed
C     that the aerosol is entirely aqueous (i.e., metastable assumption),
C     irrespective of ambient relative humidity.

C *** Compute sulfate production rate [ug/m3 s] for coarse mode

            RATE = CONDSO4( IMODE )
            SO4  = CINORG( KSO4,IMODE ) + RATE * TSTEP * H2SO4RATM1

            IF ( TrustIso ) THEN

C *** Double Precision vars for ISORROPIA [mole/m3]
               WI( 1 ) = CINORG( KNA, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) )
               WI( 2 ) =                  SO4 * ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) )
               WI( 3 ) = CINORG( KNH4,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) )
               WI( 4 ) = CINORG( KNO3,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) )
               WI( 5 ) = CINORG( KCL, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ACL_IDX ) )

               CNTRL( 1 ) = 1.0D0 ! reverse problem
               CNTRL( 2 ) = 1.0D0 ! aerosol in metastable state

C *** Set flags to account for mass conservation violations in ISRP3F

               TRUSTCL  = .TRUE.
               TRUSTNH4 = .TRUE.
               IF ( WI( 1 ) + WI( 5 ) .LT. 1.0D-20 ) THEN
                  TRUSTCL = .FALSE.
               ELSE
                  IF ( WI( 3 ) .LT. 1.0D-10 ) TRUSTNH4 = .FALSE.
                  IF ( WI( 5 ) .LT. 1.0D-10 ) TRUSTCL  = .FALSE.
               END IF

C     WI( 1 ) = NA    [mol/m3]
C     WI( 2 ) = SO4      "
C     WI( 3 ) = NH4      "
C     WI( 4 ) = NO3      "
C     WI( 5 ) = CL       "
C     GAS(1) = NH3, GAS(2) = HNO3, GAS(3) = HCl

               ! Checkpointing
               WI_CHECK2( COUNTER, : ) = WI
               CNTRL_CHECK1( COUNTER, : ) = CNTRL

               CALL ISOROPIA( WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,  
     &                        AERSLD, SCASI, OTHER, TrustIso )

               WT_CHECK1( COUNTER, : ) = WT
               GAS_CHECK2( COUNTER, : ) = GAS
               AERLIQ_CHECK1( COUNTER, : ) = AERLIQ
               AERSLD_CHECK1( COUNTER, : ) = AERSLD
               OTHER_CHECK1( COUNTER, : ) = OTHER

               IF ( GAS( 1 ) .LT. 0.0D0 .OR. GAS( 2 ) .LT. 0.0D0 .OR.
     &              GAS( 3 ) .LT. 0.0D0 ) THEN
                  IF ( FIRSTWRITE ) THEN
                     FIRSTWRITE = .FALSE.
                  END IF
                  TrustIso = .FALSE.
               END IF

            END IF   ! TrustIso

C *** Change in volatile inorganic PM concentration to achieve
C     equilibrium, calculated as initial-gas-phase concentration minus 
C     equilibrium gas-phase concentration. DVOLINORG is positive for
C     condensation and negative for evaporation.

            DVOLINORG( KNH4 ) = GNH3R8 * ( 1.0D-6 / PRECURSOR_MW( NH3_IDX ) ) - GAS( 1 )  
            DVOLINORG( KNO3 ) = GNO3R8 * ( 1.0D-6 / PRECURSOR_MW( HNO3_IDX ) ) - GAS( 2 ) 
            DVOLINORG( KCL )  = GCLR8  * ( 1.0D-6 / PRECURSOR_MW( HCL_IDX ) ) - GAS( 3 ) 
           
C *** Calculate condensation/evaporation flux for this time step and update 
C     volatile species concentrations.  final aerosol conc set to be  
C     no less than the minimum aerosol conc

            IF ( TrustIso ) THEN
               ! Checkpointing
               RATE_CHECK( COUNTER ) = RATE
               HPLUS_CHECK( COUNTER, : ) = HPLUS
               GRFAC2_CHECK1( COUNTER, : ) = GRFAC2
               GAS_CHECK1( COUNTER, : ) = GAS
               GCLR8_CHECK( COUNTER ) = GCLR8
               GNO3R8_CHECK( COUNTER ) = GNO3R8
               GNH3R8_CHECK( COUNTER ) = GNH3R8

               CALL COMPUTE_FLUX( NVOLINORG, GNH3R8, GNO3R8, GCLR8, KNH4,
     &              KNO3, KCL, GAS( 1:3 ), GRFAC2( IMODE ), AERLIQ( 1 ), RATE, J )

              ! Checkpointing
              J_CHECK1( COUNTER, : ) = J
            ELSE
              J( : ) = 0.0D0
            END IF 

            IF ( J( KNH4 ) * TSTEP * DF( KNH4 ) * NH3RAT .GT. GNH3R8 ) THEN
               J( KNH4 ) = GNH3R8 / ( TSTEP * DF( KNH4 ) * NH3RAT )
            END IF
            IF ( J( KNO3 ) * TSTEP * DF( KNO3 ) * HNO3RAT .GT. GNO3R8 ) THEN
               J( KNO3 ) = GNO3R8 / ( TSTEP * DF( KNO3 ) * HNO3RAT )
            ENDIF
            IF ( J( KCL ) * TSTEP * DF(KCL) * HCLRAT .GT. GCLR8 ) THEN
               J( KCL ) = GCLR8 / ( TSTEP * DF( KCL ) * HCLRAT )
            ENDIF

            ! Checkpointing
            J_CHECK2( COUNTER, : ) = J

C *** Integrate mass transfer equation, convert flux from molar to mass

            DO ISP = 1, NVOLINORG
               CFINAL( ISP,IMODE ) = MAX( 0.0D0,
     &                                    CINORG( ISP,IMODE )
     &                                    + J( ISP ) * TSTEP * DF( ISP ) )
            END DO               

            ! Checkpointing
            CINORG_CHECK( COUNTER, :, : ) = CINORG
            J_CHECK( COUNTER, : ) = J

C *** Calculate updated H+ concentration 

            HPLUS( IMODE ) = 2.0D0 * SO4          / AEROSPC_MW( ASO4_IDX )
     &                     + CFINAL( KNO3,IMODE ) / AEROSPC_MW( ANO3_IDX )
     &                     + CFINAL( KCL, IMODE ) / AEROSPC_MW( ACL_IDX )
     &                     - CFINAL( KNH4,IMODE ) / AEROSPC_MW( ANH4_IDX )
     &                     - CINORG( KNA, IMODE ) / AEROSPC_MW( ANA_IDX )

C *** Equilibrate aerosol LWC with CFINAL by calling CALC_H2O

            WI( 1 ) = CINORG( KNA, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) )
            WI( 2 ) =                  SO4 * ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) )
            WI( 3 ) = CFINAL( KNH4,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) )
            WI( 4 ) = CFINAL( KNO3,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) )
            WI( 5 ) = CFINAL( KCL, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ACL_IDX ) )

            CALC_H2O_STEP = CALC_H2O_STEP + 1

            ! Checkpointing
            WI_CHECK1( CALC_H2O_STEP, : ) = WI
!
            CALL CALC_H2O( WI, RHI, TEMPI, H2O_NEW )

            H2O_NEW_CHECK( CALC_H2O_STEP ) = H2O_NEW

            H2O = H2O_NEW * DFH2OR8 

C *** Compute net change in 3rd moment due to dry inorganic mass transfer

            DDRYM3DT = (
     &                   ( CFINAL( KNH4,IMODE ) - CINORG( KNH4,IMODE ) )
     &                     * ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY)
     &                 + ( CFINAL( KNO3,IMODE ) - CINORG( KNO3,IMODE ) )
     &                     * ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY)
     &                 + ( CFINAL( KCL, IMODE ) - CINORG( KCL,IMODE ) )
     &                     * ( 1.0D-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY)
     &                 + ( SO4                  - CINORG( KSO4,IMODE ) )
     &                     * ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY) ) / TSTEP

C *** Compute net change in 2nd moment due to dry inorganic mass transfer
C     (including nucleation) using equation A7 of Binkowski & Shankar (1995)
            DDRYM2DT = D_TWOTHIRDS * GRFAC1( IMODE ) / GRFAC2( IMODE ) * DDRYM3DT   

            ! Checkpointing
            DDRYM2DT_CHECK( COUNTER ) = DDRYM2DT
            GRFAC1_CHECK( COUNTER, : ) = GRFAC1
            GRFAC2_CHECK( COUNTER, : ) = GRFAC2
            DDRYM3DT_CHECK( COUNTER ) = DDRYM3DT

C *** Update dry 2nd moment for condensation/evaporation based on whether
C     net change in dry 2nd moment is production or loss
            IF ( DDRYM2DT .LT. 0.0D0 ) THEN
               LOSS = DDRYM2DT / DRYM20
               Y = DRYM20 * EXP( LOSS * TSTEP )
               ! Checkpointing
               LOSS_CHECK( COUNTER ) = LOSS
               DRYM20_CHECK( COUNTER ) = DRYM20
            ELSE
               Y = DRYM20 + DDRYM2DT * TSTEP
            ENDIF

C *** Add water and SOA to 2nd moment while preserving standard deviation

            DRYM3 = ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY ) * SO4
     &            + ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY )
     &                               * CFINAL( KNH4,IMODE )
     &            + ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY )
     &                               * CFINAL( KNO3,IMODE )
     &            + ( 1.0D-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY )
     &                               * CFINAL( KCL,IMODE )
     &            + ( 1.0D-9 * F6DPI / AEROSPC( ANA_IDX )%DENSITY )
     &                               * CINORG( KNA,IMODE )
     &            + M3OTHR                   

            WETM3 = DRYM3 + M3SOA + ( 1.0D-9 * F6DPI / AEROSPC( AH2O_IDX )%DENSITY ) * H2O  
            DRYM2 = MAX( DBLE( AEROMODE( IMODE )%MIN_M2 ), Y )
            WETM2 = DRYM2 * ( WETM3 / DRYM3 ) ** D_TWOTHIRDS

            ! Checkpointing
            DRYM2_CHECK2( COUNTER ) = DRYM2
            DRYM3_CHECK2( COUNTER ) = DRYM3
            WETM3_CHECK2( COUNTER ) = WETM3
            Y_CHECK1( COUNTER ) = Y


            PRECURSOR_CONC( NH3_IDX ) = GNH3R8 + ( CINORG( KNH4,IMODE )
     &                                - CFINAL( KNH4,IMODE ) ) * NH3RAT 
 
            PRECURSOR_CONC( HNO3_IDX ) = GNO3R8 + ( CINORG( KNO3,IMODE )
     &                                 - CFINAL( KNO3,IMODE) ) * HNO3RAT  
   
            PRECURSOR_CONC( HCL_IDX ) = GCLR8 + ( CINORG( KCL,IMODE )
     &                                - CFINAL( KCL,IMODE) ) * HCLRAT 

            AEROSPC_CONC( ANH4_IDX,IMODE ) = CFINAL( KNH4,IMODE )
            AEROSPC_CONC( ANO3_IDX,IMODE ) = CFINAL( KNO3,IMODE )
            AEROSPC_CONC( ACL_IDX,IMODE )  = CFINAL( KCL, IMODE )
            AEROSPC_CONC( ASO4_IDX,IMODE ) = SO4
            AEROSPC_CONC( AH2O_IDX,IMODE ) = H2O
            MOMENT2_CONC( IMODE ) = WETM2


C *** Update the third moments, geometric mean diameters, geometric 
C     standard deviations, modal mass totals, and modal particle 
C     densities.
            
            ! Checkpointing
            aeromode_mass_check( COUNTER, : ) = aeromode_mass
            aeromode_sdev_check( COUNTER, : ) = aeromode_sdev
            aeromode_diam_check( COUNTER, : ) = aeromode_diam
            aerospc_conc_check ( COUNTER, :, : ) = aerospc_conc
            moment0_conc_check( COUNTER, : ) = moment0_conc
            moment2_conc_check( COUNTER, : ) = moment2_conc
            moment3_conc_check( COUNTER, : ) = moment3_conc
   
            CALL GETPAR( M3_WET_FLAG, LIMIT_Sg )

         END DO                 ! end mass transfer time step loop
         
      END IF                     ! endif for 'Hybrid' method

C *** End of dynamic mass transfer calculations for coarse mode

         
C *** Call ISORROPIA in forward mode to calculate gas-particle equilibrium in 
C     the fine aerosol modes 

      GNH3R8  = PRECURSOR_CONC( NH3_IDX )
      GNO3R8  = PRECURSOR_CONC( HNO3_IDX )
      GCLR8   = PRECURSOR_CONC( HCL_IDX )

C *** Compute sulfate from total sulfate production rate [ug/m3-s] for fine 
C     modes; add in H2SO4 nucleated in model timestep

      SO4 = AEROSPC_CONC( ASO4_IDX,1 ) + AEROSPC_CONC( ASO4_IDX,2 )
      SO4 = SO4 + ( DMDT_SO4 + CONDSO4( 1 ) + CONDSO4( 2 ) )
     &              * DELT * H2SO4RATM1 

      WI( 1 ) = ( AEROSPC_CONC( ANA_IDX,1 ) + AEROSPC_CONC( ANA_IDX,2 ) )
     &        * ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) )
      WI( 2 ) = SO4 * ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) )
      WI( 3 ) = PRECURSOR_CONC( NH3_IDX ) * ( 1.0D-6 / PRECURSOR_MW( NH3_IDX ) )
     &        + ( AEROSPC_CONC( ANH4_IDX,1 ) + AEROSPC_CONC( ANH4_IDX,2 ) )
     &        * ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) )
      WI( 4 ) = PRECURSOR_CONC( HNO3_IDX ) * ( 1.0D-6 / PRECURSOR_MW( HNO3_IDX ) )
     &        + ( AEROSPC_CONC( ANO3_IDX,1 ) + AEROSPC_CONC( ANO3_IDX,2 ) )
     &        * ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) )

      WI( 5 ) = PRECURSOR_CONC(HCL_IDX) * ( 1.0D-6 / PRECURSOR_MW(HCL_IDX) ) 
     &        + ( AEROSPC_CONC(ACL_IDX,1) + AEROSPC_CONC(ACL_IDX,2) )
     &        * ( 1.0D-6 / AEROSPC_MW(ACL_IDX) )


      CNTRL( 1 ) = 0.0D0   ! Forward Problem
      CNTRL( 2 ) = 1.0D0   ! Aerosol in Metastable State

C *** Set flags to account for mass conservation violations in ISRP3F

      TRUSTCL  = .TRUE.
      TRUSTNH4 = .TRUE.
      IF ( WI( 1 ) + WI( 5 ) .LT. 1.0D-20 ) THEN
         TRUSTCL = .FALSE.
      ELSE
         IF ( WI( 3 ) .LT. 1.0D-10 ) TRUSTNH4 = .FALSE.
         IF ( WI( 5 ) .LT. 1.0D-10 ) TRUSTCL  = .FALSE.
      END IF
     
      ! Checkpointing
      WI_CHECK( : ) = WI

      CALL ISOROPIA( WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,
     &               AERSLD, SCASI, OTHER, TrustIso )

      GAS_CHECK( : ) = GAS
      WT_CHECK( : ) = WT
      AERLIQ_CHECK( : ) = AERLIQ
      AERSLD_CHECK( : ) = AERSLD
      OTHER_CHECK( : ) = OTHER

C *** Update gas-phase concentrations

      PRECURSOR_CONC( NH3_IDX )  = GAS( 1 ) * PRECURSOR_MW( NH3_IDX ) * 1.0E6
      PRECURSOR_CONC( HNO3_IDX ) = GAS( 2 ) * PRECURSOR_MW( HNO3_IDX ) * 1.0E6
      PRECURSOR_CONC( HCL_IDX )  = GAS( 3 ) * PRECURSOR_MW( HCL_IDX ) * 1.0E6

C *** Change in volatile inorganic PM concentration to achieve
C     equilibrium, calculated as initial-gas-phase concentration minus 
C     equilibrium gas-phase concentration. DVOLINORG is positive for
C     condensation and negative for evaporation.

      DVOLINORG( KNH4 ) = GNH3R8 * 1.0D-6 / DBLE( PRECURSOR_MW( NH3_IDX ) ) - GAS( 1 )  
      DVOLINORG( KNO3 ) = GNO3R8 * 1.0D-6 / DBLE( PRECURSOR_MW( HNO3_IDX ) ) - GAS( 2 ) 

      ! Checkpointing
      TRUSTCL_CHECK = TRUSTCL

      IF ( TRUSTCL ) THEN  
         DVOLINORG( KCL )  = GCLR8  * 1.0D-6 / DBLE( PRECURSOR_MW( HCL_IDX ) ) - GAS( 3 ) 
      ELSE
         DVOLINORG( KCL ) = 0.0D0
      END IF

      ! Checkpointing
      DVOLINORG_CHECK = DVOLINORG

      IF ( DVOLINORG( KNH4 ) .LT. 0.0D0 ) THEN
         DVOLMAX = -( DBLE( AEROSPC_CONC( ANH4_IDX,1 ) )
     &              + DBLE( AEROSPC_CONC( ANH4_IDX,2 ) ) ) / DF(KNH4)
         DVOLINORG( KNH4 ) = MAX( DVOLINORG( KNH4 ), DVOLMAX )
      END IF

      IF ( DVOLINORG( KNO3 ) .LT. 0.0D0 ) THEN
         DVOLMAX = -( DBLE( AEROSPC_CONC( ANO3_IDX,1 ) )
     &              + DBLE( AEROSPC_CONC( ANO3_IDX,2 ) ) ) / DF( KNO3 )
         DVOLINORG( KNO3 ) = MAX( DVOLINORG( KNO3 ), DVOLMAX)
      END IF

      IF ( DVOLINORG( KCL ) .LT. 0.0D0 ) THEN
         DVOLMAX = -( DBLE( AEROSPC_CONC( ACL_IDX,1 ) )
     &              + DBLE( AEROSPC_CONC( ACL_IDX,2 ) ) ) / DF( KCL )
         DVOLINORG( KCL ) = MAX( DVOLINORG( KCL ), DVOLMAX )
      END IF

      OMEGA( 1 ) = GRFAC2( 1 ) / ( GRFAC2( 1 ) + GRFAC2( 2 ) ) ! partitioning cof
      OMEGA( 2 ) = 1.0D0 - OMEGA( 1 )

C *** Initialize excess evaporated mass array

      DMASS = 0.0D0

      DO IMODE = 1, 2  ! modal partitioning of equilibrium aerosol mass

         CINORG( KSO4,IMODE ) = AEROSPC_CONC( ASO4_IDX, IMODE )
         CINORG( KNH4,IMODE ) = AEROSPC_CONC( ANH4_IDX, IMODE )
         CINORG( KNO3,IMODE ) = AEROSPC_CONC( ANO3_IDX, IMODE )
         CINORG( KNA, IMODE ) = AEROSPC_CONC( ANA_IDX, IMODE )
         CINORG( KCL, IMODE ) = AEROSPC_CONC( ACL_IDX, IMODE )

         WETM3 = MOMENT3_CONC( IMODE )
         WETM2 = MOMENT2_CONC( IMODE )

         IF ( IMODE .EQ. 1 ) THEN
           RATE = DMDT_SO4 + CONDSO4( IMODE )
         ELSE
           RATE = CONDSO4( IMODE )
         ENDIF

         M3OTHR = 0.0
         M3SOA = 0.0
         DO I = 1, N_AEROSPC

            ! Skip organic species
            IF ( I .EQ. ASO4_IDX .OR. I .EQ. ANH4_IDX .OR. I .EQ. ANO3_IDX .OR.
     &           I .EQ. ANA_IDX  .OR. I .EQ. ACL_IDX  .OR. I .EQ. AH2O_IDX ) CYCLE

            IF ( AEROSPC(I)%isWet ) THEN

               ! Checkpointing
               WET_CHECK( IMODE, I ) = .TRUE.

               M3SOA = M3SOA 
     &               + ( 1.0D-9 * F6DPI / AEROSPC( I )%DENSITY )
     &               * AEROSPC_CONC( I,IMODE )
            ELSE

               ! Checkpointing
               WET_CHECK( IMODE, I ) = .FALSE.

               M3OTHR = M3OTHR
     &                + ( 1.0D-9 * F6DPI / AEROSPC( I )%DENSITY )
     &                * AEROSPC_CONC( I,IMODE )
            ENDIF   


         END DO    ! N_AEROSPC loop
   
         DRYM3 = WETM3 - DBLE( H2OFAC * AEROSPC_CONC( AH2O_IDX,IMODE ) ) - M3SOA
         DRYM20 = WETM2 * ( DRYM3 / WETM3 ) ** D_TWOTHIRDS

         ! Checkpointing
         WETM3_CHECK1( IMODE ) = WETM3
         DRYM3_CHECK1( IMODE ) = DRYM3
         WETM2_CHECK1( IMODE ) = WETM2

         DO ISP = 1, NVOLINORG

            CFINAL( ISP,IMODE ) = CINORG( ISP,IMODE )
     &                          + OMEGA( IMODE ) * DVOLINORG( ISP )
     &                          * DF( ISP )

            ! Checkpointing
            CFINAL_CHECK( ISP, IMODE ) = CFINAL( ISP, IMODE )

            IF ( IMODE .EQ. 1 ) THEN

               IF ( CFINAL( ISP,IMODE ) .LT. 0.0D0 ) THEN
                  DMASS( ISP) = CFINAL( ISP,IMODE )
                  CFINAL( ISP,IMODE ) = 0.0D0
               END IF

            ELSE

               CFINAL( ISP,IMODE ) = CFINAL( ISP,IMODE ) + DMASS( ISP )

               ! Checkpointing
               CFINAL_CHECK1( ISP, IMODE ) = CFINAL( ISP, IMODE )

               IF ( CFINAL( ISP,IMODE ) .LT. 0.0D0 ) THEN
                  CFINAL( ISP,1 ) = CFINAL( ISP,1 ) + CFINAL( ISP,IMODE ) 
                  CFINAL( ISP,IMODE ) = 0.0D0

                  IF ( CFINAL( ISP, 1 ) .LT. 0.0D0 ) THEN
                     IF ( ABS( CFINAL( ISP,1 ) ) .LT. 1.0D-15 ) THEN
                        CFINAL( ISP,1 ) = 0.0D0
                     ELSE
                        PRINT *, 'Too much evaporation: aero_subs.f'
                        CFINAL( ISP,1 ) = 0.0D0
C                       STOP
                     END IF
                  END IF 
               END IF

            END IF   ! IMODE = 1
            
         END DO 
                        
         SO4 = CINORG( KSO4,IMODE ) + RATE * CONVFAC

C *** Double precision vars for ISORROPIA

         WI( 1 ) = CINORG( KNA, IMODE ) * 1.0D-6 / AEROSPC_MW( ANA_IDX )
         WI( 2 ) =                  SO4 * 1.0D-6 / AEROSPC_MW( ASO4_IDX )
         WI( 3 ) = CFINAL( KNH4,IMODE ) * 1.0D-6 / AEROSPC_MW( ANH4_IDX )
         WI( 4 ) = CFINAL( KNO3,IMODE ) * 1.0D-6 / AEROSPC_MW( ANO3_IDX )
         WI( 5 ) = CFINAL( KCL, IMODE ) * 1.0D-6 / AEROSPC_MW( ACL_IDX )

         CALC_H2O_STEP = CALC_H2O_STEP + 1

         ! Checkpointing
         WI_CHECK4( IMODE, : ) = WI
  
         CALL CALC_H2O( WI, RHI, TEMPI, H2O_NEW )
 
         H2O = H2O_NEW * DFH2OR8 


C *** Compute net change in 3rd moment due to dry inorganic mass transfer
C     (includes nucleated sulfate mass)

         DDRYM3DT = (
     &                ( CFINAL( KNH4,IMODE ) - CINORG( KNH4,IMODE ) )
     &                  * ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY )
     &              + ( CFINAL( KNO3,IMODE ) - CINORG( KNO3,IMODE ) )
     &                  * ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY )
     &              + ( CFINAL( KCL, IMODE ) - CINORG( KCL,IMODE ) )
     &                  * ( 1.0D-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY )
     &              + ( SO4                  - CINORG( KSO4,IMODE ) )
     &                  * ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY ) ) / DELT 

C *** Compute net change in 2nd moment due to dry inorganic mass transfer
C     (including nucleation) using equation A7 of Binkowski & Shankar (1995)
         DDRYM2DT = D_TWOTHIRDS * GRFAC1( IMODE ) / GRFAC2( IMODE ) * DDRYM3DT

      ! Checkpointing
         DDRYM3DT_CHECK1( IMODE ) = DDRYM3DT

C *** Update dry 2nd moment for condensation/evaporation based on whether
C     net change in dry 2nd moment is production or loss

         ! Checkpointing
         DDRYM2DT_CHECK1( IMODE ) = DDRYM2DT
         DRYM20_CHECK1( IMODE ) = DRYM20

         IF ( DDRYM2DT .LT. 0.0D0 ) THEN
            LOSS = DDRYM2DT / DRYM20
            Y = DRYM20 * EXP( LOSS * DELT )

            ! Checkpointing
            LOSS_CHECK1( IMODE ) = LOSS
 
         ELSE
            Y = DRYM20 + DDRYM2DT * DELT
         END IF

C *** Add water and SOA to 2nd moment while preserving standard deviation

         DRYM3 = ( 1.0E-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY ) * SO4
     &         + ( 1.0E-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY ) * CFINAL( KNH4,IMODE )
     &         + ( 1.0E-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY ) * CFINAL( KNO3,IMODE )
     &         + ( 1.0E-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY )  * CFINAL( KCL, IMODE )
     &         + ( 1.0D-9 * F6DPI / AEROSPC( ANA_IDX )%DENSITY )  * CINORG( KNA, IMODE )
     &         + M3OTHR

         WETM3 = DRYM3 + M3SOA + DBLE( H2OFAC ) * H2O 
         DRYM2 = MAX( DBLE( AEROMODE(IMODE)%MIN_M2 ), Y )
         WETM2 = DRYM2 * ( WETM3 / DRYM3 ) ** D_TWOTHIRDS

         ! Checkpointing
         DRYM2_CHECK( IMODE ) = DRYM2
         WETM3_CHECK( IMODE ) = WETM3
         DRYM3_CHECK( IMODE ) = DRYM3
         Y_CHECK( IMODE ) = Y

         AEROSPC_CONC( ANH4_IDX, IMODE ) = CFINAL( KNH4,IMODE )
         AEROSPC_CONC( ANO3_IDX, IMODE ) = CFINAL( KNO3,IMODE )
         AEROSPC_CONC( ACL_IDX, IMODE )  = CFINAL( KCL ,IMODE )
         AEROSPC_CONC( ASO4_IDX, IMODE ) = SO4
         AEROSPC_CONC( AH2O_IDX, IMODE ) = H2O
         if(IMODE.eq.1) MOMENT0_CONC( IMODE ) = AM0(IMODE) + DNDT * DELT
         MOMENT2_CONC( IMODE ) = WETM2


         HPLUS( IMODE ) = 0.0
         DO I = 1, N_AEROSPC
           HPLUS( IMODE ) = HPLUS( IMODE )
     &                    + AEROSPC( I )%CHARGE * AEROSPC_CONC( I,IMODE ) / AEROSPC_MW( I )
         END DO

      END DO    ! end fine mode mass transfer calculations

C *** Update the third moments, geometric mean diameters, geometric 
C     standard deviations, modal mass totals, and modal particle 
C     densities.

      AEROMODE_MASS_CHK = AEROMODE_MASS
      AEROMODE_SDEV_CHK = AEROMODE_SDEV
      AEROMODE_DIAM_CHK = AEROMODE_DIAM
      AEROSPC_CONC_CHK  = AEROSPC_CONC
      MOMENT3_CONC_CHK = MOMENT3_CONC
      MOMENT2_CONC_CHK = MOMENT2_CONC
      MOMENT0_CONC_CHK = MOMENT0_CONC

      CALL GETPAR( M3_WET_FLAG, LIMIT_Sg )

! Adjoint Code:

      AEROMODE_MASS = AEROMODE_MASS_CHK
      AEROMODE_SDEV = AEROMODE_SDEV_CHK
      AEROMODE_DIAM = AEROMODE_DIAM_CHK
      MOMENT3_CONC = MOMENT3_CONC_CHK
      MOMENT2_CONC = MOMENT2_CONC_CHK
      MOMENT0_CONC = MOMENT0_CONC_CHK
      AEROSPC_CONC = AEROSPC_CONC_CHK

      !------
      ! fwd code:
      ! CALL GETPAR( M3_WET_FLAG, LIMIT_Sg )
      ! adj code:
      CALL GETPAR_B( M3_WET_FLAG, LIMIT_Sg )

      DO IMODE = 2, 1, -1

         DO I = N_AEROSPC, 1, -1
            !-----
            ! fwd code:
            ! HPLUS( IMODE ) = HPLUS( IMODE ) + AEROSPC( I )%CHARGE
            !                * AEROSPC_CONC( I, IMODE ) / AEROSPC_MW( I )  
            ! adj code:
            aerospc_concb( I, IMODE ) = aerospc_concb( I, IMODE ) +
     &             ( AEROSPC( I )%CHARGE / AEROSPC_MW( I ) ) * 
     &             HPLUS_ADJ ( IMODE )
         END DO

         !------
         ! fwd code:
         ! HPLUS( IMODE ) = 0.0
         ! adj code:
         HPLUS_ADJ ( IMODE ) = 0.0

         !------
         ! fwd code:
         ! MOMENT2_CONC( IMODE ) = WETM2
         ! adj code:
         WETM2_ADJ = moment2_concb( IMODE )
         moment2_concb( IMODE ) = 0.0
 
         IF ( IMODE == 1 ) THEN
            !------
            ! fwd code:
            ! MOMENT0_CONC( IMODE ) = AM0(IMODE) + DNDT * DELT
            ! adj code:
            AM0_ADJ( IMODE ) = AM0_ADJ( IMODE ) + moment0_concb( IMODE )
            DNDT_ADJ = DNDT_ADJ + DELT * moment0_concb( IMODE )
            moment0_concb( IMODE ) = 0.0
         END IF

         !------
         ! fwd code:
         ! AEROSPC_CONC( AH2O_IDX, IMODE ) = H2O
         ! adj code:
         H2O_ADJ = aerospc_concb( AH2O_IDX, IMODE )
         aerospc_concb( AH2O_IDX, IMODE ) = 0.0

         !------ 
         ! fwd code:
         ! AEROSPC_CONC( ASO4_IDX, IMODE ) = SO4
         ! adj code:
         SO4_ADJ = aerospc_concb( ASO4_IDX, IMODE ) 
         aerospc_concb( ASO4_IDX, IMODE ) = 0.0

         !------
         ! fwd code:
         ! AEROSPC_CONC( ACL_IDX, IMODE )  = CFINAL( KCL ,IMODE )
         ! adj code:
         CFINAL_ADJ ( KCL, IMODE ) = CFINAL_ADJ ( KCL, IMODE ) + 
     &                               aerospc_concb( ACL_IDX, IMODE )
         aerospc_concb( ACL_IDX, IMODE ) = 0.0

         !------
         ! fwd code:
         ! AEROSPC_CONC( ANO3_IDX, IMODE ) = CFINAL( KNO3,IMODE )
         ! adj code:
         CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) 
     &                              + aerospc_concb( ANO3_IDX, IMODE )
         aerospc_concb( ANO3_IDX, IMODE ) = 0.0

         !------
         ! fwd code:
         ! AEROSPC_CONC( ANH4_IDX, IMODE ) = CFINAL( KNH4,IMODE )
         ! adj code:
         CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE )
     &                              + aerospc_concb( ANH4_IDX, IMODE )
         aerospc_concb( ANH4_IDX, IMODE ) = 0.0

         ! Checkpointing
         WETM3 = WETM3_CHECK( IMODE )
         DRYM3 = DRYM3_CHECK( IMODE )
         DRYM2 = DRYM2_CHECK( IMODE )

         !------
         ! fwd code:
         !  WETM2 = DRYM2 * ( WETM3 / DRYM3 ) ** D_TWOTHIRDS
         ! adj code:
         TEMP4 = WETM3 / DRYM3
!         IF ( TEMP4 .LE. 0.0 ) THEN
!            TEMP4B5 = 0.0
!         ELSE
            TEMP4B5 = DRYM2 * D_TWOTHIRDS * TEMP4 ** ( D_TWOTHIRDS - 1 )
     &              * WETM2_ADJ / DRYM3
!         END IF
         DRYM2_ADJ = TEMP4 ** D_TWOTHIRDS * WETM2_ADJ 
         WETM3_ADJ = TEMP4B5 
         DRYM3_ADJ = -( TEMP4 * TEMP4B5 )

         ! Checkpointing
         Y = Y_CHECK( IMODE )

         !------
         ! fwd code:
         ! DRYM2 = MAX( DBLE( AEROMODE(IMODE)%MIN_M2 ), Y )
         ! adj code:
         IF ( DBLE( AEROMODE(IMODE)%MIN_M2 ) < Y ) THEN
            Y_ADJ = DRYM2_ADJ
            DRYM2_ADJ = 0.0
         ELSE
            DRYM2_ADJ = 0.0
         END IF

         !------
         ! fwd code:
         ! WETM3 = DRYM3 + M3SOA + DBLE( H2OFAC ) * H2O
         ! adj code:
         DRYM3_ADJ = DRYM3_ADJ + WETM3_ADJ 
         M3SOA_ADJ = WETM3_ADJ
         H2O_ADJ = H2O_ADJ + DBLE( H2OFAC ) * WETM3_ADJ
         WETM3_ADJ = 0.0

         !------
         ! fwd code: 
         ! DRYM3 = ( 1.0E-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY ) * SO4
         ! &     + ( 1.0E-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY ) * CFINAL( KNH4,IMODE )
         ! &     + ( 1.0E-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY ) * CFINAL( KNO3,IMODE )
         ! &     + ( 1.0E-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY )  * CFINAL( KCL, IMODE )
         ! &     + ( 1.0D-9 * F6DPI / AEROSPC( ANA_IDX )%DENSITY )  * CINORG( KNA, IMODE )
         ! &     + M3OTHR
         ! adj code:
         SO4_ADJ = SO4_ADJ + ( 1.0E-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY ) 
     &           * DRYM3_ADJ
         CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ( KNH4, IMODE ) +
     &           ( 1.0E-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY ) * DRYM3_ADJ
         CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) +
     &           ( 1.0E-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY ) * DRYM3_ADJ
         CFINAL_ADJ( KCL, IMODE ) = CFINAL_ADJ( KCL, IMODE ) +
     &           ( 1.0E-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY ) * DRYM3_ADJ
         CINORG_ADJ ( KNA, IMODE ) = CINORG_ADJ( KNA, IMODE ) +
     &           ( 1.0D-9 * F6DPI / AEROSPC( ANA_IDX )%DENSITY ) * DRYM3_ADJ
         M3OTHR_ADJ = DRYM3_ADJ

         ! Checkpointing
         DRYM20 = DRYM20_CHECK1( IMODE )
         DDRYM2DT = DDRYM2DT_CHECK1( IMODE )

         IF ( DDRYM2DT .LT. 0.0D0 ) THEN

            ! Checkpointing
            LOSS = LOSS_CHECK1( IMODE )

            !------
            ! fwd code:
            ! LOSS = DDRYM2DT / DRYM20
            ! Y = DRYM20 * EXP( LOSS * DELT )
            ! adj code:
            DRYM20_ADJ = EXP ( LOSS * DELT ) * Y_ADJ
            LOSS_ADJ = DELT * DRYM20 * EXP ( LOSS * DELT )
     &               * Y_ADJ
            Y_ADJ = 0.0
            DDRYM2DT_ADJ = ( 1.0 / DRYM20 ) * LOSS_ADJ
            DRYM20_ADJ = DRYM20_ADJ - ( DDRYM2DT / ( DRYM20 ** 2.0 ) )
     &                 * LOSS_ADJ
            LOSS_ADJ = 0.0
         ELSE
            !------
            ! fwd code:
            ! Y = DRYM20 + DDRYM2DT * DELT
            ! adj code:
            DRYM20_ADJ = Y_ADJ
            DDRYM2DT_ADJ = DELT * Y_ADJ
            Y_ADJ = 0.0
         END IF

         ! Checkpointing
         DDRYM3DT = DDRYM3DT_CHECK1( IMODE )

         !------
         ! fwd code:
         ! DDRYM2DT = D_TWOTHIRDS * GRFAC1( IMODE ) / GRFAC2( IMODE ) * DDRYM3DT
         ! adj code:
         GRFAC1_ADJ ( IMODE ) = GRFAC1_ADJ ( IMODE ) + D_TWOTHIRDS
     &           * ( 1.0 / GRFAC2( IMODE ) ) * DDRYM3DT * DDRYM2DT_ADJ
         DDRYM3DT_ADJ = D_TWOTHIRDS * ( GRFAC1 ( IMODE ) 
     &           / GRFAC2( IMODE ) ) * DDRYM2DT_ADJ
         GRFAC2_ADJ ( IMODE ) = GRFAC2_ADJ ( IMODE ) - D_TWOTHIRDS *
     &           DDRYM3DT * GRFAC1( IMODE ) / ( GRFAC2 ( IMODE ) ** 2.0 )
     &           * DDRYM2DT_ADJ
         DDRYM2DT_ADJ = 0.0

         !------
         ! fwd code:
         ! DDRYM3DT = (
         ! &     ( CFINAL( KNH4,IMODE ) - CINORG( KNH4,IMODE ) )
         ! &       * ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY )
         ! &   + ( CFINAL( KNO3,IMODE ) - CINORG( KNO3,IMODE ) )
         ! &       * ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY )
         ! &   + ( CFINAL( KCL, IMODE ) - CINORG( KCL,IMODE ) )
         ! &       * ( 1.0D-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY )
         ! &   + ( SO4                  - CINORG( KSO4,IMODE ) )
         ! &       * ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY ) ) / DELT
         ! adj code:
         CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) + 
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ
         CINORG_ADJ ( KNH4, IMODE ) = CINORG_ADJ ( KNH4, IMODE ) -
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ
         CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) +
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ
         CINORG_ADJ ( KNO3, IMODE ) = CINORG_ADJ ( KNO3, IMODE ) -
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ
         CFINAL_ADJ ( KCL , IMODE ) = CFINAL_ADJ ( KCL , IMODE ) +
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC(  ACL_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ
         CINORG_ADJ ( KCL , IMODE ) = CINORG_ADJ ( KCL , IMODE ) -
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC(  ACL_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ
         SO4_ADJ = SO4_ADJ + 
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ
         CINORG_ADJ ( KSO4, IMODE ) = CINORG_ADJ ( KSO4, IMODE ) -
     &             ( 1.0 / DELT ) * ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX
     &             )%DENSITY ) * DDRYM3DT_ADJ         
         DDRYM3DT_ADJ = 0.0

         !------
         ! fwd code:
         ! H2O = H2O_NEW * DFH2OR8
         ! adj code:
         H2O_NEW_ADJ = H2O_ADJ * DFH2OR8
         H2O_ADJ = 0.0

         ! Checkpointing
         WI = WI_CHECK4( IMODE, : )

         WI_0 = WI

         WI_ADJ_SAVE = WI_ADJ

         CALL CALC_H2O_ADJ( WI, WI_ADJ, RHI, TEMPI, H2O_NEW, H2O_NEW_ADJ )

         CALC_H2O_STEP = CALC_H2O_STEP - 1

         !------
         ! fwd code:
         ! WI( 1 ) = CINORG( KNA, IMODE ) * 1.0D-6 / AEROSPC_MW( ANA_IDX )
         ! WI( 2 ) =                  SO4 * 1.0D-6 / AEROSPC_MW( ASO4_IDX )
         ! WI( 3 ) = CFINAL( KNH4,IMODE ) * 1.0D-6 / AEROSPC_MW( ANH4_IDX )
         ! WI( 4 ) = CFINAL( KNO3,IMODE ) * 1.0D-6 / AEROSPC_MW( ANO3_IDX )
         ! WI( 5 ) = CFINAL( KCL, IMODE ) * 1.0D-6 / AEROSPC_MW( ACL_IDX )
         ! adj code:
         CFINAL_ADJ ( KCL, IMODE ) = CFINAL_ADJ ( KCL, IMODE ) +
     &         ( 1.0D-6 / AEROSPC_MW( ACL_IDX ) ) * WI_ADJ ( 5 )
         WI_ADJ ( 5 ) = 0.0    
         CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) +
     &         ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) ) * WI_ADJ ( 4 )
         WI_ADJ ( 4 ) = 0.0
         CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) +
     &         ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) ) * WI_ADJ ( 3 )
         WI_ADJ ( 3 ) = 0.0
         SO4_ADJ = SO4_ADJ +
     &         ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) ) * WI_ADJ( 2 )
         WI_ADJ ( 2 ) = 0.0
         CINORG_ADJ( KNA, IMODE ) = CINORG_ADJ( KNA, IMODE ) +
     &         ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) ) * WI_ADJ ( 1 )
         WI_ADJ ( 1 ) = 0.0

         !------
         ! fwd code:
         ! SO4 = CINORG( KSO4,IMODE ) + RATE * CONVFAC
         ! adj code:
         CINORG_ADJ ( KSO4, IMODE ) = CINORG_ADJ( KSO4, IMODE ) +
     &         SO4_ADJ 
         RATE_ADJ = RATE_ADJ + CONVFAC * SO4_ADJ
         SO4_ADJ = 0.0

         DO ISP = NVOLINORG, 1, -1
           
            ! Checkpointing 
            CFINAL( ISP, IMODE ) = CFINAL_CHECK( ISP, IMODE )

            IF ( IMODE .EQ. 1 ) THEN


               IF ( CFINAL( ISP,IMODE ) .LT. 0.0D0 ) THEN
                  !------
                  ! fwd code:
                  ! DMASS( ISP) = CFINAL( ISP,IMODE )
                  ! CFINAL( ISP,IMODE ) = 0.0D0
                  ! adj code:
                  CFINAL_ADJ ( ISP, IMODE ) =  DMASS_ADJ ( ISP )
                  DMASS_ADJ ( ISP ) = 0.0
               END IF
            ELSE

               ! Checkpointing
               CFINAL( ISP, IMODE ) = CFINAL_CHECK1( ISP, IMODE ) 

               IF ( CFINAL( ISP, IMODE ) .LT. 0.0D0 ) THEN
                  IF ( CFINAL( ISP,1 ) .LT. 0.0D0 ) THEN
                     IF ( ABS( CFINAL( ISP,1 ) ) .LT. 1.0D-15 ) THEN
                        !------
                        ! fwd code:
                        ! CFINAL( ISP,1 ) = 0.0D0
                        ! adj code:
                        CFINAL_ADJ ( ISP, 1 ) = 0.0
                     ELSE
                        !------
                        ! fwd code:
                        ! CFINAL( ISP,1 ) = 0.0D0
                        ! adj code:
                        PRINT *, 'Too much evaporation: aero_subs.f'
                        CFINAL_ADJ ( ISP, 1 ) = 0.0
                        STOP
                     END IF
                  END IF

                  !------
                  ! fwd code:
                  ! CFINAL( ISP,1 ) = CFINAL( ISP,1 ) + CFINAL( ISP,IMODE ) 
                  ! CFINAL( ISP, IMODE ) = 0.0
                  ! adj code:
                  CFINAL_ADJ( ISP, IMODE ) = CFINAL_ADJ ( ISP, 1 ) 

               END IF

               !------
               ! fwd code:
               ! CFINAL( ISP,IMODE ) = CFINAL( ISP,IMODE ) + DMASS( ISP )
               ! adj code:
               DMASS_ADJ ( ISP ) = DMASS_ADJ ( ISP ) + CFINAL_ADJ 
     &              ( ISP, IMODE )

            END IF   ! IMODE = 1

            !------
            ! fwd code:
            ! CFINAL( ISP,IMODE ) = CINORG( ISP,IMODE )
            ! &             + OMEGA( IMODE ) * DVOLINORG( ISP )
            ! &             * DF( ISP )
            ! adj code:

            CINORG_ADJ ( ISP, IMODE ) = CINORG_ADJ ( ISP, IMODE ) + 
     &            CFINAL_ADJ ( ISP, IMODE )
            OMEGA_ADJ ( IMODE ) = OMEGA_ADJ( IMODE ) + DVOLINORG( ISP )
     &           * DF( ISP ) * CFINAL_ADJ ( ISP, IMODE )
            DVOLINORG_ADJ ( ISP ) = DVOLINORG_ADJ ( ISP ) + OMEGA( IMODE )
     &           * DF( ISP ) * CFINAL_ADJ ( ISP, IMODE )
            CFINAL_ADJ ( ISP, IMODE ) = 0.0

         END DO

         ! Checkpointing
         WETM3 = WETM3_CHECK1( IMODE )
         DRYM3 = DRYM3_CHECK1( IMODE )
         WETM2 = WETM2_CHECK1( IMODE )

         !------
         ! fwd code:
         ! DRYM3 = WETM3 - DBLE( H2OFAC * AEROSPC_CONC( AH2O_IDX,IMODE ) ) - M3SOA
         ! DRYM20 = WETM2 * ( DRYM3 / WETM3 ) ** D_TWOTHIRDS
         ! adj code: 
         TEMP3 = DRYM3 / WETM3
!         IF ( TEMP3 .LE. 0.0 ) THEN
!            TEMP3B9 = 0.0
!         ELSE
            TEMP3B9 = WETM2 * D_TWOTHIRDS * TEMP3 ** ( D_TWOTHIRDS - 1 )
     &              * DRYM20_ADJ / WETM3 
!         END IF
         WETM2_ADJ = WETM2_ADJ + TEMP3 ** D_TWOTHIRDS * DRYM20_ADJ
         DRYM3_ADJ = DRYM3_ADJ + TEMP3B9
         WETM3_ADJ = WETM3_ADJ - TEMP3 * TEMP3B9
         WETM3_ADJ = WETM3_ADJ + DRYM3_ADJ
         M3SOA_ADJ = M3SOA_ADJ - DRYM3_ADJ
         aerospc_concb( AH2O_IDX, IMODE ) = aerospc_concb( AH2O_IDX, IMODE )
     &          - DBLE( H2OFAC ) * DRYM3_ADJ  

         DO I = N_AEROSPC, 1, -1

            ! Skip organic species
            IF ( I .EQ. ASO4_IDX .OR. I .EQ. ANH4_IDX .OR. I .EQ. ANO3_IDX .OR.
     &           I .EQ. ANA_IDX  .OR. I .EQ. ACL_IDX  .OR. I .EQ. AH2O_IDX ) CYCLE
     
            IF ( AEROSPC(I)%isWet ) THEN
               
               !------
               ! fwd code:
               ! M3SOA = M3SOA
               !     + ( 1.0D-9 * F6DPI / AEROSPC( I )%DENSITY )
               !     * AEROSPC_CONC( I,IMODE )
               ! adj code:
               aerospc_concb( I, IMODE ) = aerospc_concb( I, IMODE ) + 
     &                ( 1.0D-9 * F6DPI / AEROSPC( I )%DENSITY ) *
     &                M3SOA_ADJ
            ELSE
          
               !------
               ! fwd code:
               ! M3OTHR = M3OTHR
               !      + ( 1.0D-9 * F6DPI / AEROSPC( I )%DENSITY )
               !      * AEROSPC_CONC( I,IMODE )
               ! adj code:
               aerospc_concb( I, IMODE ) = aerospc_concb( I, IMODE ) +
     &                ( 1.0D-9 * F6DPI / AEROSPC( I )%DENSITY ) *
     &                M3OTHR_ADJ
            END IF

         END DO

         !------
         ! fwd code:
         ! M3OTHR = 0.0
         ! M3SOA = 0.0
         ! adj code:
         M3SOA_ADJ = 0.0
         M3OTHR_ADJ = 0.0

         IF ( IMODE .EQ. 1 ) THEN
            !------
            ! fwd code:
            ! RATE = DMDT_SO4 + CONDSO4( 1 )
            ! adj code:
            DMDT_SO4_ADJ = DMDT_SO4_ADJ + RATE_ADJ
            CONDSO4_ADJ( IMODE ) = CONDSO4_ADJ( IMODE ) + RATE_ADJ
         ELSE
            !------
            ! fwd code:
            ! RATE = CONDSO4( IMODE )
            ! adj code:
            CONDSO4_ADJ( IMODE ) = CONDSO4_ADJ ( IMODE ) + RATE_ADJ
         END IF

         !------
         ! fwd code:
         ! WETM3 = MOMENT3_CONC( IMODE )
         ! WETM2 = MOMENT2_CONC( IMODE )
         ! adj code:
         moment2_concb( IMODE ) = moment2_concb( IMODE ) + WETM2_ADJ
         WETM2_ADJ = 0.0
         moment3_concb( IMODE ) = moment3_concb( IMODE ) + WETM3_ADJ
         WETM3_ADJ = 0.0

         !------
         ! fwd code:
         ! CINORG( KSO4,IMODE ) = AEROSPC_CONC( ASO4_IDX, IMODE )
         ! CINORG( KNH4,IMODE ) = AEROSPC_CONC( ANH4_IDX, IMODE )
         ! CINORG( KNO3,IMODE ) = AEROSPC_CONC( ANO3_IDX, IMODE )
         ! CINORG( KNA, IMODE ) = AEROSPC_CONC( ANA_IDX, IMODE )
         ! CINORG( KCL, IMODE ) = AEROSPC_CONC( ACL_IDX, IMODE )
         ! adj code:
         aerospc_concb( ACL_IDX, IMODE ) = aerospc_concb( ACL_IDX, IMODE )
     &          + CINORG_ADJ ( KCL, IMODE )
         CINORG_ADJ ( KCL, IMODE ) = 0.0
         aerospc_concb( ANA_IDX, IMODE ) = aerospc_concb( ANA_IDX, IMODE )
     &          + CINORG_ADJ ( KNA, IMODE )
         CINORG_ADJ ( KNA, IMODE ) = 0.0
         aerospc_concb( ANO3_IDX, IMODE ) = aerospc_concb( ANO3_IDX, IMODE )
     &          + CINORG_ADJ ( KNO3, IMODE )
         CINORG_ADJ ( KNO3, IMODE ) = 0.0
         aerospc_concb( ANH4_IDX, IMODE ) = aerospc_concb( ANH4_IDX, IMODE )
     &          + CINORG_ADJ ( KNH4, IMODE )
         CINORG_ADJ ( KNH4, IMODE ) = 0.0
         aerospc_concb( ASO4_IDX, IMODE ) = aerospc_concb( ASO4_IDX, IMODE )
     &          + CINORG_ADJ ( KSO4, IMODE )
         CINORG_ADJ ( KSO4, IMODE ) = 0.0

      END DO

      IMODE = 3

      !------
      ! fwd code:
      ! DMASS = 0.0
      ! adj code:
      DMASS_ADJ = 0.0

      !------
      ! fwd code:
      ! OMEGA( 1 ) = GRFAC2( 1 ) / ( GRFAC2( 1 ) + GRFAC2( 2 ) ) ! partitioning cof
      ! OMEGA( 2 ) = 1.0D0 - OMEGA( 1 )
      ! adj code:
      OMEGA_ADJ ( 1 ) = OMEGA_ADJ ( 1 ) - OMEGA_ADJ ( 2 )
      OMEGA_ADJ ( 2 ) = 0.0
      GRFAC2_ADJ ( 1 ) = GRFAC2_ADJ ( 1 ) + ( GRFAC2( 2 ) / ( GRFAC2(1)
     &                 + GRFAC2(2) ) ** 2.0 ) * OMEGA_ADJ( 1 )
      GRFAC2_ADJ ( 2 ) = GRFAC2_ADJ ( 2 ) - ( GRFAC2( 1 ) / ( ( 
     &      GRFAC2( 1 ) + GRFAC2( 2 ) ) ** 2.0 ) ) * OMEGA_ADJ ( 1 )
      OMEGA_ADJ ( 1 ) = 0.0

      ! Checkpointing
      DVOLINORG = DVOLINORG_CHECK

      IF ( DVOLINORG( KCL ) .LT. 0.0D0 ) THEN
         !------
         ! fwd code:
         ! DVOLMAX = -( DBLE( AEROSPC_CONC( ACL_IDX,1 ) )
         ! &       + DBLE( AEROSPC_CONC( ACL_IDX,2 ) ) ) / DF( KCL )
         ! DVOLINORG( KCL ) = MAX( DVOLINORG( KCL ), DVOLMAX )
         ! adj code:
         IF ( DVOLINORG( KCL ) < DVOLMAX ) THEN
            DVOLMAX_ADJ = DVOLMAX_ADJ + DVOLINORG_ADJ ( KCL )
            DVOLINORG_ADJ ( KCL ) = 0.0
         END IF
         aerospc_concb( ACL_IDX, 1 ) = aerospc_concb( ACL_IDX, 1 )
     &          - ( 1.0 / DF( KCL ) ) * DVOLMAX_ADJ
         aerospc_concb( ACL_IDX, 2 ) = aerospc_concb( ACL_IDX, 2 ) 
     &          - ( 1.0 / DF( KCL ) ) * DVOLMAX_ADJ
         DVOLMAX_ADJ = 0.0
      END IF

      IF ( DVOLINORG( KNO3 ) .LT. 0.0D0 ) THEN
         !------
         ! fwd code:
         ! DVOLMAX = -( DBLE( AEROSPC_CONC( ANO3_IDX,1 ) )
         ! &         + DBLE( AEROSPC_CONC( ANO3_IDX,2 ) ) ) / DF( KNO3 )
         ! DVOLINORG( KNO3 ) = MAX( DVOLINORG( KNO3 ), DVOLMAX)
         ! adj code:
         IF ( DVOLINORG( KNO3 ) < DVOLMAX ) THEN
            DVOLMAX_ADJ = DVOLMAX_ADJ + DVOLINORG_ADJ ( KNO3 )
            DVOLINORG_ADJ ( KNO3 ) = 0.0
         END IF
         aerospc_concb( ANO3_IDX, 1 ) = aerospc_concb( ANO3_IDX, 1 )
     &          - ( 1.0 / DF( KNO3 ) ) * DVOLMAX_ADJ
         aerospc_concb( ANO3_IDX, 2 ) = aerospc_concb( ANO3_IDX, 2 )
     &          - ( 1.0 / DF( KNO3 ) ) * DVOLMAX_ADJ
         DVOLMAX_ADJ = 0.0
      END IF

      IF ( DVOLINORG( KNH4 ) .LT. 0.0D0 ) THEN
         !------
         ! fwd code:
         ! DVOLMAX = -( DBLE( AEROSPC_CONC( ANH4_IDX,1 ) )
         ! &         + DBLE( AEROSPC_CONC( ANH4_IDX,2 ) ) ) / DF( KNH4 )
         ! DVOLINORG( KNH4 ) = MAX( DVOLINORG( KNH4 ), DVOLMAX)
         ! adj code:
         IF ( DVOLINORG( KNH4 ) < DVOLMAX ) THEN
            DVOLMAX_ADJ = DVOLMAX_ADJ + DVOLINORG_ADJ ( KNH4 )
            DVOLINORG_ADJ ( KNH4 ) = 0.0
         END IF
         aerospc_concb( ANH4_IDX, 1 ) = aerospc_concb( ANH4_IDX, 1 )
     &          - ( 1.0 / DF( KNH4 ) ) * DVOLMAX_ADJ
         aerospc_concb( ANH4_IDX, 2 ) = aerospc_concb( ANH4_IDX, 2 )
     &          - ( 1.0 / DF( KNH4 ) ) * DVOLMAX_ADJ
         DVOLMAX_ADJ = 0.0
      END IF

      IF ( TRUSTCL_CHECK ) THEN
         !------
         ! fwd code:
         ! DVOLINORG( KCL )  = GCLR8  * 1.0D-6 / DBLE( PRECURSOR_MW( HCL_IDX ) ) - GAS( 3 ) 
         ! adj code:
         GCLR8_ADJ = GCLR8_ADJ + ( 1.0D-6 / DBLE( PRECURSOR_MW( HCL_IDX ) ) )
     &             * DVOLINORG_ADJ ( KCL )
         GAS_ADJ ( 3 ) = GAS_ADJ ( 3 ) - DVOLINORG_ADJ ( KCL )
         DVOLINORG_ADJ ( KCL ) = 0.0
      ELSE  
         DVOLINORG_ADJ ( KCL ) = 0.0 
      END IF

      !------
      ! fwd code:
      ! DVOLINORG( KNH4 ) = GNH3R8 * 1.0D-6 / DBLE( PRECURSOR_MW( NH3_IDX ) ) - GAS( 1 )  
      ! DVOLINORG( KNO3 ) = GNO3R8 * 1.0D-6 / DBLE( PRECURSOR_MW( HNO3_IDX ) ) - GAS( 2 )
      ! adj code:
      GNO3R8_ADJ = GNO3R8_ADJ + ( 1.0D-6 / DBLE( PRECURSOR_MW( HNO3_IDX ) ) )
     &           * DVOLINORG_ADJ ( KNO3 )
      GAS_ADJ ( 2 ) = GAS_ADJ ( 2 ) - DVOLINORG_ADJ ( KNO3 )
      DVOLINORG_ADJ ( KNO3 ) = 0.0
      GNH3R8_ADJ = GNH3R8_ADJ + ( 1.0D-6 / DBLE( PRECURSOR_MW( NH3_IDX ) ) )
     &           * DVOLINORG_ADJ ( KNH4 ) 
      GAS_ADJ ( 1 ) = GAS_ADJ ( 1 ) - DVOLINORG_ADJ ( KNH4 )
      DVOLINORG_ADJ ( KNH4 ) = 0.0

      !------
      ! fwd code:
      ! PRECURSOR_CONC( NH3_IDX )  = GAS( 1 ) * PRECURSOR_MW( NH3_IDX ) * 1.0E6
      ! PRECURSOR_CONC( HNO3_IDX ) = GAS( 2 ) * PRECURSOR_MW( HNO3_IDX ) * 1.0E6
      ! PRECURSOR_CONC( HCL_IDX )  = GAS( 3 ) * PRECURSOR_MW( HCL_IDX ) * 1.0E6
      ! adj code:
      GAS_ADJ ( 3 ) = GAS_ADJ ( 3 ) + ( PRECURSOR_MW( HCL_IDX ) * 1.0E6 )
     &              * precursor_concb( HCL_IDX )
      precursor_concb( HCL_IDX ) = 0.0
      GAS_ADJ ( 2 ) = GAS_ADJ ( 2 ) + ( PRECURSOR_MW( HNO3_IDX ) * 1.0E6 )
     &              * precursor_concb( HNO3_IDX )
      precursor_concb( HNO3_IDX ) = 0.0
      GAS_ADJ ( 1 ) = GAS_ADJ ( 1 ) + ( PRECURSOR_MW( NH3_IDX ) * 1.0E6 )
     &              * precursor_concb( NH3_IDX )
      precursor_concb( NH3_IDX ) = 0.0

      ! Checkpointing
      WI = WI_CHECK( : )
      GAS = GAS_CHECK( : )
      WT = WT_CHECK( : )
      AERLIQ = AERLIQ_CHECK( : )
      AERSLD = AERSLD_CHECK( : )
      OTHER = OTHER_CHECK( : )

      CNTRL( 1 ) = 0.0D0
      CNTRL( 2 ) = 1.0D0

      !------
      ! fwd code:
      ! CALL ISOROPIA( WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,
      ! &    AERSLD, SCASI, OTHER )
      ! adj code:
      CALL ISOROPIA_B( WI, WI_ADJ, RHI, TEMPI, CNTRL, WT, GAS, 
     &                 GAS_ADJ, AERLIQ, AERLIQ_ADJ, AERSLD, SCASI, OTHER,
     &                 TRUSTISO )

        WI_ADJ = 0.0     !debug, mdt

      !------
      ! fwd code:
      ! WI( 5 ) = PRECURSOR_CONC(HCL_IDX) * ( 1.0D-6 / PRECURSOR_MW(HCL_IDX) ) 
      ! &        + ( AEROSPC_CONC(ACL_IDX,1) + AEROSPC_CONC(ACL_IDX,2) )
      ! &        * ( 1.0D-6 / AEROSPC_MW(ACL_IDX) )
      ! adj code:
      precursor_concb( HCL_IDX ) = precursor_concb( HCL_IDX ) + 
     &         ( 1.0D-6 / PRECURSOR_MW(HCL_IDX) ) * WI_ADJ( 5 )
      aerospc_concb( ACL_IDX, 1 ) = aerospc_concb( ACL_IDX, 1 ) +
     &       ( 1.0D-6 / AEROSPC_MW(ACL_IDX) ) * WI_ADJ ( 5 )
      aerospc_concb( ACL_IDX, 2 ) = aerospc_concb( ACL_IDX, 2 ) +
     &       ( 1.0D-6 / AEROSPC_MW(ACL_IDX) ) * WI_ADJ ( 5 )
      WI_ADJ ( 5 ) = 0.0

      !------
      ! fwd code:
      ! WI( 4 ) = PRECURSOR_CONC( HNO3_IDX ) * ( 1.0D-6 / PRECURSOR_MW( HNO3_IDX ) )
      ! &        + ( AEROSPC_CONC( ANO3_IDX,1 ) + AEROSPC_CONC( ANO3_IDX,2 ) )
      ! &        * ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) )
      ! adj code:
      precursor_concb( HNO3_IDX ) = precursor_concb( HNO3_IDX ) +
     &         ( 1.0D-6 / PRECURSOR_MW( HNO3_IDX ) ) * WI_ADJ ( 4 )
      aerospc_concb( ANO3_IDX, 1 ) = aerospc_concb( ANO3_IDX, 1 ) +
     &       ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) ) * WI_ADJ ( 4 )
      aerospc_concb( ANO3_IDX, 2 ) = aerospc_concb( ANO3_IDX, 2 ) +
     &       ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) ) * WI_ADJ ( 4 )
      WI_ADJ( 4 ) = 0.0

      !------
      ! fwd code:
      ! WI( 3 ) = PRECURSOR_CONC( NH3_IDX ) * ( 1.0D-6 / PRECURSOR_MW( NH3_IDX ) )
      ! &        + ( AEROSPC_CONC( ANH4_IDX,1 ) + AEROSPC_CONC( ANH4_IDX,2 ) )
      ! &        * ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) )
      ! adj code:
      precursor_concb( NH3_IDX ) = precursor_concb( NH3_IDX ) +
     &         ( 1.0D-6 / PRECURSOR_MW( NH3_IDX ) ) * WI_ADJ ( 3 )
      aerospc_concb( ANH4_IDX, 1 ) = aerospc_concb( ANH4_IDX, 1 ) +
     &       ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) ) * WI_ADJ ( 3 )
      aerospc_concb( ANH4_IDX, 2 ) = aerospc_concb( ANH4_IDX, 2 ) +
     &       ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) ) * WI_ADJ ( 3 )
      WI_ADJ ( 3 ) = 0.0

      !------
      ! fwd code:
      ! WI( 2 ) = SO4 * ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) )
      ! adj code:
      SO4_ADJ = ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) ) * 
     &          WI_ADJ ( 2 )
      WI_ADJ ( 2 ) = 0.0

      !------
      ! fwd code:
      ! WI( 1 ) = ( AEROSPC_CONC( ANA_IDX,1 ) + AEROSPC_CONC( ANA_IDX,2 ) )
      ! &        * ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) )
      ! adj code:
      aerospc_concb( ANA_IDX, 1 ) = aerospc_concb( ANA_IDX, 1 ) +
     &       ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) ) * WI_ADJ ( 1 )
      aerospc_concb( ANA_IDX, 2 ) = aerospc_concb( ANA_IDX, 2 ) +
     &       ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) ) * WI_ADJ ( 1 )
      WI_ADJ ( 1 ) = 0.0

      !------
      ! fwd code:
      ! SO4 = AEROSPC_CONC( ASO4_IDX,1 ) + AEROSPC_CONC( ASO4_IDX,2 )
      ! SO4 = SO4 + ( DMDT_SO4 + CONDSO4( 1 ) + CONDSO4( 2 ) )
      ! &              * DELT * H2SO4RATM1
      ! adj code:
      DMDT_SO4_ADJ = DMDT_SO4_ADJ + DELT * H2SO4RATM1 * SO4_ADJ
      CONDSO4_ADJ( 1 ) = CONDSO4_ADJ ( 1 ) + DELT * H2SO4RATM1 * SO4_ADJ
      CONDSO4_ADJ( 2 ) = CONDSO4_ADJ ( 2 ) + DELT * H2SO4RATM1 * SO4_ADJ
      aerospc_concb( ASO4_IDX, 1 ) = aerospc_concb( ASO4_IDX, 1 ) + SO4_ADJ
      aerospc_concb( ASO4_IDX, 2 ) = aerospc_concb( ASO4_IDX, 2 ) + SO4_ADJ
      SO4_ADJ = 0.0

      !------
      ! fwd code:
      ! GNH3R8  = PRECURSOR_CONC( NH3_IDX )
      ! GNO3R8  = PRECURSOR_CONC( HNO3_IDX )
      ! GCLR8   = PRECURSOR_CONC( HCL_IDX )
      ! adj code:
      precursor_concb( HCL_IDX ) = precursor_concb( HCL_IDX ) +
     &         GCLR8_ADJ
      GCLR8_ADJ = 0.0
      precursor_concb( HNO3_IDX ) = precursor_concb( HNO3_IDX ) +
     &         GNO3R8_ADJ
      GNO3R8_ADJ = 0.0
      precursor_concb( NH3_IDX ) = precursor_concb( NH3_IDX ) +
     &         GNH3R8_adj
      GNH3R8_ADJ = 0.0

      IF ( HYBRID ) THEN

         GAS_ADJ = 0.0
         AERLIQ_ADJ = 0.0
         HPLUS_ADJ = 0.0
         AM1_ADJ = 0.0
         AM2_ADJ = 0.0
         FCONC_SO4_ADJ = 0.0

         DO B = COUNTER, 1, -1

            ! Checkpointing
            aeromode_diam = aeromode_diam_check( B, : )
            aeromode_sdev = aeromode_sdev_check( B, : )
            aeromode_mass = aeromode_mass_check( B, : )
            aerospc_conc  = aerospc_conc_check ( B, :, : )
            moment0_conc = moment0_conc_check( B, : )
            moment2_conc = moment2_conc_check( B, : )
            moment3_conc = moment3_conc_check( B, : )

            !------
            ! fwd code:
            ! CALL GETPAR( M3_WET_FLAG, LIMIT_Sg )
            ! adj code:
            CALL GETPAR_B( M3_WET_FLAG, LIMIT_Sg )

            !------
            ! fwd code:
            ! AEROSPC_CONC( ANH4_IDX,IMODE ) = CFINAL( KNH4,IMODE )
            ! AEROSPC_CONC( ANO3_IDX,IMODE ) = CFINAL( KNO3,IMODE )
            ! AEROSPC_CONC( ACL_IDX,IMODE )  = CFINAL( KCL, IMODE )
            ! AEROSPC_CONC( ASO4_IDX,IMODE ) = SO4
            ! AEROSPC_CONC( AH2O_IDX,IMODE ) = H2O
            ! MOMENT2_CONC( IMODE ) = WETM2
            ! adj code:
            WETM2_ADJ = WETM2_ADJ + moment2_concb( IMODE )
            moment2_concb( IMODE ) = 0.0
            H2O_ADJ = H2O_ADJ + aerospc_concb( AH2O_IDX, IMODE )
            aerospc_concb( AH2O_IDX, IMODE ) = 0.0
            SO4_ADJ = aerospc_concb( ASO4_IDX, IMODE )
            aerospc_concb( ASO4_IDX, IMODE ) = 0.0
            CFINAL_ADJ ( KCL, IMODE ) = CFINAL_ADJ ( KCL, IMODE ) +
     &            aerospc_concb( ACL_IDX, IMODE )
            aerospc_concb( ACL_IDX, IMODE ) = 0.0
            CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) +
     &            aerospc_concb( ANO3_IDX, IMODE )
            aerospc_concb( ANO3_IDX, IMODE ) = 0.0
            CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) +
     &            aerospc_concb( ANH4_IDX, IMODE )
            aerospc_concb( ANH4_IDX, IMODE ) = 0.0

            !------
            ! fwd code:
            ! PRECURSOR_CONC( HCL_IDX ) = GCLR8 + ( CINORG( KCL,IMODE )
            ! &                        - CFINAL( KCL,IMODE) ) * HCLRAT
            ! adj code:
            GCLR8_ADJ = GCLR8_ADJ + precursor_concb( HCL_IDX )
            CINORG_ADJ ( KCL, IMODE ) = CINORG_ADJ( KCL, IMODE ) +
     &            HCLRAT * precursor_concb( HCL_IDX )
            CFINAL_ADJ ( KCL, IMODE ) = CFINAL_ADJ ( KCL, IMODE ) -
     &            HCLRAT * precursor_concb( HCL_IDX )
            precursor_concb( HCL_IDX ) = 0.0

            !------
            ! fwd code:
            ! PRECURSOR_CONC( HNO3_IDX ) = GNO3R8 + ( CINORG( KNO3,IMODE )
            ! &                          - CFINAL( KNO3,IMODE) ) * HNO3RAT
            ! adj code:
            GNO3R8_ADJ = GNO3R8_ADJ + precursor_concb( HNO3_IDX )
            CINORG_ADJ ( KNO3, IMODE ) = CINORG_ADJ ( KNO3, IMODE ) +
     &            HNO3RAT * precursor_concb( HNO3_IDX )
            CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) - 
     &            HNO3RAT * precursor_concb( HNO3_IDX )
            precursor_concb( HNO3_IDX ) = 0.0

            !------
            ! fwd code: 
            ! PRECURSOR_CONC( NH3_IDX ) = GNH3R8 + ( CINORG( KNH4,IMODE )
            ! &                         - CFINAL( KNH4,IMODE ) ) * NH3RAT
            ! adj code:
            GNH3R8_ADJ = GNH3R8_ADJ + precursor_concb( NH3_IDX )
            CINORG_ADJ ( KNH4, IMODE ) = CINORG_ADJ ( KNH4, IMODE ) +
     &            NH3RAT * precursor_concb( NH3_IDX )
            CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) -
     &            NH3RAT * precursor_concb( NH3_IDX )
            precursor_concb( NH3_IDX ) = 0.0

            ! Checkpointing
            WETM3 = WETM3_CHECK2( B )
            DRYM3 = DRYM3_CHECK2( B )
            DRYM2 = DRYM2_CHECK2( B )
            Y = Y_CHECK1( B )

            !------
            ! fwd code:
            ! WETM3 = DRYM3 + M3SOA + ( 1.0D-9 * F6DPI / AEROSPC( AH2O_IDX )%DENSITY ) * H2O  
            ! DRYM2 = MAX( DBLE( AEROMODE( IMODE )%MIN_M2 ), Y )
            ! WETM2 = DRYM2 * ( WETM3 / DRYM3 ) ** D_TWOTHIRDS
            ! adj code:
            TEMP2 = WETM3 / DRYM3
!            IF ( TEMP2 .LE. 0.0 ) THEN
!               TEMP2B5 = 0.0
!            ELSE
               TEMP2B5 = DRYM2 * D_TWOTHIRDS * TEMP2 ** ( D_TWOTHIRDS - 1 )
     &                 * WETM2_ADJ / DRYM3 
!            END IF
            DRYM2_ADJ = TEMP2 ** D_TWOTHIRDS * WETM2_ADJ 
            WETM3_ADJ = TEMP2B5
            DRYM3_ADJ = -( TEMP2 * TEMP2B5 )
            IF ( Y > DBLE ( AEROMODE( IMODE )%MIN_M2 ) ) THEN
               Y_ADJ = DRYM2_ADJ 
               DRYM2_ADJ = 0.0
            END IF
            DRYM3_ADJ = DRYM3_ADJ + WETM3_ADJ
            M3SOA_ADJ = M3SOA_ADJ + WETM3_ADJ
            H2O_ADJ = H2O_ADJ + ( 1.0D-9 * F6DPI / AEROSPC( 
     &                AH2O_IDX )%DENSITY ) * WETM3_ADJ
            WETM3_ADJ = 0.0

            !------
            ! fwd code:
            ! DRYM3 = ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY ) * SO4
            ! &     + ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY )
            ! &                        * CFINAL( KNH4,IMODE )
            ! &     + ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY )
            ! &                        * CFINAL( KNO3,IMODE )
            ! &     + ( 1.0D-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY )
            ! &                        * CFINAL( KCL,IMODE )
            ! &     + ( 1.0D-9 * F6DPI / AEROSPC( ANA_IDX )%DENSITY )
            ! &                        * CINORG( KNA,IMODE )
            ! &     + M3OTHR
            ! adj code:
            SO4_ADJ = SO4_ADJ + ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX
     &                )%DENSITY ) * DRYM3_ADJ
            CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) +
     &            ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY ) *
     &            DRYM3_ADJ
            CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) +
     &            ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY ) *
     &            DRYM3_ADJ
            CFINAL_ADJ ( KCL, IMODE ) = CFINAL_ADJ ( KCL, IMODE ) +
     &            ( 1.0D-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY ) *
     &            DRYM3_ADJ
            CINORG_ADJ ( KNA, IMODE ) = CINORG_ADJ ( KNA, IMODE ) +
     &            ( 1.0D-9 * F6DPI / AEROSPC( ANA_IDX )%DENSITY ) *
     &            DRYM3_ADJ
            M3OTHR_ADJ = M3OTHR_ADJ + DRYM3_ADJ
            DRYM3_ADJ = 0.0

            ! Checkpointing
            DDRYM2DT = DDRYM2DT_CHECK( B )

            IF ( DDRYM2DT .LT. 0.0D0 ) THEN

               ! Checkpointing
               LOSS = LOSS_CHECK( B )
               DRYM20 = DRYM20_CHECK( B )

               !------
               ! fwd code:
               ! LOSS = DDRYM2DT / DRYM20
               ! Y = DRYM20 * EXP( LOSS * TSTEP )
               ! adj code:
               DRYM20_ADJ = DRYM20_ADJ + EXP ( LOSS * TSTEP ) * Y_ADJ
               LOSS_ADJ = LOSS_ADJ + DRYM20 * TSTEP * EXP ( LOSS * 
     &                    TSTEP ) * Y_ADJ
               Y_ADJ = 0.0
               DDRYM2DT_ADJ = DDRYM2DT_ADJ + ( 1.0 / DRYM20 ) * LOSS_ADJ
               DRYM20_ADJ = DRYM20_ADJ - ( DDRYM2DT / ( DRYM20 ** 2.0 ) )
     &                    * LOSS_ADJ
               LOSS_ADJ = 0.0
            ELSE
               !------
               ! fwd code:
               ! Y = DRYM20 + DDRYM2DT * TSTEP
               ! adj code:
               DRYM20_ADJ = DRYM20_ADJ + Y_ADJ 
               DDRYM2DT_ADJ = DDRYM2DT_ADJ + TSTEP * Y_ADJ
               Y_ADJ = 0.0
            ENDIF

            ! Checkpointing
            GRFAC2 = GRFAC2_CHECK( B, : )
            GRFAC1 = GRFAC1_CHECK( B, : )
            DDRYM3DT = DDRYM3DT_CHECK( B )

            !------
            ! fwd code:
            ! DDRYM2DT = D_TWOTHIRDS * GRFAC1( IMODE ) / GRFAC2( IMODE ) * DDRYM3DT
            ! adj code:
            DDRYM3DT_ADJ = DDRYM3DT_ADJ + ( D_TWOTHIRDS * GRFAC1( IMODE ) 
     &                   / GRFAC2( IMODE ) ) * DDRYM2DT_ADJ 
            GRFAC1_ADJ ( IMODE ) = GRFAC1_ADJ ( IMODE ) + D_TWOTHIRDS * 
     &            ( 1.0 / GRFAC2( IMODE ) ) * DDRYM3DT * DDRYM2DT_ADJ
            GRFAC2_ADJ ( IMODE ) = GRFAC2_ADJ ( IMODE ) - D_TWOTHIRDS *
     &            ( GRFAC1( IMODE ) / ( GRFAC2( IMODE ) ** 2.0 ) ) * 
     &            DDRYM3DT * DDRYM2DT_ADJ
            DDRYM2DT_ADJ = 0.0

            !------
            ! fwd code:
            ! DDRYM3DT = (
            ! &            ( CFINAL( KNH4,IMODE ) - CINORG( KNH4,IMODE ) )
            !              * ( 1.0D-9 * F6DPI / AEROSPC( ANH4_IDX )%DENSITY)
            !          + ( CFINAL( KNO3,IMODE ) - CINORG( KNO3,IMODE ) )
            !              * ( 1.0D-9 * F6DPI / AEROSPC( ANO3_IDX )%DENSITY)
            !          + ( CFINAL( KCL, IMODE ) - CINORG( KCL,IMODE ) )
            !              * ( 1.0D-9 * F6DPI / AEROSPC( ACL_IDX )%DENSITY)
            !          + ( SO4                  - CINORG( KSO4,IMODE ) )
            !                  * ( 1.0D-9 * F6DPI / AEROSPC( ASO4_IDX )%DENSITY) ) / TSTEP
            ! adj code:
            CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) +
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC( 
     &            ANH4_IDX )%DENSITY) * DDRYM3DT_ADJ
            CINORG_ADJ ( KNH4, IMODE ) = CINORG_ADJ ( KNH4, IMODE ) - 
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC(
     &            ANH4_IDX )%DENSITY) * DDRYM3DT_ADJ
            CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) +
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC( 
     &            ANO3_IDX )%DENSITY) * DDRYM3DT_ADJ
            CINORG_ADJ ( KNO3, IMODE ) = CINORG_ADJ ( KNO3, IMODE ) - 
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC(
     &            ANO3_IDX )%DENSITY) * DDRYM3DT_ADJ
            CFINAL_ADJ ( KCL , IMODE ) = CFINAL_ADJ ( KCL , IMODE ) +
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC( 
     &             ACL_IDX )%DENSITY) * DDRYM3DT_ADJ
            CINORG_ADJ ( KCL , IMODE ) = CINORG_ADJ ( KCL , IMODE ) - 
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC(
     &             ACL_IDX )%DENSITY) * DDRYM3DT_ADJ
            SO4_ADJ = SO4_ADJ + 
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC( 
     &            ASO4_IDX )%DENSITY) * DDRYM3DT_ADJ
            CINORG_ADJ ( KSO4, IMODE ) = CINORG_ADJ ( KSO4, IMODE ) - 
     &            ( 1.0 / TSTEP ) * ( 1.0D-9 * F6DPI / AEROSPC(
     &            ASO4_IDX )%DENSITY) * DDRYM3DT_ADJ
     
            !------
            ! fwd code:
            ! H2O = H2O_NEW * DFH2OR8
            ! adj code:
            H2O_NEW_ADJ = H2O_NEW_ADJ + DFH2OR8 * H2O_ADJ
            H2O_ADJ = 0.0

            ! Checkpointing
            WI = WI_CHECK1( CALC_H2O_STEP, : )
            H2O_NEW = H2O_NEW_CHECK( CALC_H2O_STEP )

            !------
            ! fwd code:
            ! CALL CALC_H2O( WI, RHI, TEMPI, H2O_NEW )
            ! adj code:
            CALL CALC_H2O_ADJ( WI, WI_ADJ, RHI, TEMPI, H2O_NEW, H2O_NEW_ADJ )

            CALC_H2O_STEP = CALC_H2O_STEP - 1

            !------
            ! fwd code:
            ! WI( 1 ) = CINORG( KNA, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) )
            ! WI( 2 ) =                  SO4 * ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) )
            ! WI( 3 ) = CFINAL( KNH4,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) )
            ! WI( 4 ) = CFINAL( KNO3,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) )
            ! WI( 5 ) = CFINAL( KCL, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ACL_IDX ) )
            ! adj code:
            CFINAL_ADJ ( KCL, IMODE ) = CFINAL_ADJ ( KCL, IMODE ) +
     &            ( 1.0D-6 / AEROSPC_MW( ACL_IDX ) ) * WI_ADJ ( 5 )
            WI_ADJ ( 5 ) = 0.0
            CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) +
     &            ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) ) * WI_ADJ ( 4 )
            WI_ADJ ( 4 ) = 0.0
            CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) +
     &            ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) ) * WI_ADJ ( 3 )
            WI_ADJ ( 3 ) = 0.0
            SO4_ADJ = SO4_ADJ + ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) ) *
     &                WI_ADJ ( 2 )
            WI_ADJ ( 2 ) = 0.0
            CINORG_ADJ ( KNA, IMODE ) = CINORG_ADJ ( KNA, IMODE ) +
     &            ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) ) * WI_ADJ ( 1 )

            !------
            ! fwd code:
            ! HPLUS( IMODE ) = 2.0D0 * SO4          / AEROSPC_MW( ASO4_IDX )
            !              + CFINAL( KNO3,IMODE ) / AEROSPC_MW( ANO3_IDX )
            !              + CFINAL( KCL, IMODE ) / AEROSPC_MW( ACL_IDX )
            !              - CFINAL( KNH4,IMODE ) / AEROSPC_MW( ANH4_IDX )
            !              - CINORG( KNA, IMODE ) / AEROSPC_MW( ANA_IDX )
            ! adj code:
            SO4_ADJ = SO4_ADJ + 2.0 / AEROSPC_MW( ASO4_IDX ) *
     &                HPLUS_ADJ ( IMODE )
            CFINAL_ADJ ( KNO3, IMODE ) = CFINAL_ADJ ( KNO3, IMODE ) + 
     &            ( 1.0 / AEROSPC_MW( ANO3_IDX ) ) * HPLUS_ADJ ( IMODE )
            CFINAL_ADJ ( KCL, IMODE ) = CFINAL_ADJ ( KCL, IMODE ) +
     &            ( 1.0 / AEROSPC_MW( ACL_IDX ) ) * HPLUS_ADJ ( IMODE )
            CFINAL_ADJ ( KNH4, IMODE ) = CFINAL_ADJ ( KNH4, IMODE ) - 
     &            ( 1.0 / AEROSPC_MW( ANH4_IDX ) ) * HPLUS_ADJ ( IMODE )
            CINORG_ADJ ( KNA, IMODE ) = CINORG_ADJ ( KNA, IMODE ) - 
     &            ( 1.0 / AEROSPC_MW( ANA_IDX ) ) * HPLUS_ADJ ( IMODE )
            HPLUS_ADJ ( IMODE ) = 0.0
            J_ADJ = 0.0

            DO ISP = NVOLINORG, 1, -1
       
               ! Checkpointing
               CINORG = CINORG_CHECK( B, :, : )
               J = J_CHECK( B, : )

               !------
               ! fwd code:
               ! CFINAL( ISP,IMODE ) = MAX( 0.0D0,
               !                          CINORG( ISP,IMODE )
               !                          + J( ISP ) * TSTEP * DF( ISP ) )
               ! adj code:
               IF ( CINORG( ISP, IMODE ) + J( ISP ) * TSTEP * DF( ISP ) 
     &              > 0 ) THEN
                  CINORG_ADJ ( ISP, IMODE ) = CINORG_ADJ ( ISP, IMODE ) 
     &                  + CFINAL_ADJ ( ISP, IMODE )
                  J_ADJ ( ISP ) = J_ADJ ( ISP ) + TSTEP * DF( ISP ) *
     &                  CFINAL_ADJ ( ISP, IMODE )
                  CFINAL_ADJ ( ISP, IMODE ) = 0.0
               ELSE 
                  CFINAL_ADJ ( ISP, IMODE ) = 0.0
               END IF

            END DO

            ! Checkpointing
            GCLR8 = GCLR8_CHECK( B )
            GNO3R8 = GNO3R8_CHECK( B ) 
            GNH3R8 = GNH3R8_CHECK( B )
            J = J_CHECK2( B, : )

            IF ( J( KCL ) * TSTEP * DF(KCL) * HCLRAT .GT. GCLR8 ) THEN
               !------
               ! fwd code:
               ! J( KCL ) = GCLR8 / ( TSTEP * DF( KCL ) * HCLRAT )
               ! adj code:
               GCLR8_ADJ = GCLR8_ADJ + ( 1.0 / ( TSTEP * DF( KCL ) * 
     &                     HCLRAT ) ) * J_ADJ ( KCL )
               J_ADJ ( KCL ) = 0.0
            END IF
  
            IF ( J( KNO3 ) * TSTEP * DF( KNO3 ) * HNO3RAT .GT. GNO3R8 ) THEN
               !------
               ! fwd code:
               ! J( KNO3 ) = GNO3R8 / ( TSTEP * DF( KNO3 ) * HNO3RAT )
               ! adj code:
               GNO3R8_ADJ = GNO3R8_ADJ + ( 1.0 / ( TSTEP * DF( KNO3 ) * 
     &                      HNO3RAT ) ) * J_ADJ ( KNO3 )
               J_ADJ ( KNO3 ) = 0.0
            ENDIF

            IF ( J( KNH4 ) * TSTEP * DF( KNH4 ) * NH3RAT .GT. GNH3R8 ) THEN
               !------
               ! fwd code:
               ! J( KNH4 ) = GNH3R8 / ( TSTEP * DF( KNH4 ) * NH3RAT )
               ! adj code:
               GNH3R8_ADJ = GNH3R8_ADJ + ( 1.0 / ( TSTEP * DF( KNH4 ) * 
     &                      NH3RAT ) ) * J_ADJ ( KNH4 ) 
               J_ADJ ( KNH4 ) = 0.0
            END IF

            IF ( TrustIso ) THEN
               ! Checkpointing
               GAS = GAS_CHECK1( B, : )
               GRFAC2 = GRFAC2_CHECK1( B, : )
               HPLUS = HPLUS_CHECK( B , : )
               RATE = RATE_CHECK( B )
               J = J_CHECK1( B, : )

               !------
               ! fwd code:
               ! CALL COMPUTE_FLUX( NVOLINORG, GNH3R8, GNO3R8, GCLR8, KNH4,
               ! KNO3, KCL, GAS( 1:3 ), GRFAC2( IMODE ), AERLIQ( 1 ), RATE, J )
               ! adj code:
               CALL COMPUTE_FLUX_ADJ ( NVOLINORG, GNH3R8, GNH3R8_ADJ, 
     &                              GNO3R8, GNO3R8_ADJ, GCLR8, GCLR8_ADJ,
     &                              KNH4, KNO3, KCL, GAS(1:3), 
     &                              GAS_ADJ(1:3), GRFAC2( IMODE ),
     &                              GRFAC2_ADJ( IMODE ), AERLIQ( 1 ), 
     &                              AERLIQ_ADJ( 1 ), RATE,
     &                              RATE_ADJ, J, J_ADJ )

           ELSE

               !------
               ! fwd code:
               ! J ( : ) = 0.0
               ! adj code:
               J_ADJ ( : ) = 0.0

            END IF

            !------
            ! fwd code:
            ! DVOLINORG( KNH4 ) = GNH3R8 * ( 1.0D-6 / PRECURSOR_MW( NH3_IDX ) ) - GAS( 1 )  
            ! DVOLINORG( KNO3 ) = GNO3R8 * ( 1.0D-6 / PRECURSOR_MW( HNO3_IDX ) ) - GAS( 2 ) 
            ! DVOLINORG( KCL )  = GCLR8  * ( 1.0D-6 / PRECURSOR_MW( HCL_IDX ) ) - GAS( 3 )
            ! adj code:
            GCLR8_ADJ = GCLR8_ADJ + ( 1.0D-6 / PRECURSOR_MW( HCL_IDX ) )
     &                * DVOLINORG_ADJ ( KCL )
            GAS_ADJ ( 3 ) = GAS_ADJ ( 3 ) - DVOLINORG_ADJ ( KCL )
            DVOLINORG_ADJ ( KCL ) = 0.0
            GNO3R8_ADJ = GNO3R8_ADJ + ( 1.0D-6 / PRECURSOR_MW( HNO3_IDX ) )
     &                 * DVOLINORG_ADJ ( KNO3 ) 
            GAS_ADJ ( 2 ) = GAS_ADJ ( 2 ) - DVOLINORG_ADJ ( KNO3 )
            DVOLINORG_ADJ ( KNO3 ) = 0.0
            GNH3R8_ADJ = GNH3R8_ADJ + ( 1.0D-6 / PRECURSOR_MW( NH3_IDX ) )
     &                 * DVOLINORG_ADJ ( KNH4 )
            GAS_ADJ ( 1 ) = GAS_ADJ ( 1 ) - DVOLINORG_ADJ ( KNH4 )
            DVOLINORG_ADJ ( KNH4 ) = 0.0

            IF ( TrustIso ) THEN

               ! Checkpointing
               WI = WI_CHECK2( B, : )
               WT = WT_CHECK1( B, : )
               CNTRL = CNTRL_CHECK1( B, : )
               GAS = GAS_CHECK2( B, : )
               AERLIQ = AERLIQ_CHECK1( B, : )
               AERSLD = AERSLD_CHECK1( B, : )
               OTHER = OTHER_CHECK1( B, : )

               !------
               ! fwd code:
               ! CALL ISOROPIA( WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,  
               !                  AERSLD, SCASI, OTHER )
               ! adj code:
               CALL ISOROPIA_B  ( WI, WI_ADJ, RHI, TEMPI, CNTRL, WT,
     &                            GAS, GAS_ADJ, AERLIQ, AERLIQ_ADJ,
     &                            AERSLD, SCASI, OTHER, TRUSTISO )

               !------
               ! fwd code:
               ! WI( 1 ) = CINORG( KNA, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) )
               ! WI( 2 ) =                  SO4 * ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) )
               ! WI( 3 ) = CINORG( KNH4,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) )
               ! WI( 4 ) = CINORG( KNO3,IMODE ) * ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) )
               ! WI( 5 ) = CINORG( KCL, IMODE ) * ( 1.0D-6 / AEROSPC_MW( ACL_IDX ) )
               ! adj code:
               CINORG_ADJ ( KCL, IMODE ) = CINORG_ADJ ( KCL, IMODE ) +
     &               ( 1.0D-6 / AEROSPC_MW( ACL_IDX ) ) * WI_ADJ ( 5 )
               WI_ADJ ( 5 ) = 0.0
               CINORG_ADJ ( KNO3, IMODE ) = CINORG_ADJ ( KNO3, IMODE ) +
     &               ( 1.0D-6 / AEROSPC_MW( ANO3_IDX ) ) * WI_ADJ ( 4 )
               WI_ADJ ( 4 ) = 0.0
               CINORG_ADJ ( KNH4, IMODE ) = CINORG_ADJ ( KNH4, IMODE ) +
     &               ( 1.0D-6 / AEROSPC_MW( ANH4_IDX ) ) * WI_ADJ ( 3 )
               WI_ADJ ( 3 ) = 0.0
               SO4_ADJ = SO4_ADJ + ( 1.0D-6 / AEROSPC_MW( ASO4_IDX ) ) *
     &                 WI_ADJ ( 2 )
               WI_ADJ ( 2 ) = 0.0
               CINORG_ADJ ( KNA, IMODE ) = CINORG_ADJ ( KNA, IMODE ) +
     &               ( 1.0D-6 / AEROSPC_MW( ANA_IDX ) ) * WI_ADJ ( 1 )
               WI_ADJ ( 1 ) = 0.0

               GAS_ADJ = 0.0
               AERLIQ_ADJ = 0.0
 
            END IF

            !------
            ! fwd code:
            ! RATE = CONDSO4( IMODE )
            ! SO4  = CINORG( KSO4,IMODE ) + RATE * TSTEP * H2SO4RATM1
            ! adj code:
            CINORG_ADJ ( KSO4, IMODE ) = CINORG_ADJ ( KSO4, IMODE ) +
     &            SO4_ADJ 
            RATE_ADJ = RATE_ADJ + TSTEP * H2SO4RATM1 * SO4_ADJ
            CONDSO4_ADJ ( IMODE ) = CONDSO4_ADJ ( IMODE ) + RATE_ADJ

            !------
            ! fwd code:
            ! GNO3R8 = PRECURSOR_CONC( HNO3_IDX )
            ! GNH3R8 = PRECURSOR_CONC( NH3_IDX )
            ! GCLR8  = PRECURSOR_CONC( HCL_IDX )
            ! adj code:
            precursor_concb( HCL_IDX ) = precursor_concb( HCL_IDX ) +
     &               GCLR8_ADJ
            GCLR8_ADJ = 0.0
            precursor_concb( NH3_IDX ) = precursor_concb( NH3_IDX ) +
     &               GNH3R8_ADJ
            GNH3R8_ADJ = 0.0
            precursor_concb( HNO3_IDX ) = precursor_concb( HNO3_IDX ) +
     &               GNO3R8_ADJ
            GNO3R8_ADJ = 0.0

            ! Checkpointing
            WETM3 = WETM3_CHECK3( B )
            DRYM3 = DRYM3_CHECK3( B )
            WETM2 = WETM2_CHECK3( B )

            !------
            ! fwd code:
            ! WETM3  = MOMENT3_CONC( IMODE )
            ! WETM2  = MOMENT2_CONC( IMODE )
            ! DRYM3  = WETM3 - DBLE( H2OFAC * AEROSPC_CONC( AH2O_IDX, IMODE ) )
            ! DRYM20 = WETM2 * ( DRYM3 / WETM3 ) ** D_TWOTHIRDS
            ! adj code:
            TEMP1 = DRYM3 / WETM3
!            IF ( TEMP1 .LE. 0.0 ) THEN
!               TEMP1_ADJ = 0.0
!            ELSE
               TEMP1_ADJ = WETM2 * D_TWOTHIRDS * TEMP1 ** ( D_TWOTHIRDS - 1 )
     &                   * DRYM20_ADJ / WETM3
!            END IF
            WETM2_ADJ = TEMP1 ** D_TWOTHIRDS * DRYM20_ADJ
            DRYM3_ADJ = TEMP1_ADJ
            WETM3_ADJ = DRYM3_ADJ - TEMP1 * TEMP1_ADJ
            WETM3_ADJ = WETM3_ADJ + DRYM3_ADJ
            aerospc_concb( AH2O_IDX, IMODE ) = aerospc_concb( AH2O_IDX, 
     &             IMODE ) - H2OFAC * DRYM3_ADJ
            DRYM3_ADJ = 0.0
            moment2_concb( IMODE ) = moment2_concb( IMODE ) + WETM2_ADJ
            WETM2_ADJ = 0.0
            moment3_concb( IMODE ) = moment3_concb( IMODE ) + WETM3_ADJ
            WETM3_ADJ = 0.0

            !------
            ! fwd code:
            ! M3SOA = 0.0d0
            ! adj code:
            M3SOA_ADJ = 0.0

            !------
            ! fwd code:
            ! M3OTHR = SOILFAC * AEROSPC_CONC( ASOIL_IDX,IMODE )
            !      + ANTHFAC * AEROSPC_CONC( ACORS_IDX,IMODE )
            ! adj code:
            aerospc_concb( ASOIL_IDX, IMODE ) = aerospc_concb( ASOIL_IDX,
     &                     IMODE ) + SOILFAC * M3OTHR_ADJ
            aerospc_concb( ACORS_IDX, IMODE ) = aerospc_concb( ACORS_IDX,
     &                     IMODE ) + ANTHFAC * M3OTHR_ADJ
            M3OTHR_ADJ = 0.0

            !------
            ! fwd code:
            ! CINORG( KNH4,IMODE ) = AEROSPC_CONC( ANH4_IDX,IMODE )
            ! CINORG( KNO3,IMODE ) = AEROSPC_CONC( ANO3_IDX,IMODE )
            ! CINORG( KCL, IMODE ) = AEROSPC_CONC( ACL_IDX,IMODE )
            ! CINORG( KSO4,IMODE ) = AEROSPC_CONC( ASO4_IDX,IMODE )
            ! CINORG( KNA, IMODE ) = AEROSPC_CONC( ANA_IDX,IMODE )
            ! CINORG( KHP, IMODE ) = HPLUS( IMODE )
            ! adj code:
            HPLUS_ADJ ( IMODE ) = HPLUS_ADJ ( IMODE ) + CINORG_ADJ ( 
     &                            KHP, IMODE )
            CINORG_ADJ ( KHP, IMODE ) = 0.0
            aerospc_concb( ANA_IDX, IMODE ) = aerospc_concb( ANA_IDX,  
     &                     IMODE ) + CINORG_ADJ ( KNA, IMODE )
            CINORG_ADJ ( KNA, IMODE ) = 0.0
            aerospc_concb( ASO4_IDX, IMODE ) = aerospc_concb( ASO4_IDX,
     &                     IMODE ) + CINORG_ADJ ( KSO4, IMODE )
            CINORG_ADJ ( KSO4, IMODE ) = 0.0
            aerospc_concb( ACL_IDX, IMODE ) = aerospc_concb( ACL_IDX, 
     &                     IMODE ) + CINORG_ADJ ( KCL, IMODE )
            CINORG_ADJ ( KCL, IMODE ) = 0.0
            aerospc_concb( ANO3_IDX, IMODE ) = aerospc_concb( ANO3_IDX,
     &                     IMODE ) + CINORG_ADJ( KNO3, IMODE )
            CINORG_ADJ ( KNO3, IMODE ) = 0.0
            aerospc_concb( ANH4_IDX, IMODE ) = aerospc_concb( ANH4_IDX,
     &                     IMODE ) + CINORG_ADJ( KNH4, IMODE )
            CINORG_ADJ ( KNH4, IMODE ) = 0.0

            IF ( ISTEP > 1 ) THEN
 
               !------
               ! fwd code:
               ! GRFAC1( 3 ) = FCONC_SO4( 3,1 ) 
               ! GRFAC2( 3 ) = FCONC_SO4( 3,2 ) 
               ! adj code:
               FCONC_SO4_ADJ ( 3, 2 ) = FCONC_SO4_ADJ ( 3, 2 ) +
     &                                  GRFAC2_ADJ ( 3 ) 
               GRFAC2_ADJ ( 3 ) = 0.0
               FCONC_SO4_ADJ ( 3, 1 ) = FCONC_SO4_ADJ ( 3, 1 ) +
     &                                  GRFAC1_ADJ ( 1 )
               GRFAC1_ADJ ( 1 ) = 0.0

               ! Checkpointing
               AM0( 3 ) = AM0_CHECK( B )
               AM1( 3 ) = AM1_CHECK( B )
               AM2( 3 ) = AM2_CHECK( B )
               FCONC_SO4( 3, : ) = FCONC_SO4_CHECK( B, : )

               CALL HCOND3_ADJ ( AM0( 3 ), AM0_ADJ( 3 ), AM1( 3 ),
     &                        AM1_ADJ( 3 ), AM2( 3 ), AM2_ADJ( 3 ),
     &                        DV_SO4, ALPHSULF, CBAR_SO4, FCONC_SO4(3,:),
     &                        FCONC_SO4_ADJ ( 3, : ) )

               !------
               ! fwd code:
               ! AM0( 3 ) = MOMENT0_CONC( 3 )
               ! AM1( 3 ) = MOMENT0_CONC( 3 ) * aeromode_DIAM( 3 )
               ! &        * EXP( 0.5 * aeromode_SDEV( 3 ) * aeromode_SDEV( 3 ) )
               ! AM2( 3 ) = MOMENT2_CONC( 3 )
               ! adj code:
               moment2_concb( 3 ) = moment2_concb( 3 ) + AM2_ADJ ( 3 )
               AM2_ADJ ( 3 ) = 0.0
               moment0_concb( 3 ) = moment0_concb( 3 ) + aeromode_DIAM( 3 )
     &                * EXP( 0.5 * aeromode_SDEV( 3 ) * aeromode_SDEV( 3 ) )
     &                * AM1_ADJ ( 3 )
               aeromode_diamb    ( 3 ) = aeromode_diamb    ( 3 ) + 
     &                 MOMENT0_CONC( 3 ) * EXP( 0.5 * aeromode_SDEV( 3 ) * 
     &                 aeromode_SDEV( 3 ) ) * AM1_ADJ ( 3 )
               aeromode_sdevb    ( 3 ) = aeromode_sdevb    ( 3 ) +
     &                 aeromode_DIAM( 3 ) * MOMENT0_CONC( 3 ) *
     &                 aeromode_SDEV ( 3 ) * EXP( 0.5 * aeromode_SDEV( 3 ) 
     &                 * aeromode_SDEV( 3 ) ) * AM1_ADJ ( 3 )
               AM1_ADJ ( 3 ) = 0.0
               moment0_concb( 3 ) = moment0_concb( 3 ) + AM0_ADJ ( 3 )
               AM0_ADJ ( 3 ) = 0.0

            END IF
     
            ISTEP = ISTEP - 1

         END DO

      END IF  ! if ( HYBRID )

      IF ( HYBRID ) THEN
         !------
         ! fwd code:
         ! CONDSO4( 2 ) = OMEGA_AC_SO4 * SCONDRATE
         ! CONDSO4( 3 ) = SCONDRATE - ( CONDSO4( 1 ) + CONDSO4( 2 ) )
         ! adj code:
         SCONDRATE_ADJ = SCONDRATE_ADJ + CONDSO4_ADJ ( 3 )
         CONDSO4_ADJ ( 1 ) = CONDSO4_ADJ ( 1 ) - CONDSO4_ADJ ( 3 )
         CONDSO4_ADJ ( 2 ) = CONDSO4_ADJ ( 2 ) - CONDSO4_ADJ ( 3 )
         CONDSO4_ADJ ( 3 ) = 0.0
         OMEGA_AC_SO4_ADJ = OMEGA_AC_SO4_ADJ + SCONDRATE * CONDSO4_ADJ ( 2 )
         SCONDRATE_ADJ = SCONDRATE_ADJ + OMEGA_AC_SO4 * CONDSO4_ADJ ( 2 )
         CONDSO4_ADJ ( 2 ) = 0.0
      ELSE                                  ! fine equilibrium
         !------
         ! fwd code:
         ! CONDSO4( 2 ) = SCONDRATE - CONDSO4( 1 )
         ! CONDSO4( 3 ) = 0.0D0               ! no coarse mode chemistry
         ! adj code:
         CONDSO4_ADJ ( 3 ) = 0.0
         SCONDRATE_ADJ = SCONDRATE_ADJ + CONDSO4_ADJ ( 2 )
         CONDSO4_ADJ ( 1 ) = CONDSO4_ADJ ( 1 ) - CONDSO4_ADJ ( 2 )
         CONDSO4_ADJ ( 2 ) = 0.0
      END IF

      !------
      ! fwd code:
      ! CONDSO4( 1 ) = OMEGA_AT_SO4 * SCONDRATE
      ! adj code:
      OMEGA_AT_SO4_ADJ = OMEGA_AT_SO4_ADJ + SCONDRATE * CONDSO4_ADJ ( 1 )
      SCONDRATE_ADJ = SCONDRATE_ADJ + OMEGA_AT_SO4 * CONDSO4_ADJ ( 1 )
      CONDSO4_ADJ ( 1 ) = 0.0

      !------
      ! fwd code;
      ! CGR( 1 ) = SO4FAC * SCONDRATE * OMEGA_AT_SO4
      ! CGR( 2 ) = SO4FAC * SCONDRATE * OMEGA_AC_SO4
      ! adj code:
      SCONDRATE_ADJ = SCONDRATE_ADJ + SO4FAC * OMEGA_AC_SO4 * CGR_ADJ ( 2 )
      OMEGA_AC_SO4_ADJ = OMEGA_AC_SO4_ADJ + SO4FAC * SCONDRATE *
     &                   CGR_ADJ ( 2 )
      CGR_ADJ ( 2 ) = 0.0
      SCONDRATE_ADJ = SCONDRATE_ADJ + SO4FAC * OMEGA_AT_SO4 * CGR_ADJ ( 1 )
      OMEGA_AT_SO4_ADJ = OMEGA_AT_SO4_ADJ + SO4FAC * SCONDRATE *
     &                   CGR_ADJ ( 1 ) 
      CGR_ADJ ( 1 ) = 0.0

      !------
      ! fwd code:
      ! FCONCM1_SO4  = 1.0D0 / TMP
      ! OMEGA_AT_SO4 = FCONCM1_SO4 * FCONC_SO4( 1,2 )
      ! OMEGA_AC_SO4 = FCONCM1_SO4 * FCONC_SO4( 2,2 )
      ! adj code:
      FCONCM1_SO4_ADJ = FCONCM1_SO4_ADJ + FCONC_SO4( 2, 2 ) *
     &                  OMEGA_AC_SO4_ADJ
      FCONC_SO4_ADJ ( 2, 2 ) = FCONC_SO4_ADJ ( 2, 2 ) + FCONCM1_SO4 *
     &                         OMEGA_AC_SO4_ADJ
      OMEGA_AC_SO4_ADJ = 0.0
      FCONCM1_SO4_ADJ = FCONCM1_SO4_ADJ + FCONC_SO4( 1, 2 ) *
     &                  OMEGA_AT_SO4_ADJ
      FCONC_SO4_ADJ ( 1, 2 ) = FCONC_SO4_ADJ( 1, 2 ) + FCONCM1_SO4 *
     &                         OMEGA_AT_SO4_ADJ
      OMEGA_AT_SO4_ADJ = 0.0
      TMP_ADJ = TMP_ADJ - ( 1.0 / ( TMP ** 2.0 ) ) * FCONCM1_SO4_ADJ
      FCONCM1_SO4_ADJ = 0.0

      DO I = N_MODE, 1, -1
 
         !------
         ! fwd code:
         ! TMP = TMP + FCONC_SO4( I,2 )
         ! adj code:
         FCONC_SO4_ADJ ( I, 2 ) = FCONC_SO4_ADJ ( I, 2 ) + TMP_ADJ
         
      END DO

      !------
      ! fwd code:
      ! TMP = 0.0
      ! adj code:
      TMP_ADJ = 0.0

      IF ( SO4RATE /= 0 ) THEN

         !------
         ! fwd code:
         ! SCONDRATE = MAX( SO4RATE - DMDT_SO4, 0.0D0 )
         ! adj code:
         IF ( DBLE(SO4RATE) - DMDT_SO4 > 0 ) THEN
            SO4RATE_ADJ = SO4RATE_ADJ + SCONDRATE_ADJ
            DMDT_SO4_ADJ = DMDT_SO4_ADJ - SCONDRATE_ADJ
            SCONDRATE_ADJ = 0.0
         ELSE
            SCONDRATE_ADJ = 0.0
         END IF

         DNDT = DNDT_CHECK
         DMDT_SO4 = DMDT_SO4_CHECK
         DM2DT = DM2DT_CHECK

         !------
         ! fwd code:
         ! CALL NEWPART3 ( AIRRH, AIRTEMP, XH2SO4, SO4RATE,
         !               DNDT, DMDT_SO4, DM2DT )
         ! adj code:
         CALL NEWPART3_ADJ ( AIRRH, AIRTEMP, XH2SO4, XH2SO4_ADJ, SO4RATE,
     &                       SO4RATE_ADJ, DNDT, DNDT_ADJ, DMDT_SO4,
     &                       DMDT_SO4_ADJ, DM2DT, DM2DT_ADJ )

         !------
         ! fwd code:
         ! XH2SO4 = SO4RATE / TMP
         ! XH2SO4 = MAX( XH2SO4, CONMIN )
         ! PRECURSOR_CONC( SULF_IDX ) = XH2SO4
         ! adj code:
         XH2SO4_ADJ = XH2SO4_ADJ + precursor_concb( SULF_IDX ) 
         precursor_concb( SULF_IDX ) = 0.0
         IF ( XH2SO4 < CONMIN ) THEN
            XH2SO4_ADJ = 0.0
         END IF
         SO4RATE_ADJ = SO4RATE_ADJ + ( 1.0 / TMP ) * XH2SO4_ADJ
         TMP_ADJ = TMP_ADJ - ( SO4RATE / ( TMP ** 2.0 ) ) * XH2SO4_ADJ
         XH2SO4_ADJ = 0.0
      
         DO I = N_MODE, 1, -1 
       
            !------
            ! fwd code:
            ! TMP = TMP + FCONC_SO4( I,2 )
            ! adj code:
            FCONC_SO4_ADJ ( I, 2 ) = FCONC_SO4_ADJ ( I, 2 ) + TMP_ADJ

         END DO

         !------
         ! fwd code:
         ! TMP = 0.0
         ! adj code:
         TMP_ADJ = 0.0

      END IF  ! SO4RATE /= 0

      !------
      ! fwd code:
      ! DMDT_SO4  = 0.0D0
      ! DNDT      = 0.0D0
      ! DM2DT     = 0.0D0
      ! SCONDRATE = 0.0D0
      ! adj code:
      SCONDRATE_ADJ = 0.0
      DM2DT_ADJ = 0.0
      DNDT_ADJ = 0.0
      DMDT_SO4_ADJ = 0.0

      DO I = N_MODE, 1, -1 
   
         !------
         ! fwd code:
         ! GRFAC1( I ) = FCONC_SO4( I,1 )
         ! GRFAC2( I ) = FCONC_SO4( I,2 )
         ! adj code:
          FCONC_SO4_ADJ ( I, 2 ) = FCONC_SO4_ADJ ( I, 2 ) + 
     &             GRFAC2_ADJ ( I )
          GRFAC2_ADJ ( I ) = 0.0
          FCONC_SO4_ADJ ( I, 1 ) = FCONC_SO4_ADJ ( I, 1 ) +
     &             GRFAC1_ADJ ( I )
          GRFAC1_ADJ ( I ) = 0.0

      END DO

      IF ( .NOT. HYBRID ) THEN
         !------
         ! fwd code:
         ! FCONC_SO4( N_MODE,1 ) = 0.D0
         ! FCONC_SO4( N_MODE,2 ) = 0.D0
         ! adj code:
         FCONC_SO4_ADJ ( N_MODE, 2 ) = 0.0
         FCONC_SO4_ADJ ( N_MODE, 1 ) = 0.0
      END IF

      DO I = N_MODE, 1, -1

         ! Checkpointing
         AM0( I ) = AM0_CHECK1( I )
         AM1( I ) = AM1_CHECK1( I )
         AM2( I ) = AM2_CHECK1( I )
         FCONC_SO4( I, : ) = FCONC_SO4_CHECK1( I, : )

         !------               
         ! fwd code:
         ! CALL HCOND3( AM0( I ), AM1( I ), AM2( I ), DV_SO4, ALPHSULF, 
         !               CBAR_SO4, FCONC_SO4( I,: ) )
         ! adj code:
         CALL HCOND3_ADJ ( AM0( I ) , AM0_ADJ( I ), AM1( I ),
     &               AM1_ADJ( I ), AM2( I ), AM2_ADJ( I ),
     &               DV_SO4, ALPHSULF, CBAR_SO4, FCONC_SO4( I, : ),
     &               FCONC_SO4_ADJ ( I, : ) )

      END DO

      DO I = N_MODE, 1, -1
 
         !------
         ! fwd code:
         ! AM0( I ) = MOMENT0_CONC( I ) 
         ! AM1( I ) = MOMENT0_CONC( I ) * aeromode_DIAM( I )
         !        * EXP(0.5 * aeromode_SDEV( I ) * aeromode_SDEV( I ) )
         ! AM2( I ) = MOMENT2_CONC( I )
         ! adj code:
         moment2_concb( I ) = moment2_concb( I ) + AM2_ADJ ( I )
         AM2_ADJ ( I ) = 0.0
         moment0_concb( I ) = moment0_concb( I ) + aeromode_DIAM( I )
     &          * EXP(0.5 * aeromode_SDEV( I ) * aeromode_SDEV( I ) ) *
     &          AM1_ADJ ( I )
         aeromode_diamb    ( I ) = aeromode_diamb    ( I ) +
     &           MOMENT0_CONC( I ) * EXP(0.5 * aeromode_SDEV( I ) * 
     &           aeromode_SDEV( I ) ) * AM1_ADJ ( I )
         aeromode_sdevb    ( I ) = aeromode_sdevb    ( I ) +
     &           MOMENT0_CONC( I ) * aeromode_DIAM ( I ) * aeromode_SDEV( I )
     &           * EXP(0.5 * aeromode_SDEV( I ) * aeromode_SDEV( I ) ) *
     &           AM1_ADJ ( I )
         AM1_ADJ ( I ) = 0.0
         moment0_concb( I ) = moment0_concb ( I ) + AM0_ADJ ( I )
         AM0_ADJ ( I ) = 0.0
   
      END DO

      DO I = N_MODE, 1, -1

         DO N = N_AEROSPC, 1, -1

            !------
            ! fwd code:
            ! HPLUS( I ) = HPLUS( I )
            !          + AEROSPC( N )%CHARGE * AEROSPC_CONC( N,I ) / AEROSPC_MW( N )
            ! adj code:
            aerospc_concb( N, I ) = aerospc_concb( N, I ) +
     &                   AEROSPC( N )%CHARGE / AEROSPC_MW( N ) * 
     &                   HPLUS_ADJ ( I )

         END DO 

      END DO

      !------
      ! fwd code:
      ! HPLUS = 0.0
      ! adj code:
      HPLUS_ADJ = 0.0

      DO I = N_AEROSPC, 1, -1

         IF ( AEROSPC( I )%CHARGE /= 0 ) THEN
            !------
            ! fwd code:
            ! TMP = TMP + AEROSPC_CONC( I,N_MODE )
            ! adj code:
            aerospc_concb( I, N_MODE ) = aerospc_concb( I, N_MODE ) +
     &             TMP_ADJ
         END IF

      END DO

      !------
      ! fwd code:
      ! TMP = 0.0
      ! adj code:
      TMP_ADJ = 0.0

      END SUBROUTINE

