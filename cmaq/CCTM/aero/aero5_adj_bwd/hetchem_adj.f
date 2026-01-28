
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
      SUBROUTINE HETCHEM_ADJ   ( GAMMA, DT )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of subroutine that calculates the heterogeneous conversion 
!   of N2O5 to HNO3.
!
! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      USE AERO_DATA_B
      USE PRECURSOR_DATA_B
      USE MET_DATA

      IMPLICIT NONE

! *** ARGUMENTS
!
! N2O5->NO3 rxn probability
      REAL, INTENT(OUT) :: gamma
! Synchronization time step
      REAL, INTENT(IN) :: dt
!
!
! *** PARAMETERS
!
! molecular diffusivity
! of N2O5 at 101325 Pa
! and 273.15 K [m2/sec]
      REAL, PARAMETER :: std_diff_n2o5=0.1e-4
! g/kg unit conversion
      REAL, PARAMETER :: gpkg=1.0e+03
!
!
! *** LOCAL VARIABLES
!
! *** chemical species concentrations
!
! gas-phase nitric acid [ug/m3]
      REAL :: ghno3
! gas-phase dinitrogen pentoxide [ug/m3]
      REAL :: gn2o5
! gas-phase NO2 [ug/m3]
      REAL :: gno2
! gas-phase HONO [ug/m3]
      REAL :: ghono
!
! *** 2nd and 3rd moments before equilibration (without H2O)
!
      REAL :: old_m3_i, old_m3_j
      REAL :: old_m2_i, old_m2_j
!
! *** variables for N2O5 + H2O -> 2 HNO3 conversion
!
! M3 before equilibration w.H2O
      REAL :: wet_m3_i, wet_m3_j
! M2 before equilibration w.H2O
      REAL :: wet_m2_i, wet_m2_j
! Initial median diameter w.H2O
      REAL :: dg_at_wet, dg_ac_wet
! Initial effective diameter w.H2O
      REAL :: de_at_wet, de_ac_wet
! modal factors to calculate KN2O5
      REAL :: xxf_at, xxf_ac
! molecular velocity (m/s)
      REAL :: cbar
! ambient molecular diffusivity [m2/s]
      REAL :: diff_n2o5
! function to compute GAMMA
      REAL :: N2O5PROB
! pseudo-first order rate constant
      REAL :: kn2o5
! fraction of N2O5 left after chemical rxn
      REAL :: expdt_n2o5
!
! *** saved variables for conversion
! converts ug -> mol
      REAL, SAVE :: faerno2
! converts ug -> mol
      REAL, SAVE :: faern2o5
! converts mol -> ug
      REAL, SAVE :: dfhno3
! converts mol -> ug
      REAL, SAVE :: dfhono
!
! *** variables for 2 NO2 + H2O -> HONO + HNO3 conversion
!
! pseudo-first order rate constant
      REAL :: kno2
! fraction of NO2 left after chemical rxn
      REAL :: expdt_no2
! aerosol surface area (m**2/m**3)
      REAL :: totsurfa
!
! *** first time switch
      LOGICAL, SAVE :: firstime=.true.
      INTRINSIC EXP
      INTRINSIC MAX
      INTRINSIC LOG
      INTRINSIC SQRT

! Adjoint variables
      
      REAL       :: TOTSURFA_ADJ
      REAL       :: GHONO_ADJ
      REAL       :: KNO2_ADJ
      REAL       :: WET_M2_I_ADJ, WET_M2_J_ADJ
      REAL       :: WET_M3_I_ADJ, WET_M3_J_ADJ
      REAL       :: XXF_AT_ADJ
      REAL       :: XXF_AC_ADJ
      REAL       :: DG_AC_WET_ADJ, DG_AT_WET_ADJ
      REAL       :: KN2O5_ADJ
      REAL       :: EXPDT_N2O5_ADJ
      REAL       :: EXPDT_NO2_ADJ
      REAL       :: DE_AC_WET_ADJ
      REAL       :: DE_AT_WET_ADJ
      REAL       :: OLD_M2_I_ADJ
      REAL       :: OLD_M2_J_ADJ
      REAL       :: OLD_M3_I_ADJ
      REAL       :: OLD_M3_J_ADJ
      REAL       :: GNO2_ADJ
      REAL       :: GN2O5_ADJ
      REAL       :: GHNO3_ADJ
      REAL       :: GAMMA_ADJ

! Checkpointing variables
      REAL       :: GN2O5_CHECK, GNO2_CHECK, aeromode_sdev_check( 3 )
      REAL       :: aeromode_diam_check( 3 ), KN2O5_CHK

!-------------------------- Begin Execution -----------------

      TOTSURFA_ADJ = 0.0
      GHONO_ADJ = 0.0
      KNO2_ADJ = 0.0
      WET_M2_I_ADJ = 0.0
      WET_M2_J_ADJ = 0.0
      WET_M3_I_ADJ = 0.0
      WET_M3_J_ADJ = 0.0
      XXF_AT_ADJ = 0.0
      XXF_AC_ADJ = 0.0
      DG_AC_WET_ADJ = 0.0
      DG_AT_WET_ADJ = 0.0
      KN2O5_ADJ = 0.0 
      EXPDT_N2O5_ADJ = 0.0
      EXPDT_NO2_ADJ = 0.0
      DE_AC_WET_ADJ = 0.0
      DE_AT_WET_ADJ = 0.0
      OLD_M2_I_ADJ = 0.0 
      OLD_M2_J_ADJ = 0.0
      OLD_M3_I_ADJ = 0.0 
      OLD_M3_J_ADJ = 0.0
      GNO2_ADJ = 0.0
      GN2O5_ADJ = 0.0
      GHNO3_ADJ = 0.0
      GAMMA_ADJ = 0.0

! Forward Code

C *** compute only on first pass

      If ( firstime ) Then
        firstime = .false.
        FAERNO2 = 1.0E-6 / precursor_mw( NO2_IDX )
        FAERN2O5 = 1.0E-6 / precursor_mw( N2O5_IDX )
        DFHONO = precursor_mw( HONO_IDX ) / 1.0E-6
        DFHNO3 = precursor_mw( HNO3_IDX ) / 1.0E-6
      Endif   ! first time condition

c *** fetch vapor-phase concentrations [ug/m3]

      GHNO3 = precursor_conc( HNO3_IDX )
      GN2O5 = precursor_conc( N2O5_IDX )
      GNO2  = precursor_conc( NO2_IDX )
      GHONO = precursor_conc( HONO_IDX )

c *** set up variables needed for calculating KN2O5

c *** capture values of "dry" 2nd and 3rd moments before equilibration
c     the folowing code assumes that GETPAR has been called with
c     M3_WET_FLAG set to .FALSE. and that the 2nd and 3rd moments have
c     been adjusted for the new SOA.

      OLD_M3_I = moment3_conc( 1 )
      OLD_M3_J = moment3_conc( 2 )
      OLD_M2_I = moment2_conc( 1 )
      OLD_M2_J = moment2_conc( 2 )

c *** compute GAMMA as function of TEMP, RH, & particle composition
c     Note: the last argument to this function can be changed to use 
c     a different parameterization of GAMMA.

      GAMMA = N2O5PROB( AIRTEMP, AIRRH, 0 )

c *** calculate molecular speed (m/s) using Eq 4 of Pleim et al (1995)

      CBAR = SQRT( 8.0 * RGASUNIV * AIRTEMP * GPKG
     &     / ( PI * precursor_mw(N2O5_IDX) ))

c *** correct molecular diffusivity for ambient conditions

      DIFF_N2O5 = STD_DIFF_N2O5
     &          * ( ( AIRTEMP / STDTEMP ) ** 1.75 )
     &          * ( STDATMPA / AIRPRS )

c *** estimate the "wet third moments" by adding aerosol water
c      Note: this is the H2O concentration from previous time step

      WET_M3_I = OLD_M3_I + H2OFAC * aerospc_conc( AH2O_IDX,1 )
      WET_M3_J = OLD_M3_J + H2OFAC * aerospc_conc( AH2O_IDX,2 )

c *** calculate "wet second moment" assuming that H2O does not
c     affect the geometric standard deviation

      WET_M2_I = OLD_M2_I * ( WET_M3_I / OLD_M3_I ) ** ( 2.0/3.0 )
      WET_M2_J = OLD_M2_J * ( WET_M3_J / OLD_M3_J ) ** ( 2.0/3.0 )

c *** calculate "wet" geometric mean (same as median) diameters

      DG_AT_WET = aeromode_diam( 1 ) * SQRT( WET_M2_I / OLD_M2_I )
      DG_AC_WET = aeromode_diam( 2 ) * SQRT( WET_M2_J / OLD_M2_J )

C *** calculate effective diameters using Eq 3 of Pleim et al (1995)

      DE_AT_WET = DG_AT_WET * EXP( 1.5
     &          * ( LOG( EXP( aeromode_sdev( 1 ) ) ) ** 2.0 ) )
      DE_AC_WET = DG_AC_WET * EXP( 1.5
     &          * ( LOG( EXP( aeromode_sdev( 2 ) ) ) ** 2.0 ) )

      ! Checkpointing
      aeromode_sdev_check = aeromode_sdev
      aeromode_diam_check = aeromode_diam

c *** calculate pseudo-first order rate constant using Eq 2 of
c     Pleim et al (1995)

      XXF_AT = WET_M2_I /
     &         ( 4.0 + 0.5 * DE_AT_WET * GAMMA * CBAR / DIFF_N2O5 )
      XXF_AC = WET_M2_J /
     &         ( 4.0 + 0.5 * DE_AC_WET * GAMMA * CBAR / DIFF_N2O5 )
      KN2O5 =   GAMMA * CBAR * PI * ( XXF_AT + XXF_AC )

c *** calculate fraction of N2O5 remaining after chemical reaction

      ! Checkpointing
      KN2O5_CHK = KN2O5

      EXPDT_N2O5 = EXP( - KN2O5 * DT )

c *** set up variables needed for calculating KNO2

c *** calculate aerosol surface area

      TOTSURFA  = ( WET_M2_I + WET_M2_J ) * PI

c *** calculate pseudo-first order rate constant using Eq 1 of Vogel
c     et al. (2003). Units of KNO2 is in 1/min in the paper; divide it
c     by 60 to convert it into 1/sec

      KNO2 = MAX ( 0.0, 5.0E-5 * TOTSURFA )

c *** calculate fraction of NO2 remaining after chemical reaction

      EXPDT_NO2 = EXP( -2.0 * KNO2 * DT )

c *** compute new gas-phase concs after heterogeneous reactions occur

c *** adjust nitrous acid for contribution from NO2

      GHONO = GHONO
     &      + ( 0.5 * GNO2  * FAERNO2  * DFHONO ) * ( 1.0 - EXPDT_NO2 )

c *** adjust nitric acid for contributions from N2O5 and NO2

      GHNO3 = GHNO3
     &      + ( 2.0 * GN2O5 * FAERN2O5 * DFHNO3 ) * ( 1.0 - EXPDT_N2O5 )
     &      + ( 0.5 * GNO2  * FAERNO2  * DFHNO3 ) * ( 1.0 - EXPDT_NO2 )

      ! Checkpointing
      GN2O5_CHECK = GN2O5
      GNO2_CHECK = GNO2

c *** adjust N2O5 for heterogeneous loss

      GN2O5 = GN2O5 * EXPDT_N2O5

c *** adjust NO2 for heterogeneous loss

      GNO2  = GNO2  * EXPDT_NO2

! Adjoint Code

      !------
      ! fwd code:
      ! precursor_conc( HNO3_IDX ) = MAX( GHNO3, CONMIN )
      ! precursor_conc( N2O5_IDX ) = MAX( GN2O5, CONMIN )
      ! precursor_conc( NO2_IDX )  = MAX( GNO2, CONMIN )
      ! precursor_conc( HONO_IDX ) = MAX( GHONO, CONMIN )
      ! adj code:
      IF ( GHONO > CONMIN ) THEN
         GHONO_ADJ = precursor_concb( HONO_IDX )
         precursor_concb( HONO_IDX ) = 0.0
      ELSE
         precursor_concb( HONO_IDX ) = 0.0
      END IF
      IF ( GNO2 > CONMIN ) THEN
         GNO2_ADJ = precursor_concb( NO2_IDX ) 
         precursor_concb( NO2_IDX ) = 0.0
      ELSE
         precursor_concb( NO2_IDX ) = 0.0
      END IF
      IF ( GN2O5 > CONMIN ) THEN
         GN2O5_ADJ = precursor_concb( N2O5_IDX ) 
         precursor_concb( N2O5_IDX ) = 0.0
      ELSE
         precursor_concb( N2O5_IDX ) = 0.0
      END IF
      IF ( GHNO3 > CONMIN ) THEN
         GHNO3_ADJ = precursor_concb( HNO3_IDX ) 
         precursor_concb( HNO3_IDX ) = 0.0
      ELSE
         precursor_concb( HNO3_IDX ) = 0.0
      END IF

      ! Checkpointing
      GN2O5 = GN2O5_CHECK
      GNO2 = GNO2_CHECK

      !------
      ! fwd code:
      ! GNO2  = GNO2  * EXPDT_NO2
      ! adj code:
      EXPDT_NO2_ADJ = GNO2 * GNO2_ADJ
      GNO2_ADJ = EXPDT_NO2 * GNO2_ADJ

      !------
      ! fwd code:
      ! GN2O5 = GN2O5 * EXPDT_N2O5
      ! adj code:
      EXPDT_N2O5_ADJ = GN2O5 * GN2O5_ADJ
      GN2O5_ADJ = EXPDT_N2O5 * GN2O5_ADJ
 
      !------
      ! fwd code:
      ! GHNO3 = GHNO3
      !     + ( 2.0 * GN2O5 * FAERN2O5 * DFHNO3 ) * ( 1.0 - EXPDT_N2O5 )
      !     + ( 0.5 * GNO2  * FAERNO2  * DFHNO3 ) * ( 1.0 - EXPDT_NO2 )
      ! adj code:
      GN2O5_ADJ = GN2O5_ADJ + 2.0 * FAERN2O5 * DFHNO3 * 
     &            ( 1.0 - EXPDT_N2O5 ) * GHNO3_ADJ
      EXPDT_N2O5_ADJ = EXPDT_N2O5_ADJ - 2.0 * FAERN2O5 * DFHNO3 * GN2O5 
     &                 * GHNO3_ADJ
      GNO2_ADJ = GNO2_ADJ + 0.5 * FAERNO2 * DFHNO3 * 
     &           ( 1.0 - EXPDT_NO2 ) * GHNO3_ADJ
      EXPDT_NO2_ADJ = EXPDT_NO2_ADJ - 0.5 * FAERNO2 * DFHNO3 * GNO2 
     &                * GHNO3_ADJ

      !-----
      ! fwd code:
      ! GHONO = GHONO
      ! &  + ( 0.5 * GNO2  * FAERNO2  * DFHONO ) * ( 1.0 - EXPDT_NO2 )
      ! adj code:
      GNO2_ADJ = GNO2_ADJ + 0.5 * FAERNO2 * DFHONO * ( 1.0 - EXPDT_NO2 )
     &           * GHONO_ADJ
      EXPDT_NO2_ADJ = EXPDT_NO2_ADJ - 0.5 * FAERNO2 * DFHONO 
     &                * GNO2 * GHONO_ADJ

      !-----
      ! fwd code:
      ! EXPDT_NO2 = EXP( -2.0 * KNO2 * DT )
      ! adj code:
      KNO2_ADJ = - 2.0 * DT * EXP( -2.0 * KNO2 * DT ) 
     &            * EXPDT_NO2_ADJ
      EXPDT_NO2_ADJ = 0.0

      !-----
      ! fwd code:
      ! KNO2 = MAX ( 0.0, 5.0E-5 * TOTSURFA )
      ! adj code:
      IF ( 5.0E-5 * TOTSURFA > 0.0 ) THEN
         TOTSURFA_ADJ = 5.0E-5 * KNO2_ADJ
         KNO2_ADJ = 0.0
      ELSE
         KNO2_ADJ = 0.0
      END IF

      !-----
      ! fwd code:
      ! TOTSURFA = ( WET_M2_I + WET_M2_J ) * PI
      ! adj code:
      WET_M2_I_ADJ = WET_M2_I_ADJ + PI * TOTSURFA_ADJ
      WET_M2_J_ADJ = WET_M2_J_ADJ + PI * TOTSURFA_ADJ
      TOTSURFA_ADJ = 0.0

      ! Checkpointing
      KN2O5 = KN2O5_CHK

      !-----
      ! fwd code:
      ! EXPDT_N2O5 = EXP( - KN2O5 * DT )
      ! adj code:
      KN2O5_ADJ = KN2O5_ADJ - DT * EXP( - KN2O5 * DT ) 
     &           * EXPDT_N2O5_ADJ
      EXPDT_N2O5_ADJ = 0.0

      !-----
      ! fwd code:
      ! KN2O5 = GAMMA * CBAR * PI * ( XXF_AT + XXF_AC )
      ! adj code:
      GAMMA_ADJ = CBAR * PI * ( XXF_AT + XXF_AC ) * KN2O5_ADJ
      XXF_AT_ADJ = XXF_AT_ADJ + GAMMA * CBAR * PI * KN2O5_ADJ
      XXF_AC_ADJ = XXF_AC_ADJ + GAMMA * CBAR * PI * KN2O5_ADJ
      KN2O5_ADJ = 0.0

      !-----
      ! fwd code:
      ! XXF_AC = WET_M2_J /
      ! &     ( 4.0 + 0.5 * DE_AC_WET * GAMMA * CBAR / DIFF_N2O5 )
      ! adj code:
      WET_M2_J_ADJ = WET_M2_J_ADJ + ( 1.0 / 
     &          ( 4.0 + 0.5 * DE_AC_WET * GAMMA * CBAR / DIFF_N2O5 ) ) 
     &          * XXF_AC_ADJ
      DE_AC_WET_ADJ = DE_AC_WET_ADJ - 0.5 * WET_M2_J * ( GAMMA * CBAR 
     &              / DIFF_N2O5 ) / ( ( 4.0 + 0.5 * DE_AC_WET * GAMMA 
     &              * CBAR / DIFF_N2O5 ) ** 2.0 ) * XXF_AC_ADJ
      GAMMA_ADJ = GAMMA_ADJ - 0.5 * WET_M2_J * ( DE_AC_WET * CBAR 
     &          / DIFF_N2O5 ) / ( ( 4.0 + 0.5 * DE_AC_WET * GAMMA 
     &          * CBAR / DIFF_N2O5 ) ** 2.0 ) * XXF_AC_ADJ
      XXF_AC_ADJ = 0.0

      !-----
      ! fwd code:
      ! XXF_AT = WET_M2_I /
      ! &  ( 4.0 + 0.5 * DE_AT_WET * GAMMA * CBAR / DIFF_N2O5 )
      ! adj code:
      WET_M2_I_ADJ = WET_M2_I_ADJ + ( 1.0 / 
     &          ( 4.0 + 0.5 * DE_AT_WET * GAMMA * CBAR / DIFF_N2O5 ) ) 
     &          * XXF_AT_ADJ
      DE_AT_WET_ADJ = DE_AT_WET_ADJ - 0.5 * WET_M2_I * ( GAMMA * CBAR 
     &              / DIFF_N2O5 ) / ( ( 4.0 + 0.5 * DE_AT_WET * GAMMA 
     &              * CBAR / DIFF_N2O5 ) ** 2.0 ) * XXF_AT_ADJ
      GAMMA_ADJ = GAMMA_ADJ - 0.5 * WET_M2_I * ( DE_AT_WET * CBAR 
     &          / DIFF_N2O5 ) / ( ( 4.0 + 0.5 * DE_AT_WET * GAMMA 
     &          * CBAR / DIFF_N2O5 ) ** 2.0 ) * XXF_AT_ADJ
      XXF_AT_ADJ = 0.0

      ! Checkpointing
      aeromode_diam = aeromode_diam_check
      aeromode_sdev = aeromode_sdev_check

      !-----
      ! fwd code:
      ! DE_AC_WET = DG_AC_WET * EXP( 1.5
      ! &    * ( LOG( EXP( aerosize_sdev( 2 ) ) ) ** 2.0 ) )
      ! adj code:
      DG_AC_WET_ADJ = DG_AC_WET_ADJ + EXP( 1.5 
     &          * ( LOG( EXP( aeromode_sdev( 2 )  ) ) ** 2.0 ) ) 
     &          * DE_AC_WET_ADJ
      aeromode_sdevb   ( 2 )  = aeromode_sdevb( 2 ) +
     &            3.0 * DG_AC_WET * EXP( 1.5   
     &          * ( LOG( EXP( aeromode_sdev( 2 )  ) ) ** 2.0 ) ) 
     &          * LOG( EXP( aeromode_sdev( 2 )  ) ) * DE_AC_WET_ADJ
      DE_AC_WET_ADJ = 0.0

      !-----
      ! fwd code:
      ! DE_AT_WET = DG_AT_WET * EXP( 1.5
      ! &     * ( LOG( EXP( aerosize_sdev( 1 ) ) ) ** 2.0 ) )
      ! adj code:
      DG_AT_WET_ADJ = DG_AT_WET_ADJ + EXP( 1.5 
     &     * ( LOG( EXP( aeromode_sdev( 1 )  ) ) ** 2.0 ) )
     &     * DE_AT_WET_ADJ
      aeromode_sdevb   ( 1 )  = aeromode_sdevb( 1 ) +
     &            3.0 * DG_AT_WET * EXP( 1.5   
     &          * ( LOG( EXP( aeromode_sdev( 1 ) ) ) ** 2.0 ) ) 
     &          * LOG( EXP( aeromode_sdev( 1 ) ) ) * DE_AT_WET_ADJ
      DE_AT_WET_ADJ = 0.0

      !-----
      ! fwd code:
      ! DG_AC_WET = aerosize_diam( 2 ) * SQRT( WET_M2_J / OLD_M2_J )
      ! adj code:
      WET_M2_J_ADJ = WET_M2_J_ADJ + aeromode_diam( 2 ) * ( 0.5 
     &             * ( 1 / OLD_M2_J ) / SQRT( WET_M2_J / OLD_M2_J ) ) 
     &             * DG_AC_WET_ADJ
      OLD_M2_J_ADJ = OLD_M2_J_ADJ - aeromode_diam( 2 ) * ( 0.5 
     &             * ( WET_M2_J / ( OLD_M2_J ** 2.0 ) ) 
     &             / SQRT( WET_M2_J / OLD_M2_J ) ) * DG_AC_WET_ADJ
      aeromode_diamb   ( 2 ) = aeromode_diamb( 2 ) +
     &              SQRT( WET_M2_J / OLD_M2_J ) * DG_AC_WET_ADJ
      DG_AC_WET_ADJ = 0.0

      !-----
      ! fwd code:
      ! DG_AT_WET = aerosize_diam( 1 ) * SQRT( WET_M2_I / OLD_M2_I )
      ! adj code:
      WET_M2_I_ADJ = WET_M2_I_ADJ + aeromode_diam( 1 ) * ( 0.5 
     &             * ( 1 / OLD_M2_I ) / SQRT( WET_M2_I / OLD_M2_I ) ) 
     &             * DG_AT_WET_ADJ
      OLD_M2_I_ADJ = OLD_M2_I_ADJ - aeromode_diam( 1 ) * ( 0.5 
     &             * ( WET_M2_I / ( OLD_M2_I ** 2.0 ) ) 
     &             / SQRT( WET_M2_I / OLD_M2_I ) ) * DG_AT_WET_ADJ
      aeromode_diamb   ( 1 ) = aeromode_diamb( 1 ) +
     &            SQRT( WET_M2_I / OLD_M2_I ) * DG_AT_WET_ADJ
      DG_AT_WET_ADJ = 0.0

      !-----
      ! fwd code:
      ! WET_M2_J = OLD_M2_J * ( WET_M3_J / OLD_M3_J ) ** ( 2.0/3.0 )
      ! adj code:
      OLD_M2_J_ADJ = OLD_M2_J_ADJ + ( WET_M3_J / OLD_M3_J ) ** ( 2.0/3.0 )
     &             * WET_M2_J_ADJ
      WET_M3_J_ADJ = WET_M3_J_ADJ + ( 2.0/3.0 ) * OLD_M2_J * ( 1.0 / 
     &          OLD_M3_J ) * ( ( WET_M3_J / OLD_M3_J ) ** ( -1.0/3.0 ) 
     &          ) * WET_M2_J_ADJ
      OLD_M3_J_ADJ = OLD_M3_J_ADJ - ( 2.0/3.0 ) * OLD_M2_J 
     &             * ( WET_M3_J / ( OLD_M3_J ** 2.0 ) ) 
     &             * ( ( WET_M3_J / OLD_M3_J ) ** ( -1.0/3.0 ) ) 
     &             * WET_M2_J_ADJ
      WET_M2_J_ADJ = 0.0

      !-----
      ! fwd code:
      ! WET_M2_I = OLD_M2_I * ( WET_M3_I / OLD_M3_I ) ** ( 2.0/3.0 )
      ! adj code:
      OLD_M2_I_ADJ = OLD_M2_I_ADJ + ( WET_M3_I / OLD_M3_I ) ** ( 2.0/3.0 )
     &             * WET_M2_I_ADJ
      WET_M3_I_ADJ = WET_M3_I_ADJ + ( 2.0/3.0 ) * OLD_M2_I * ( 1.0 / 
     &          OLD_M3_I ) * ( ( WET_M3_I / OLD_M3_I ) ** ( -1.0/3.0 ) 
     &          ) * WET_M2_I_ADJ
      OLD_M3_I_ADJ = OLD_M3_I_ADJ - ( 2.0/3.0 ) * OLD_M2_I 
     &             * ( WET_M3_I / ( OLD_M3_I ** 2.0 ) ) 
     &             / ( ( WET_M3_I / OLD_M3_I ) ** ( 1.0/3.0 ) ) 
     &             * WET_M2_I_ADJ
      WET_M2_I_ADJ = 0.0

      !-----
      ! fwd code:
      ! WET_M3_J = OLD_M3_J + H2OFAC * aerospc_conc( AH2O_IDX,2 )
      ! adj code:
      OLD_M3_J_ADJ = OLD_M3_J_ADJ + WET_M3_J_ADJ
      aerospc_concb( AH2O_IDX, 2 ) = aerospc_concb( AH2O_IDX, 2 ) +
     &             H2OFAC * WET_M3_J_ADJ
      WET_M3_J_ADJ = 0.0

      !-----
      ! fwd code:
      ! WET_M3_I = OLD_M3_I + H2OFAC * aerospc_conc( AH2O_IDX,1 )
      ! adj code:
      OLD_M3_I_ADJ = OLD_M3_I_ADJ + WET_M3_I_ADJ
      aerospc_concb( AH2O_IDX, 1 ) = aerospc_concb( AH2O_IDX, 1 ) +
     &             H2OFAC * WET_M3_I_ADJ
      WET_M3_I_ADJ = 0.0

      !-----
      ! fwd code:
      ! GAMMA = N2O5PROB( AIRTEMP, AIRRH, 0 )
      ! adj code:
      CALL N2O5PROB_ADJ( AIRTEMP, AIRRH, 0, GAMMA_ADJ )

      !------
      ! fwd code:
      ! OLD_M3_I = moment3_conc( 1 )
      ! OLD_M3_J = moment3_conc( 2 )
      ! OLD_M2_I = moment2_conc( 1 )
      ! OLD_M2_J = moment2_conc( 2 )
      ! adj code:
      moment2_concb( 2 ) = moment2_concb( 2 ) + OLD_M2_J_ADJ
      moment2_concb( 1 ) = moment2_concb( 1 ) + OLD_M2_I_ADJ
      moment3_concb( 2 ) = moment3_concb( 2 ) + OLD_M3_J_ADJ
      moment3_concb( 1 ) = moment3_concb( 1 ) + OLD_M3_I_ADJ
      OLD_M2_I_ADJ= 0.0
      OLD_M2_J_ADJ = 0.0
      OLD_M3_I_ADJ = 0.0
      OLD_M3_J_ADJ = 0.0

      !------
      ! fwd code:
      ! GHNO3 = precursor_conc( HNO3_IDX )
      ! GN2O5 = precursor_conc( N2O5_IDX )
      ! GNO2  = precursor_conc( NO2_IDX )
      ! GHONO = precursor_conc( HONO_IDX )
      ! adj code:
      precursor_concb( HONO_IDX ) = precursor_concb( HONO_IDX ) + GHONO_ADJ
      precursor_concb( NO2_IDX  ) = precursor_concb( NO2_IDX  ) + GNO2_ADJ
      precursor_concb( N2O5_IDX ) = precursor_concb( N2O5_IDX ) + GN2O5_ADJ
      precursor_concb( HNO3_IDX ) = precursor_concb( HNO3_IDX ) + GHNO3_ADJ
      GHNO3_ADJ = 0.0
      GN2O5_ADJ = 0.0
      GNO2_ADJ = 0.0
      GHONO_ADJ = 0.0

      END SUBROUTINE HETCHEM_ADJ
