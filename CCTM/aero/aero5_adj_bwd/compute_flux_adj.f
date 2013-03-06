
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
      SUBROUTINE COMPUTE_FLUX_ADJ ( NVOLINORG, GNH3R8, GNH3R8_ADJ, 
     &                             GNO3R8, GNO3R8_ADJ, GCLR8, GCLR8_ADJ,
     &                             KNH4, KNO3, KCL, CEQ, CEQ_ADJ, 
     &                             CONDRATE, CONDRATE_ADJ, HPLUS, 
     &                             HPLUS_ADJ, RATE, RATE_ADJ, J, J_ADJ )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of subroutine that determines the evaporative/condensational
!   flux of volatile inorganic species to aerosol modes.
!
! INPUTS:
!   nvolinorg, GNH3R8, GCLR8, KNH4, KNO3, KCL, CEQ, CONDRATE, HPLUS,
!   RATE, CEQ_ADJ, J_ADJ
!
! OUTPUTS
!   GNH3R8_ADJ, GNO3R8_ADJ, GCLR8_ADJ, CEQ_ADJ, CONDRATE_ADJ, 
!   HPLUS_ADJ, RATE_ADJ
!
! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      USE AERO_DATA_B
      USE MET_DATA
      USE PRECURSOR_DATA_B

      IMPLICIT NONE

!     Arguments:
      INTEGER      NVOLINORG
      REAL( 8 ) :: GNH3R8, GNO3R8, GCLR8 ! gas concentrations [ug/m3]
      INTEGER      KNH4, KNO3, KCL       ! Indices to species
      REAL( 8 ) :: CEQ( NVOLINORG )      ! vapor concentrations [mol/m3]
      REAL( 8 ) :: CONDRATE              ! effective condensation rate (I) for 3rd moment
      REAL( 8 ) :: HPLUS                 ! hydrogen ion concentration for mode [mol/m3]
      REAL( 8 ) :: RATE
      REAL( 8 ) :: J( NVOLINORG )        ! molar cond./evap. flux [mol/m3-s]

!     Local Variables:
      REAL( 8 ),  PARAMETER :: AFACT = 1.0D-01  ! factor for H+ limiter
      REAL( 8 ),  PARAMETER :: SMALL = 1.0D-25
      REAL( 8 ) :: CINF( NVOLINORG ) ! gas concentration in mol/m3
      REAL( 8 ) :: QK              ! factor for modifying vapor press. based on H+ limit
      REAL( 8 ) :: HFLUX           ! flux of H+ to mode from cond/evap
      REAL( 8 ) :: HLIM            ! maximum allowable H+ flux to mode
      REAL( 8 ) :: AA, BB, CC      ! terms in quadratic equation
      REAL( 8 ) :: JH2SO4          ! molar flux of H2SO4(g) [mol/m3/s]
      REAL( 8 ) :: CH2SO4          ! effective H2SO4(g) concentration [mol/m3]
      INTEGER      ISP             ! inorganic species index

! Adjoint variables
      REAL( 8 ), INTENT( OUT ) :: GNH3R8_ADJ
      REAL( 8 ), INTENT( OUT ) :: GNO3R8_ADJ
      REAL( 8 ), INTENT( OUT ) :: GCLR8_ADJ
      REAL( 8 ), INTENT( INOUT ) :: CEQ_ADJ ( NVOLINORG )
      REAL( 8 ), INTENT( OUT ) :: CONDRATE_ADJ
      REAL( 8 ), INTENT( OUT ) :: HPLUS_ADJ
      REAL( 8 ), INTENT( OUT ) :: RATE_ADJ
      REAL( 8 )                :: J_ADJ ( NVOLINORG )
      REAL( 8 ) :: CINF_ADJ ( NVOLINORG )
      REAL( 8 ) :: QK_ADJ
      REAL( 8 ) :: HFLUX_ADJ
      REAL( 8 ) :: HLIM_ADJ
      REAL( 8 ) :: AA_ADJ
      REAL( 8 ) :: BB_ADJ
      REAL( 8 ) :: CC_ADJ
      REAL( 8 ) :: JH2SO4_ADJ
      REAL( 8 ) :: CH2SO4_ADJ

      ! Checkpointing Variables
      REAL( 8 ) :: CEQ_CHECK( NVOLINORG )

! initialize adjoint variables
      GNH3R8_ADJ = 0.0
      GNO3R8_ADJ = 0.0
      GCLR8_ADJ = 0.0
      CONDRATE_ADJ = 0.0
      HPLUS_ADJ = 0.0
      RATE_ADJ = 0.0
      CINF_ADJ = 0.0
      QK_ADJ = 0.0
      HFLUX_ADJ = 0.0
      HLIM_ADJ  = 0.0
      AA_ADJ = 0.0
      BB_ADJ = 0.0    
      CC_ADJ = 0.0
      JH2SO4_ADJ = 0.0
      CH2SO4_ADJ = 0.0

!---------------------------- Begin Execution -----------------------

! Forward Code:
C     Convert gas concentration from ug/m3 to mol/m3
      Cinf( KNH4 ) = GNH3R8 * 1.0D-6 / PRECURSOR_MW( NH3_IDX )
      Cinf( KNO3 ) = GNO3R8 * 1.0D-6 / PRECURSOR_MW( HNO3_IDX )
      Cinf( KCL )  = GCLR8  * 1.0D-6 / PRECURSOR_MW( HCL_IDX )
            
C     Calculate cond/evap fluxes (no H+ limiting)
      DO isp = 1, nvolinorg
         J( isp ) = CondRate * ( Cinf( isp ) - Ceq( isp ) )
      END DO

C     Convert rate to mol/m3/s and get effective Cinf for H2SO4(g)
      JH2SO4  = rate * 1.0D-6 / PRECURSOR_MW( SULPRD_IDX )
      CH2SO4  = JH2SO4 / CondRate

C     Limit H+ flux (Pilinis et al., 2000, AS&T). Note: J is flux
C     to entire mode, not one particle
      Hlim  = Afact * Hplus
      Hflux = 2.0D0 * JH2SO4 + J( KNO3 ) + J( KCL ) - J( KNH4 )

C     If Hflux is too large, limit the flux by modifying species
C     vapor pressures with Qk factor (Pilinis et al., 2000, AS&T).
      IF ( ABS( Hflux ) .GT. Hlim ) THEN
         Hlim = SIGN( Hlim, Hflux )

C        Solve quadratic for Qk: aa*Qk^2 + bb*Qk + cc = 0
         aa = Ceq( KCL ) + Ceq( KNO3 )

         bb = Hlim / CondRate
     &      + Cinf( KNH4) - Cinf( KNO3 ) - Cinf( KCL ) - 2.0D0 * CH2SO4
         cc = -Ceq( KNH4 )

         Qk = 0.0D0 ! initialize Qk

         IF ( aa .LT. small .AND. 0.0D0 .LT. bb ) THEN ! bb*Qk + cc = 0
            Qk = -cc / bb
         ELSE IF (aa .LT. small .AND. bb .LE. 0.0D0 ) THEN
            Qk = 0.0D0
         ELSE IF (-cc .LT. small .AND. bb .LT. 0.0D0 ) THEN  ! aa*Qk^2 + bb*Qk = 0
            Qk = -bb / aa
         ELSE IF (-cc .LT. small .AND. 0.0D0 .LE. bb ) THEN
            Qk = 0.0D0
         ELSE 
            Qk = ( -bb + SQRT ( bb**2 - 4.0D0 * aa * cc ) ) / ( 2.0D0 * 
     &             aa )
            IF ( bb ** 2 - 4.0D0 * aa * cc .LT. 0.0D0 ) THEN
               PRINT *, 'Compute_Flux, sqrt<0'
               Qk = 0.0D0
            END IF
         END IF

C     Modify vapor pressures and get new fluxes
         IF ( Qk .GT. small ) THEN
            Ceq( KNH4 ) = Ceq( KNH4 ) / Qk
            Ceq( KNO3 ) = Ceq( KNO3 ) * Qk
            Ceq( KCl )  = Ceq( KCl )  * Qk

            ! Checkpointing
            CEQ_CHECK( KNH4 ) = CEQ( KNH4 )
            CEQ_CHECK( KNO3 ) = CEQ( KNO3 )
            CEQ_CHECK( KCL  ) = CEQ( KCL  )
            
            DO isp = 1, nvolinorg
               J( isp ) = CondRate * ( Cinf( isp ) - Ceq( isp ) )
            END DO
         END IF

      END IF   ! |Hflux| > Hlim

! Adjoint Code:

      IF ( ABS( HFLUX ) .GT. HLIM ) THEN

!     Modify vapor pressures and get new fluxes
         IF ( QK .GT. SMALL ) THEN
            !------
            ! fwd code:
            ! DO isp = 1, nvolinorg
            !    J( isp ) = CondRate * ( Cinf( isp ) - Ceq( isp ) )
            ! END DO
            ! adj code:
            DO ISP = NVOLINORG, 1, -1
               CONDRATE_ADJ = CONDRATE_ADJ + ( CINF ( ISP ) - CEQ( ISP )
     &                      ) * J_ADJ ( ISP )
               CINF_ADJ ( ISP ) = CINF_ADJ ( ISP ) + CONDRATE * J_ADJ (
     &                            ISP )
               CEQ_ADJ ( ISP ) = CEQ_ADJ ( ISP ) - CONDRATE * J_ADJ (
     &                           ISP )
               J_ADJ ( ISP ) = 0.0
            END DO

            ! Checkpointing
            CEQ( KNH4 ) = CEQ_CHECK( KNH4 )
            CEQ( KNO3 ) = CEQ_CHECK( KNO3 )
            CEQ( KCL ) = CEQ_CHECK( KCL )

            !------
            ! fwd code:
            ! Ceq( KNH4 ) = Ceq( KNH4 ) / Qk
            ! Ceq( KNO3 ) = Ceq( KNO3 ) * Qk
            ! Ceq( KCl )  = Ceq( KCl )  * Qk
            ! adj code:
            QK_ADJ = QK_ADJ + CEQ( KCL ) * CEQ_ADJ ( KCL )
            CEQ_ADJ ( KCL ) = CEQ_ADJ ( KCL ) * QK 
            QK_ADJ = QK_ADJ + CEQ( KNO3 ) * CEQ_ADJ ( KNO3 )
            CEQ_ADJ ( KNO3 ) = CEQ_ADJ ( KNO3 ) * QK
            QK_ADJ = QK_ADJ - ( CEQ( KNH4 ) / ( QK ** 2.0 ) ) * CEQ_ADJ
     &               ( KNH4 )
            CEQ_ADJ ( KNH4 ) = CEQ_ADJ ( KNH4 ) / QK

         END IF

         IF ( AA .LT. SMALL .AND. 0.0D0 .LT. BB ) THEN ! bb*Qk + cc = 0
            !------
            ! fwd code:
            ! Qk = -cc / bb
            ! adj code:
            CC_ADJ = CC_ADJ - ( 1.0 / BB ) * QK_ADJ
            BB_ADJ = BB_ADJ + ( CC / ( BB ** 2.0 ) ) * QK_ADJ
            QK_ADJ = 0.0
         ELSE IF (AA .LT. SMALL .AND. BB .LE. 0.0D0 ) THEN
            !------
            ! fwd code:
            ! Qk = 0.0D0
            ! adj code:
            QK_ADJ = 0.0
         ELSE IF (-CC .LT. SMALL .AND. BB .LT. 0.0D0 ) THEN  ! aa*Qk^2 + bb*Qk = 0
            !------
            ! fwd code:
            ! Qk = -bb / aa
            ! adj code:
            BB_ADJ = BB_ADJ - ( 1.0 / AA ) * QK_ADJ
            AA_ADJ = AA_ADJ + ( BB / ( AA ** 2.0 ) ) * QK_ADJ
            QK_ADJ = 0.0
         ELSE IF (-CC .LT. SMALL .AND. 0.0D0 .LE. BB ) THEN
            !------
            ! fwd code:
            ! Qk = 0.0D0
            ! adj code:
            QK_ADJ = 0.0
         ELSE
            IF ( BB ** 2 - 4.0D0 * AA * CC .LT. 0.0D0 ) THEN
               !------
               ! fwd code:
               ! Qk = 0.0D0
               ! adj code:
               QK_ADJ = 0.0
            END IF
            !------
            ! fwd code:
            ! Qk = ( -bb + SQRT ( bb**2 - 4.0D0 * aa * cc ) ) / ( 2.0D0 * aa )
            ! adj code:
            AA_ADJ = AA_ADJ - ( ( CC / ( AA * SQRT( BB ** 2.0 - 4.0 * AA
     &                      * CC ) ) ) - ( 1.0 / 2.0 ) * ( ( SQRT( BB **
     &                        2.0 - 4.0 * AA * CC ) - BB ) / ( AA ** 2.0
     &                      ) ) ) * QK_ADJ
            BB_ADJ = BB_ADJ + ( ( 1.0 / 2.0 ) * ( ( BB / SQRT( BB ** 2.0
     &                      - 4.0 * AA * CC ) ) - 1.0 ) / AA ) * QK_ADJ
            CC_ADJ = CC_ADJ - ( 1.0 / SQRT( BB ** 2.0 - 4.0 * AA * CC )
     &             ) * QK_ADJ
            QK_ADJ = 0.0

         END IF

         !------
         ! fwd code:
         ! Qk = 0.0D0 ! initialize Qk
         ! adj code:
         QK_ADJ = 0.0

         !------
         ! fwd code:
         ! aa = Ceq( KCL ) + Ceq( KNO3 )
         ! bb = Hlim / CondRate
         ! &      + Cinf( KNH4) - Cinf( KNO3 ) - Cinf( KCL ) - 2.0D0 * CH2SO4
         ! cc = -Ceq( KNH4 )
         ! adj code:
         CEQ_ADJ ( KNH4 ) = CEQ_ADJ ( KNH4 ) - CC_ADJ
         CC_ADJ = 0.0
         HLIM_ADJ = HLIM_ADJ + ( 1.0D0 / CONDRATE ) * BB_ADJ
         CONDRATE_ADJ = CONDRATE_ADJ - ( HLIM / ( CONDRATE ** 2.0 ) ) * 
     &                  BB_ADJ
         CINF_ADJ ( KNH4 ) = CINF_ADJ ( KNH4 ) + BB_ADJ
         CINF_ADJ ( KNO3 ) = CINF_ADJ ( KNO3 ) - BB_ADJ
         CINF_ADJ ( KCL ) = CINF_ADJ ( KCL ) - BB_ADJ
         CH2SO4_ADJ = CH2SO4_ADJ - 2.0 * BB_ADJ
         BB_ADJ = 0.0
         CEQ_ADJ ( KCL ) = CEQ_ADJ ( KCL ) + AA_ADJ
         CEQ_ADJ ( KNO3 ) = CEQ_ADJ ( KNO3 ) + AA_ADJ
         AA_ADJ = 0.0

         !------
         ! fwd code:
         ! Hlim = SIGN( Hlim, Hflux )
         ! adj code:
         HLIM_ADJ = SIGN( 1.d0, HLIM * HFLUX ) * HLIM_ADJ

      END IF   ! |Hflux| > Hlim

      !------
      ! fwd code:
      ! Hlim  = Afact * Hplus
      ! Hflux = 2.0D0 * JH2SO4 + J( KNO3 ) + J( KCL ) - J( KNH4 )
      ! adj code:
      HPLUS_ADJ = HPLUS_ADJ + AFACT * HLIM_ADJ
      HLIM_ADJ = 0.0

      !------
      ! fwd code:
      ! JH2SO4  = rate * 1.0D-6 / PRECURSOR_MW( SULPRD_IDX )
      ! CH2SO4  = JH2SO4 / CondRate
      ! adj code:
      JH2SO4_ADJ = JH2SO4_ADJ + ( 1.0 / CONDRATE ) * CH2SO4_ADJ
      CONDRATE_ADJ = CONDRATE_ADJ - ( JH2SO4 / ( CONDRATE ** 2.0 ) ) *
     &               CH2SO4_ADJ
      CH2SO4_ADJ = 0.0
      RATE_ADJ = RATE_ADJ + ( 1.0D-6 / PRECURSOR_MW( SULPRD_IDX ) )
     &         * JH2SO4_ADJ
      JH2SO4_ADJ = 0.0

      !------
      ! fwd code:
      ! DO isp = 1, nvolinorg
      !    J( isp ) = CondRate * ( Cinf( isp ) - Ceq( isp ) )
      ! END DO
      ! adj code:
      DO ISP = NVOLINORG, 1, -1
         CONDRATE_ADJ = CONDRATE_ADJ + ( CINF( ISP ) - CEQ( ISP ) ) *
     &                  J_ADJ ( ISP )
         CINF_ADJ ( ISP ) = CINF_ADJ ( ISP ) + CONDRATE * J_ADJ ( ISP )
         CEQ_ADJ ( ISP ) = CEQ_ADJ ( ISP ) - CONDRATE * J_ADJ ( ISP )
      END DO

      !------
      ! fwd code:
      ! Cinf( KNH4 ) = GNH3R8 * 1.0D-6 / PRECURSOR_MW( NH3_IDX )
      ! Cinf( KNO3 ) = GNO3R8 * 1.0D-6 / PRECURSOR_MW( HNO3_IDX )
      ! Cinf( KCL )  = GCLR8  * 1.0D-6 / PRECURSOR_MW( HCL_IDX )
      ! adj code:
      GCLR8_ADJ = GCLR8_ADJ + ( 1.0D-6 / PRECURSOR_MW( HCL_IDX ) ) *
     &            CINF_ADJ ( KCL )
      CINF_ADJ ( KCL ) = 0.0
      GNO3R8_ADJ = GNO3R8_ADJ + ( 1.0D-6 / PRECURSOR_MW( HNO3_IDX ) ) * 
     &             CINF_ADJ ( KNO3 ) 
      CINF_ADJ ( KNO3 ) = 0.0
      GNH3R8_ADJ = GNH3R8_ADJ + ( 1.0D-6 / PRECURSOR_MW( NH3_IDX ) ) *
     &             CINF_ADJ ( KNH4 )
      CINF_ADJ ( KNH4 ) = 0.0

      END SUBROUTINE COMPUTE_FLUX_ADJ
