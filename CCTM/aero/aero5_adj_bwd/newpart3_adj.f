
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
      SUBROUTINE NEWPART3_ADJ   ( RH, TEMP, XH2SO4, XH2SO4_ADJ, SO4RATE,
     &                          SO4RATE_ADJ, DNDT, DNDT_ADJ, DMDT_SO4, 
     &                          DMDT_SO4_ADJ, DM2DT, DM2DT_ADJ )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of subroutine that calculates the new particle production
!   rate due to binary nucleation of H2SO4 and H2O vapor.
!
! INPUTS:
!   RH, Temp, XH2SO4, SO4RATE, DNDT_ADJ, DMDT_SO4_ADJ, DM2DT_ADJ
!
! OUTPUTS:
!   XH2SO4_ADJ, SO4RATE_ADJ
!
! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      USE AERO_DATA_B
      USE MET_DATA

      IMPLICIT NONE

! *** Arguments:

      REAL, INTENT( IN ) :: RH         ! fractional relative humidity
      REAL, INTENT( IN ) :: TEMP       ! ambient temperature [ K ]
      REAL, INTENT( IN ) :: XH2SO4     ! sulfuric acid concentration [ ug/m**3 ]
      REAL, INTENT( IN ) :: SO4RATE    ! gas-phase H2SO4 production rate [ ug/m**3 s ]

      REAL( 8 )                :: DNDT     ! particle number production rate [ m^-3/s ]
      REAL( 8 )                :: DMDT_so4 ! SO4 mass production rate [ ug/m**3 s ]
      REAL( 8 )                :: DM2DT    ! second moment production rate [ m**2/m**3 s ]

! *** Parameters
      
      CHARACTER( 16 ), PARAMETER :: PNAME = 'NEWPART_ADJ'

! *** Conversion Factors:

      REAL, PARAMETER :: MWH2SO4 = 98.0 ! molecular weight for H2SO4

      REAL, PARAMETER :: UGM3_NCM3 = (AVO / MWH2SO4) * 1.0e-12 ! micrograms/m**3 
! *** Particle size parameters:

      REAL, PARAMETER :: D20 = 2.0E-07            ! diameter of a new particle [cm]

      REAL, PARAMETER :: D20SQ = D20 * D20        ! new-particle diameter squared [cm**2]

      REAL, PARAMETER :: M2_20 = 1.0E-4 * D20SQ   ! new-particle diameter squared [m**2]

      REAL, PARAMETER :: V20 = PI6 * D20 * D20SQ  ! volume of a new particle [cm**3]

      REAL( 8 ) :: SULFMASS                       ! mass of a new particle [ug]
      REAL( 8 ) :: SULFMASS1                      ! inverse of sulfmass [ug**-1]

! *** Local Variables used to determine nucleation rate

      REAL NWV     ! water vapor concentration [ 1/cm**3 ]
      REAL NAV0    ! saturation vapor conc of H2SO4 [ 1/cm**3 ]
      REAL NAV     ! H2SO4 vapor concentration [ 1/cm**3 ]
      REAL NAC     ! critical H2SO4 vapor concentration [ 1/cm**3 ]
      REAL RA      ! fractional relative acidity
      REAL DELTA   ! temperature-correction term
      REAL XAL     ! H2SO4 mole fraction in the critical nucleus
      REAL NSULF   ! see usage
      REAL( 8 ) :: CHI   ! factor to calculate Jnuc
      REAL( 8 ) :: JNUCK ! nucleation rate [ cm ** -3  s ** -1 ]
                   ! (Kulmala et al.)
      REAL TT      ! dummy variable for statement functions

! *** Statement Functions

      REAL PH2SO4, PH2O ! arithmetic statement functions for saturation
                        ! vapor pressures of h2so4 and h2o [Pa] taken
                        ! from Appendix of Kulmala et al. (1998), p8306

      PH2O( TT )   = EXP( 77.34491296 -7235.424651 / TT
     &                   - 8.2 * LOG( TT ) + TT * 5.7113E-03 )

      PH2SO4( TT ) = EXP( 27.78492066 - 10156.0 / TT )

! Adjoint variables:

      REAL                :: XH2SO4_ADJ
      REAL                :: SO4RATE_ADJ
      REAL( 8 )               :: DNDT_ADJ
      REAL( 8 )               :: DMDT_SO4_ADJ
      REAL( 8 )               :: DM2DT_ADJ
      REAL :: NAV_ADJ
      REAL :: RA_ADJ
      REAL :: XAL_ADJ
      REAL :: NSULF_ADJ
      REAL( 8 ) :: CHI_ADJ
      REAL( 8 ) :: JNUCK_ADJ

! Initialize adjoint variables

      XH2SO4_ADJ = 0.0
      SO4RATE_ADJ = 0.0
      NAV_ADJ = 0.0
      RA_ADJ = 0.0
      XAL_ADJ = 0.0
      NSULF_ADJ = 0.0
      CHI_ADJ = 0.0
      JNUCK_ADJ = 0.0

!---------------------------- Begin Execution ---------------------

! Forward Code:

C *** Initialize variables
      DNDT     = 0.0D0
      DMDT_so4 = 0.0D0
      DM2DT    = 0.0D0

      sulfmass = 1.0D+3 * aerospc( ASO4_IDX )%density * v20
      sulfmass1 = 1.0D0 / sulfmass

C *** Nucleation Rate
C     Calculate the sulfuric acid/water nucleation rate.  This code
C     implements Section 3 of Kulmala et al. (1998).  Note that all
C     variables in the Kulmala parameterization are in cgs units.

C *** Calculate water vapor concentration [1/cm**3] using ambient RH
C     and the formula in Appendix of Kulmala et al. (1998), p8306.
      Nwv = RH * ph2o( Temp ) / ( RGASUNIV * Temp ) * AVO * 1.0E-6

C *** Calculate saturation vapor concentration of H2SO4 [1/cm**3]
C       using formula in the Appendix of Kulmala et al. (1998), p8306.
      Nav0 = ph2so4( Temp ) / ( RGASUNIV * Temp ) * AVO * 1.0E-6

C *** Convert ambient H2SO4 vapor concentration into [1/cm**3] units
      Nav = ugm3_ncm3 * XH2SO4

C *** Calculate critical concentration of H2SO4 vapor needed to produce
C     1 particle/(cm**3 s) using Equation 18 of Kulmala et al (1998)
      Nac = EXP( - 14.5125 + 0.1335 * Temp
     &           - 10.5462 * RH + 1958.4 * RH  / Temp )

C *** Calculate relative acidity, defined as the ambient concentration
C     divided by the saturation concentration of H2SO4 vapor
      RA = Nav / Nav0

C *** Calculate temperature correction factor using Equation 22 of
C      Kulmala et al (1998)
      delta = 1.0 + ( Temp - 273.15 ) / 273.15

C *** Calculate mole fraction of H2SO4 in the critical nucleus using
C     Equation 17 of Kulmala et al (1998)
      Xal = 1.2233 - 0.0154 * RA / ( RA + RH ) + 0.0102 * log( Nav )
     &    - 0.0415 * LOG( Nwv ) + 0.0016 * Temp

C *** Calculate Nsulf as defined in Equation 21 of Kulmala et al (1998)
      Nsulf = LOG( Nav / Nac )

C *** Calculate natural log of nucleation rate using Equation 20 of
C     Kulmala et al (1998)
      chi = 25.1289 * Nsulf - 4890.8 * Nsulf / Temp
     &    - 1743.3 / Temp - 2.2479 * delta * Nsulf * RH
     &    + 7643.4 * Xal / Temp
     &    - 1.9712 * Xal * delta / RH

C *** Calculate nucleation rate using Eq 19 of Kulmala et al (1998)
      JnucK = Exp( chi )   ! [ #/cm**3 s ]

C *** Moment Production Rates
C     Calculate the production rate of number, sulfate mass, and 2nd
C     moment, due to the H2SO4/H2O nucleation assuming that each critical
C     nucleus grows instantaneously to 2.0 nm.  This code follows Section
C     1.2 of Binkowski & Roselle (2003), except the assumed particle
C     diameter has been changed from 3.5 to 2.0 nm.

C *** Convert production rate of number to [ #/m**3 s]
      DNDT = 1.0E06 * JnucK

C *** Calculate mass production rate [ ug / (m**3 s) ] analogous to
C     Equation 6a of Binkowski & Roselle (2003).  Set the upper limit
C     of the mass production rate as the gas-phase production rate of
C     H2SO4, and adjust the number production rate accordingly.
      DMDT_so4 = sulfmass * DNDT
      IF ( DMDT_so4 .GT. SO4RATE ) THEN
         DMDT_so4 = SO4RATE
         DNDT = DMDT_SO4 * sulfmass1
      END IF

C *** Calculate the production rate of 2nd moment [ m**2 / (m**3 s) ]
C     This is similar to Equation 6c of Binkowski & Roselle (2003),
C     except the factor of PI is removed and the assumed particle
C     diameter is different.
      DM2DT = DNDT * M2_20
!
! Adjoint Code:

      !------
      ! fwd code:
      ! DM2DT = DNDT * m2_20
      ! adj code:
      DNDT_ADJ = DNDT_ADJ + M2_20 * DM2DT_ADJ

      !------
      ! fwd code:
      ! IF ( DMDT_so4 .GT. SO4RATE ) THEN
      !   DMDT_so4 = SO4RATE
      !   DNDT = DMDT_SO4 * sulfmass1
      !END IF
      ! adj code:
      IF ( DMDT_SO4 .GT. SO4RATE ) THEN
         DMDT_SO4_ADJ = DMDT_SO4_ADJ + SULFMASS1 * DNDT_ADJ 
         DNDT_ADJ = 0.0
         SO4RATE_ADJ = SO4RATE_ADJ + DMDT_SO4_ADJ
         DMDT_SO4_ADJ = 0.0
      ELSE
         SO4RATE_ADJ = 0.0
      END IF

      !------
      ! fwd code:
      ! DMDT_so4 = sulfmass * DNDT
      ! adj code:
      DNDT_ADJ = DNDT_ADJ + SULFMASS * DMDT_SO4_ADJ

      !------
      ! fwd code:
      ! DNDT = 1.0E06 * JnucK
      ! adj code:
      JNUCK_ADJ = 1.0E06 * DNDT_ADJ

      !------
      ! fwd code:
      ! JnucK = Exp( chi ) 
      ! adj code:
      CHI_ADJ = EXP( CHI ) * JNUCK_ADJ

      !------
      ! fwd code:
      ! chi = 25.1289 * Nsulf - 4890.8 * Nsulf / Temp
      ! &    - 1743.3 / Temp - 2.2479 * delta * Nsulf * RH
      ! &    + 7643.4 * Xal / Temp
      ! &    - 1.9712 * Xal * delta / RH
      ! adj code:
      NSULF_ADJ = ( 25.1289 - 2.2479 * DELTA * RH - ( 4890.8
     &          / TEMP ) ) * CHI_ADJ
      XAL_ADJ = ( ( 7643.4 / TEMP ) - ( 1.9712 * DELTA / RH ) 
     &        ) * CHI_ADJ

      !------
      ! fwd code:
      ! Nsulf = LOG( Nav / Nac )
      ! adj code:
      NAV_ADJ = ( 1.0 / NAV ) * NSULF_ADJ

      !------
      ! fwd code:
      ! Xal = 1.2233 - 0.0154 * RA / ( RA + RH ) + 0.0102 * log( Nav )
      ! &    - 0.0415 * LOG( Nwv ) + 0.0016 * Temp
      ! adj code:
      RA_ADJ = ( ( 0.0154 * RA / ( RA + RH ) ** 2.0 ) - (
     &         0.0154 / ( RA + RH ) ) ) * XAL_ADJ
      NAV_ADJ = NAV_ADJ + ( 0.0102 / NAV ) * XAL_ADJ

      !------
      ! fwd code:
      ! RA = Nav / Nav0
      ! adj code:
      NAV_ADJ = NAV_ADJ + ( 1.0 / NAV0 ) * RA_ADJ

      !------
      ! fwd code:
      ! Nav = ugm3_ncm3 * XH2SO4
      ! adj code:
      XH2SO4_ADJ = UGM3_NCM3 * NAV_ADJ

      !------
      ! fwd code:
      ! DNDT     = 0.0D0
      ! DMDT_so4 = 0.0D0
      ! DM2DT    = 0.0D0
      ! adj code:
      DM2DT_ADJ = 0.0
      DMDT_SO4_ADJ = 0.0
      DNDT_ADJ = 0.0

      END SUBROUTINE NEWPART3_ADJ
