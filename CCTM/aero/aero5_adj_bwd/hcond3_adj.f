
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
      SUBROUTINE HCOND3_ADJ   ( AM0, AM0_ADJ, AM1, AM1_ADJ, AM2 ,AM2_ADJ
     &                        , DV, ALPHA, CBAR, F, F_ADJ )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of subroutine that calculates the size-dependent term in the
!   condensational-growth rate expression for the 2nd and 3rd moments
!   of a lognormal aerosol mode using the harmonic mean method.
!
! INPUTS:
!   AM0, AM1, DV, ALPHA, CBAR, F_ADJ
!
! OUTPUTS:
!   AM0_ADJ, AM1_ADJ
!
! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      IMPLICIT NONE

! *** Arguments:

      REAL( 8 ), INTENT( IN ) :: am0   ! zeroth moment of mode  [ #/m**3 ]
      REAL( 8 ), INTENT( IN ) :: am1   ! first moment of mode   [ m/m**3 ]
      REAL( 8 ), INTENT( IN ) :: am2   ! second moment of mode  [ m**2/m**3 ]
      REAL,      INTENT( IN ) :: DV    ! molecular diffusivity of the
                                       ! condensing vapor  [ m**2/s ]
      REAL,      INTENT( IN ) :: alpha ! accommodation coefficient
      REAL,      INTENT( IN ) :: cbar  ! kinetic velocity of condensing vapor [ m/s ]

      REAL( 8 )                :: F( 2 ) ! size-dependent term in condensational-growth
                                         ! rate: F(1) = 2nd moment [ m**2/m**3 s ]
                                         !       F(2) = 3rd moment [ m**3/m**3 s ]

C *** Local Variables:

      REAL( 8 ) :: GNC2 ! integrals used to calculate F(1) [m^2 / m^3 s]
      REAL( 8 ) :: GFM2 !

      REAL( 8 ) :: GNC3 ! integrals used to calculate F(2) [m^3 / m^3 s]
      REAL( 8 ) :: GFM3 !

      REAL( 8 ), PARAMETER :: pi = 3.14159265358979324D0
      REAL( 8 ), PARAMETER :: twopi = 2.0D0 * pi
      REAL( 8 ), PARAMETER :: pi4 = 0.25D0 * pi

! Adjoint variables

      REAL( 8 )                 :: AM0_ADJ
      REAL( 8 )                 :: AM1_ADJ
      REAL( 8 )                 :: AM2_ADJ
      REAL( 8 )                 :: F_ADJ ( 2 )
      REAL( 8 )                 :: GNC2_ADJ, GNC3_ADJ, GFM2_ADJ
      REAL( 8 )                 :: GFM3_ADJ

      CHARACTER( 80 ) :: filename, jdate_str, jtime_str
      CHARACTER( 80 ) :: filename1, filename2, i_str, counter_str
      LOGICAL :: exists

! Initialize adjoint variables
      GNC2_ADJ = 0.0
      GNC3_ADJ = 0.0
      GFM2_ADJ  = 0.0
      GFM3_ADJ = 0.0

!------------------------ Begin Execution ----------------------

! Forward Code:

C *** Implement equation A15 of Binkowski & Shankar (1995) for the
C     2nd and 3rd moments of a lognormal mode of arbitrary size.

      gnc2 = twopi * Dv * am0          ! 2nd moment, near-continuum
      gnc3 = twopi * Dv * am1          ! 3rd moment, near-continuum
      gfm2 = pi4 * alpha * cbar * am1  ! 2nd moment, free-molecular
      gfm3 = pi4 * alpha * cbar * am2  ! 3rd moment, free-molecular

C *** Implement equation A13 of Binkowski & Shankar (1995) for a
C     lognormal mode of arbitrary size.  These are the size-dependent
C     terms in the condensational-growth rate expression, given in
C     equation 7a of B&S (1995).

      F( 1 ) = gnc2 * gfm2 / ( gnc2 + gfm2 )  ! 2nd moment
      F( 2 ) = gnc3 * gfm3 / ( gnc3 + gfm3 )  ! 3rd moment

! Adjoint Code:

      !------
      ! fwd code:
      ! F( 1 ) = gnc2 * gfm2 / ( gnc2 + gfm2 )  ! 2nd moment
      ! F( 2 ) = gnc3 * gfm3 / ( gnc3 + gfm3 )  ! 3rd moment
      ! adj code:
      GNC3_ADJ = ( ( GFM3 / ( GFM3 + GNC3 ) ) - ( ( GFM3 * 
     &           GNC3 ) / ( ( GFM3 + GNC3 ) ** 2.0 ) ) ) * F_ADJ ( 2 )
      GFM3_ADJ = ( ( GNC3 / ( GFM3 + GNC3 ) ) - ( ( GFM3 * 
     &           GNC3 ) / ( ( GFM3 + GNC3 ) ** 2.0 ) ) ) * F_ADJ ( 2 )
      F_ADJ ( 2 ) = 0.0
      GNC2_ADJ = ( ( GFM2 / ( GFM2 + GNC2 ) ) - ( ( GFM2 *
     &           GNC2 ) / ( ( GFM2 + GNC2 ) ** 2.0 ) ) ) * F_ADJ ( 1 )
      GFM2_ADJ = ( ( GNC2 / ( GFM2 + GNC2 ) ) - ( ( GFM2 * 
     &           GNC2 ) / ( ( GFM2 + GNC2 ) ** 2.0 ) ) ) * F_ADJ ( 1 )
      F_ADJ ( 1 ) = 0.0

      !------
      ! fwd code:
      ! gnc2 = twopi * Dv * am0          ! 2nd moment, near-continuum
      ! gnc3 = twopi * Dv * am1          ! 3rd moment, near-continuum
      ! gfm2 = pi4 * alpha * cbar * am1  ! 2nd moment, free-molecular
      ! gfm3 = pi4 * alpha * cbar * am2  ! 3rd moment, free-molecular
      ! adj code:
      AM2_ADJ = AM2_ADJ + PI4 * ALPHA * CBAR * GFM3_ADJ 
      GFM3_ADJ = 0.0
      AM1_ADJ = AM1_ADJ + PI4 * ALPHA * CBAR * GFM2_ADJ
      GFM2_ADJ = 0.0
      AM1_ADJ = AM1_ADJ + TWOPI * DV * GNC3_ADJ
      GNC3_ADJ = 0.0
      AM0_ADJ = AM0_ADJ + TWOPI * DV * GNC2_ADJ
      GNC2_ADJ = 0.0

      END SUBROUTINE HCOND3_ADJ
