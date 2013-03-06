C-----------------------------------------------------------------------
      Module met_data

C /////////////////////////////////////////////////////////////////////
C  Module met_data defines the parameters and meteorological data 
C     required in aerosol calculations.
C
C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C----------------------------------------------------------------------

      Implicit None

C fundamental constants
      Include SUBST_CONST    ! common constants

      Real( 8 ), Parameter :: f6dpi = 6.0D0/ pi
      Real( 8 ), Parameter :: pi6 = pi / 6.0D0
      Real( 8 ), Parameter :: boltzmann = rgasuniv / avo  ! [ j/K ]

C concentration lower limit used for gas species                                  
      Real,      Parameter :: min_gasconc = 1.0E-30   ! [ kg/m^3 ]

C meteorological data (from MET_CRO_3D file)
      Real :: airtemp      ! air temp [ K ]
      Real :: airprs       ! air pressure [ Pa]
      Real :: airdens      ! air density [ kg/m^3]
      Real :: airRH        ! fractional relative humidity
      Real :: airQV        ! water vapor mass mixing ratio [ kg_H2O/kg_air ]
      Real :: H2OsatVP     ! saturation water vapor pressure [ Pa ]
      Real :: H2OVP        ! ambient water vapor pressure [ Pa ]

      End Module met_data
