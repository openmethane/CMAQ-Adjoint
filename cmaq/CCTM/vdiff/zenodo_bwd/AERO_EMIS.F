
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_bwd/AERO_EMIS.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  MODULE AERO_EMIS contains emissions code required for the modal
C     aerosol module in CMAQ
C                                 Coded by Dr. Francis S. Binkowski
C                                      and Dr. Jeffrey O. Young
C
C  CONTAINS: SUBROUTINE RDEMIS_AE
C            Variable declarations needed for other subroutines in CMAQ's
C             vertical diffusion module
C
C  DEPENDENT UPON:  NONE
C
C  REVISION HISTORY:
C
C   30 Aug 01 J.Young:  dyn alloc - Use HGRD_DEFN
C   09 Oct 03 J.Gipson: added MW array for AE emis species to module contents
C   31 Jan 05 J.Young:  dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module, GRID_CONF
C   26 Apr 05 P.Bhave:  removed code supporting the "old type" of emission
C                        files that had unspeciated PM10 and PM2.5 only
C                       removed need for 'AERO_SPC.EXT' by declaring the
C                        required variables locally
C   13 Jun 05 P.Bhave:  added vars needed for sea-salt emission processing
C                       inherit N_AE_EMIS,AE_EMIS,AE_EMIS_MAP from AE_EMIS.EXT
C                       moved RHO* parameters from RDEMIS_AE to this module
C                        for use by SSEMIS routine
C   24 Aug 07 J.Young:  Modified to enable in-line plume rise calculation for
C                       3D pt source emissions. Distinguish between PM (primary,
C                       unspeciated, file data) and AE (model speciated). Re-
C                       named RDEMIS_AE to GET_AERO_EMIS.
C   11 Apr 08 J.Kelly:  added code to emit coarse surface area
C   Feb 2011 by Peter Percell at UH-IMAQS:
C   -- reduced to essentials for adjoint/4dvar backward mode
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      MODULE AERO_EMIS

      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE
      
      INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names and map table
            
C Array dimensions
      INTEGER, PARAMETER :: NSSDIAG = 12   ! number of species in sea-salt
                                           ! diagnostic emission file
      INTEGER, PARAMETER :: NSSSPC = 4     ! number of chemical species in
                                           ! fresh sea-salt aerosol
      INTEGER, PARAMETER :: NSSMOD = 3     ! number of lognormal modes in
                                           ! sea-salt aerosol


C Molar masses (molecular wt.) of each aerosol species
      REAL, SAVE :: AE_EM_MW( N_AE_EMIS )
      REAL, SAVE :: R_AE_EM_MW( N_AE_EMIS )

C     Bulk component densities [ kg/m3 ]
      REAL, PARAMETER :: RHOSO4  = 1.8E3  ! density of sulfate aerosol
      REAL, PARAMETER :: RHONO3  = 1.8E3  ! density of nitrate aerosol
      REAL, PARAMETER :: RHOORG  = 2.0E3  ! density of organic aerosol
      REAL, PARAMETER :: RHOSOIL = 2.6E3  ! density of soil dust
      REAL, PARAMETER :: RHOSEAS = 2.2E3  ! density of marine aerosol
      REAL, PARAMETER :: RHOANTH = 2.2E3  ! density of elemental carbon and
                                          !  unspeciated anthropogenic aerosol
      REAL, PARAMETER :: RHOH2O  = 1.0E3  ! density of water
                                             
C Aerosol species names in sea-salt output file
      CHARACTER( 16 ), SAVE :: WRSS_SPC( NSSDIAG )  ! in sea-salt output file

C Species names in the sea-salt-emissions diagnostic file
      DATA WRSS_SPC(  1 ) / 'ANAJ'  /   ! accumulation mode sodium
      DATA WRSS_SPC(  2 ) / 'ACLJ'  /   ! accumulation mode chloride
      DATA WRSS_SPC(  3 ) / 'ASO4J' /   ! accumulation mode sulfate
      DATA WRSS_SPC(  4 ) / 'AH2OJ' /   ! accumulation mode water
      DATA WRSS_SPC(  5 ) / 'ANAK'  /   ! coarse mode sodium
      DATA WRSS_SPC(  6 ) / 'ACLK'  /   ! coarse mode chloride
      DATA WRSS_SPC(  7 ) / 'ASO4K' /   ! coarse mode sulfate 
      DATA WRSS_SPC(  8 ) / 'AH2OK' /   ! coarse mode water   
      DATA WRSS_SPC(  9 ) / 'ANUMJ' /   ! accumulation mode number
      DATA WRSS_SPC( 10 ) / 'ANUMK' /   ! coarse mode number
      DATA WRSS_SPC( 11 ) / 'ASRFJ' /   ! accumulation mode surface area
      DATA WRSS_SPC( 12 ) / 'ASRFK' /   ! coarse mode surface area

C Indices to sea-salt species in mass-emission arrays
      INTEGER, PARAMETER :: KNA  = 1   ! position of sodium
      INTEGER, PARAMETER :: KCL  = 2   ! position of chloride
      INTEGER, PARAMETER :: KSO4 = 3   ! position of sulfate
      INTEGER, PARAMETER :: KH2O = 4   ! position of water

      END MODULE AERO_EMIS
