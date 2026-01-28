
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

      INCLUDE SUBST_AE_EMIS     ! aerosol emission surrogate names and map table

C Molar masses (molecular wt.) of each aerosol species
      REAL, SAVE :: AE_EM_MW( N_AE_EMIS )

      END MODULE AERO_EMIS
