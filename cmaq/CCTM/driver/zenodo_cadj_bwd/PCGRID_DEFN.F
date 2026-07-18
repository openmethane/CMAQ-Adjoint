
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_bwd/PCGRID_DEFN.F,v 1.1 2010/11/17 19:22:51 sergey Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE PCGRID_DEFN

C-----------------------------------------------------------------------
C Function:
C   Create CGRID target, PCGRID. CGRID becomes a pointer to enable a one-
C   cell ghost region around CGRID

C Preconditions:
C   Horizontal domain extents must be set (subroutine PAR_INIT -> HGRD_DEFN)
C   Number of species in the species groups must be available (include files
C   in CGRID_SPCS)
C   This module is referenced in driver and sciproc. driver must call HGRD_INIT,
C   CGRID_MAP, and PAR_INIT prior to calling PCGRID_INIT (below).

C Subroutines and functions called:

C Revision history:
C   Apr 00 - Al, David, Jeff - initial implementation
C   Dec 00 - Jeff - add CGRID_MAP "method"
C   Jan 01 - Jeff - remove CGRID_MAP "method" - now in module CGRID_SPCS
C   Nov 01 - Jeff - correct pointer reference for ping
C   Dec 04 - Jeff - dyn alloc - USE VGRD_DEFN
C   Jan 05 - Jeff - dyn alloc - establish both horizontal & vertical
C                   domain specifications in one module

C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications
      USE CGRID_SPCS          ! CGRID species number and offsets
      USE EMIS_DEFN, ONLY: EMLAYS
      USE PMEM_DEFN, ONLY: N_EMIS_PM

      IMPLICIT NONE

      INCLUDE SUBST_GC_EMIS   ! gas chemistry species table
      INCLUDE SUBST_AE_EMIS   ! aerosol species table
      INCLUDE SUBST_NR_EMIS   ! non-reactive species table
      INCLUDE SUBST_TR_EMIS   ! tracer species table

      INTEGER :: N_SPC_PEMIS
      REAL, ALLOCATABLE, SAVE, TARGET :: PCGRID( :,:,:,: )
      REAL, ALLOCATABLE, SAVE, TARGET :: PLGRID( :,:,:,: )
      REAL, ALLOCATABLE, SAVE, TARGET :: PLGRID_EM( :,:,:,: )

      CONTAINS
         FUNCTION PCGRID_INIT () RESULT ( SUCCESS )

!        INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
         INCLUDE SUBST_IODECL      ! I/O definitions and declarations

         LOGICAL :: SUCCESS
         INTEGER :: ALLOCSTAT
         INTEGER, SAVE :: LOGDEV
         LOGICAL, SAVE :: FIRSTIME = .TRUE.
         CHARACTER( 120 ) :: XMSG = ' '

C PCGRID extents
         INTEGER :: LCOL
         INTEGER :: HCOL
         INTEGER :: LROW
         INTEGER :: HROW
         CHARACTER( 16 ), SAVE :: CTM_EMLAYS = 'CTM_EMLAYS'
         INTEGER :: STATUS
         INTEGER, EXTERNAL :: ENVINT

         SUCCESS = .TRUE.

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            LOGDEV = INIT3()

!slz read emlays from environment
            EMLAYS = ENVINT( CTM_EMLAYS, 'Number of Emissions Layers', NLAYS, STATUS )

#ifdef parallel
            LCOL = 1 - NTHIK
            HCOL = NCOLS + NTHIK
            LROW = 1 - NTHIK
            HROW = NROWS + NTHIK
#else
            LCOL = 1
            HCOL = NCOLS
            LROW = 1
            HROW = NROWS
#endif

            N_SPC_PEMIS = N_EMIS_PM
     &                  + N_GC_EMIS   ! = N_GC_SPC + 1
     &                  + N_NR_EMIS
     &                  + N_TR_EMIS
            
            ALLOCATE ( PCGRID( LCOL:HCOL,LROW:HROW,NLAYS,NSPCSD ),
     &                 STAT = ALLOCSTAT )
            ALLOCATE ( PLGRID( LCOL:HCOL,LROW:HROW,NLAYS,NSPCSD ),
     &                 STAT = ALLOCSTAT )
            ALLOCATE ( PLGRID_EM( LCOL:HCOL,LROW:HROW,EMLAYS,N_SPC_PEMIS ),
     &                 STAT = ALLOCSTAT )

            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating PCGRID'
               CALL M3WARN ( 'PCGRID_INIT', 0, 0, XMSG )
               SUCCESS = .FALSE.; RETURN
               END IF

            PCGRID = 0.0
            PLGRID = 0.0
            PLGRID_EM = 0.0

            ELSE
            XMSG = 'PCGRID already ALLOCATED'
            CALL M3WARN ( 'PCGRID_INIT', 0, 0, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF

          RETURN
          END FUNCTION PCGRID_INIT

      END MODULE PCGRID_DEFN
