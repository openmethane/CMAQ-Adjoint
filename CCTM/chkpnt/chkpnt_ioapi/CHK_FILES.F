!**********************************************************************
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


      MODULE CHK_FILES

!-----------------------------------------------------------------------
! Function:
!   The module for manipulation of checkpoint files
!   Abstract version, for real file manipulation calls one of the implementations
!
! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version from Peter Percell
!-----------------------------------------------------------------------

      IMPLICIT NONE

      CONTAINS

!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_OPEN ( CHK_ID, SDATE, STIME, CREATENEW, TSTEP )
!-----------------------------------------------------------------------
! Function:
!     Subroutine creates checkpoint file
!     The checkpoint file remains opened

! Subroutines and functions called:
!     CHK_FILE_CREATE_IMPL, CHK_FILE_OPEN_IMPL, M3EXIT

! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!   Oct  2011 by Peter Percell at UH-IMAQS:
!     - removed unused checkpoint files
!     - made ADJ_VDIFF_CHK a full copy of the CGRID array
!   Nov  2011 by M. Turner, CU Boulder:
!     - Added phot species to chem checkpoint file
!   May  2013 by Peter Percell at UH-IMAQS:
!     - changed ADJ_RHOJ_CHK to ADJ_VA_RHOJ_CHK and added ADJ_HA_RHOJ_CHK
!     - added TSTEP as an optional calling argument to pass in a computed
!       sync step
!   July 2013 by J.Resler, ICS:
!     - fixed bug in VDESC initialization for AERO and VDIFF files
!   Aug 2013 by Peter Percell at UH-IMAQS:
!     - reduced vdiff checkpoint file to just layer 1 of cgrid
!   Sep  2013 by ShunLiu Zhao, Carleton:
!     - Removed phot species from chem checkpoint file
!     - Added TSTEP to chem checkpoint file
!   Oct 2013 by J.Resler, ICS:
!     - fixes, cleaning
!   Mar 2014 by Peter Percell, UH:
!     - replaced ADJ_XADV_CHK and ADJ_YADV_CHK with just ADJ_HADV_CHK
!     - renamed ADJ_ZADV_CHK to ADJ_VADV_CHK
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      USE CHK_FILES_IMPL
      USE ADJ_ARRAY_MOD         ! adjoint arrays (may not need) (slc.1.2014)

      IMPLICIT NONE

!     Include Files:

      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_CONC     ! gas chem conc file species and map table
      INCLUDE SUBST_AE_CONC     ! aerosol conc file species and map table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_GC_DIFF     ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF     ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF     ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF     ! tracer diffusion species and map table
      INCLUDE SUBST_GC_EMIS     ! gas chem emis surrogate names and map table
      ! Necessary to checkpoint emissions sensitivities (mdt|slc.1.2014)
      INCLUDE SUBST_AE_EMIS     ! aerosol emis surrogate names and map table 
      INCLUDE SUBST_NR_EMIS     ! non-reactive emis surrogate names and map table
      INCLUDE SUBST_TR_EMIS     ! tracer emis surrogate names and map table
      ! End of additional includes for emissions sensitivities
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC     ! file header data structure

!     Input parameters
      INTEGER, INTENT(IN)       :: CHK_ID       ! checkpoint file identification
      INTEGER, INTENT(IN)       :: SDATE        ! starting date (YYYYDDD)
      INTEGER, INTENT(IN)       :: STIME        ! starting time (HHMMSS)
      LOGICAL, INTENT(IN)       :: CREATENEW    ! create new file or open the old one
      INTEGER, INTENT(IN), OPTIONAL :: TSTEP(2) ! time step vector (HHMMSS)
                                                !   TSTEP(1) = local output step
                                                !   TSTEP(2) = sync. step

!     Local variables:
      CHARACTER(len=NAMLEN_ADJ) :: PNAME = 'CHK_FILE_OPEN'
      CHARACTER( 512 )          :: XMSG = ' '

      INTEGER, SAVE             :: LOGDEV
      LOGICAL, SAVE             :: FIRSTIME = .TRUE.

      type(CHK_FILE),pointer    :: CHK          ! checkpoint file description structure
      INTEGER                   :: INDX, L, S

!-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         LOGDEV = INIT3 ()
         FIRSTIME = .FALSE.
      END IF

!     File description structure
      CHK => CHKS(CHK_ID)

!     Check if the file is already opened
      IF (CHKS(CHK_ID)%ISUSED .EQV. .TRUE.) THEN
         RETURN
      END IF

!     chk file is to be opened
      CHK%ISUSED = .TRUE.
      CHK%FNAME = trim(CHK_NAMES(CHK_ID))
      write(logdev,*) 'Opening file '//CHK%FNAME

      CHK%SDATE = SDATE
      CHK%STIME = STIME

      CHK%NLAYS = NLAYS         ! default for almost all checkpoint files

!  Special settings for particular file
      SELECT CASE ( CHK_ID )

         CASE ( ADJ_CHEM_CHK )
            CHK%FDESC(1) = 'Chemistry checkpoint file'
            CHK%FDESC(2) = 'From CMAQ model dyn alloc version CTM'
            CHK%NVARS = N_GC_SPC+1
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            CHK%VTYPE(1:CHK%NVARS) = M3REAL
            CHK%UNITS(1:N_GC_SPC) = 'ppmV'
            CHK%VNAME(1:N_GC_SPC) = GC_SPC(1:N_GC_SPC)
            CHK%VDESC(1:N_GC_SPC) = 'Variable '//CHK%VNAME(1:N_GC_SPC)
            CHK%VNAME(N_GC_SPC+1) = 'STEP'
            CHK%UNITS(N_GC_SPC+1) = 'min'
            CHK%VDESC(N_GC_SPC+1) = 'Time step'

         CASE ( ADJ_AERO_CHK )
            CHK%FDESC(1) = 'Aerosols checkpoint file'
            CHK%FDESC(2) = 'From CMAQ model dyn alloc version CTM'
            CHK%NVARS = N_GC_SPC + 1 + N_AE_SPC + N_NR_SPC + N_TR_SPC
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            INDX = 0
            DO S = INDX+1, INDX+N_GC_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = GC_SPC(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_GC_SPC
            CHK%VTYPE(INDX+1) = M3REAL
            CHK%VNAME(INDX+1) = 'RHOJ'
            CHK%UNITS(INDX+1) = 'm*kg/m**3'
            CHK%VDESC(INDX+1) = 'Variable '//TRIM(CHK%VNAME(INDX+1))
            INDX = INDX + 1
            DO S = INDX+1, INDX+N_AE_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = AE_SPC(S-INDX)
               IF ( CHK%VNAME(S)(1:3) == 'NUM' ) THEN
                  CHK%UNITS(S) = 'number/m**3'
               ELSE IF ( CHK%VNAME(S)(1:3) == 'SRF' ) THEN
                  CHK%UNITS(S) = 'm**2/m**3'
               ELSE
                  CHK%UNITS(S) = 'micrograms/m**3'
               END IF
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_AE_SPC
            DO S = INDX+1, INDX+N_NR_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = NR_SPC(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_NR_SPC
            DO S = INDX+1, INDX+N_TR_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = TR_SPC(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_TR_SPC

         CASE ( ADJ_VDIFF_CHK )
            ! No longer want to limit the vertical diffusion checkpointing to the first layer
            !    - need to account for stack emissions (6.2015.slc)
            !CHK%NLAYS = 1
            CHK%FDESC(1) = 'Vertical diffusion checkpoint file'
            CHK%FDESC(2) = 'From CMAQ model dyn alloc version CTM'
            CHK%NVARS = N_GC_SPC + 1 + N_AE_SPC + N_NR_SPC + N_TR_SPC
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            INDX = 0
            DO S = INDX+1, INDX+N_GC_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = GC_SPC(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_GC_SPC
            CHK%VTYPE(INDX+1) = M3REAL
            CHK%VNAME(INDX+1) = 'RHOJ'
            CHK%UNITS(INDX+1) = 'm*kg/m**3'
            CHK%VDESC(INDX+1) = 'Variable '//TRIM(CHK%VNAME(INDX+1))
            INDX = INDX + 1
            DO S = INDX+1, INDX+N_AE_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = AE_SPC(S-INDX)
               IF ( CHK%VNAME(S)(1:3) == 'NUM' ) THEN
                  CHK%UNITS(S) = 'number/m**3'
               ELSE IF ( CHK%VNAME(S)(1:3) == 'SRF' ) THEN
                  CHK%UNITS(S) = 'm**2/m**3'
               ELSE
                  CHK%UNITS(S) = 'micrograms/m**3'
               END IF
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_AE_SPC
            DO S = INDX+1, INDX+N_NR_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = NR_SPC(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_NR_SPC
            DO S = INDX+1, INDX+N_TR_SPC
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = TR_SPC(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_TR_SPC

         CASE ( ADJ_EMIST_CHK )
            ! No longer want to limit the vertical diffusion checkpointing to the first layer
            !    - need to account for stack emissions (6.2015.slc)
            !CHK%NLAYS = 1
            CHK%FDESC(1) = 'Emissions checkpoint file'
            CHK%FDESC(2) = 'From CMAQ model dyn alloc version CTM'
            CHK%NVARS = N_GC_DIFF + N_AE_DIFF + N_NR_DIFF + N_TR_DIFF
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            INDX = 0
            DO S = INDX+1, INDX+N_GC_DIFF
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = GC_DIFF(S-INDX)
               CHK%UNITS(S) = 'mol'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_GC_DIFF
!            CHK%VTYPE(INDX+1) = M3REAL
!            CHK%VNAME(INDX+1) = 'RHOJ'
!            CHK%UNITS(INDX+1) = 'm*kg/m**3'
!            CHK%VDESC(INDX+1) = 'Variable '//TRIM(CHK%VNAME(INDX+1))
!            INDX = INDX + 1
            DO S = INDX+1, INDX+N_AE_DIFF
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = AE_DIFF(S-INDX)
               IF ( CHK%VNAME(S)(1:3) == 'NUM' ) THEN
                  CHK%UNITS(S) = 'number/m**3'
               ELSE IF ( CHK%VNAME(S)(1:3) == 'SRF' ) THEN
                  CHK%UNITS(S) = 'm**2/m**3'
               ELSE
                  CHK%UNITS(S) = 'micrograms/m**3'
               END IF
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_AE_DIFF
            DO S = INDX+1, INDX+N_NR_DIFF
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = NR_DIFF(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_NR_DIFF
            DO S = INDX+1, INDX+N_TR_DIFF
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = TR_DIFF(S-INDX)
               CHK%UNITS(S) = 'ppmV'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_TR_DIFF

         CASE ( ADJ_HA_RHOJ_CHK )
            CHK%FDESC(1) = 'Horizontal advection RHOJ checkpoint file'
            CHK%FDESC(2) = 'From CMAQ model dyn alloc version CTM'
            CHK%NVARS = 1
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            CHK%VTYPE(1) = M3REAL
            CHK%VNAME(1) = 'RHOJ'
            CHK%UNITS(1) = 'KG/M**2'
            CHK%VDESC(1) = 'advected air density X total Jacobian'

         CASE ( ADJ_VA_RHOJ_CHK )
            CHK%FDESC(1) = 'Vertical advection RHOJ checkpoint file'
            CHK%FDESC(2) = 'From CMAQ model dyn alloc version CTM'
            CHK%NVARS = 1
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            CHK%VTYPE(1) = M3REAL
            CHK%VNAME(1) = 'RHOJ'
            CHK%UNITS(1) = 'KG/M**2'
            CHK%VDESC(1) = 'advected air density X total Jacobian'

         CASE ( ADJ_HADV_CHK )
            CHK%FDESC(1) = 'Horizontal advection checkpoint file'
            CHK%FDESC(2) = ' '
            CHK%NVARS = N_GC_ADV + N_AE_ADV + N_NR_ADV + N_TR_ADV + 1
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            INDX = 0
            CHK%VNAME(INDX+1:INDX+N_GC_ADV) = GC_ADV(1:N_GC_ADV)
            INDX = INDX + N_GC_ADV
            CHK%VNAME(INDX+1:INDX+N_AE_ADV) = AE_ADV(1:N_AE_ADV)
            INDX = INDX + N_AE_ADV
            CHK%VNAME(INDX+1:INDX+N_NR_ADV) = NR_ADV(1:N_NR_ADV)
            INDX = INDX + N_NR_ADV
            CHK%VNAME(INDX+1:INDX+N_TR_ADV) = TR_ADV(1:N_TR_ADV)
            INDX = INDX + N_TR_ADV
            CHK%VNAME(INDX + 1) = 'RHOJ'
            CHK%VTYPE(1:CHK%NVARS) = M3REAL
            CHK%UNITS(1:CHK%NVARS) = 'Coupled CGRID'
            CHK%VDESC(1:CHK%NVARS) = 'Horizontal advection checkpoint for ' //
     &                                CHK%VNAME(1:CHK%NVARS)

         CASE ( ADJ_VADV_CHK )
            CHK%FDESC(1) = 'Vertical advection checkpoint file'
            CHK%FDESC(2) = ' '
            CHK%NVARS = N_GC_ADV + N_AE_ADV + N_NR_ADV + N_TR_ADV + 1
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            INDX = 0
            CHK%VNAME(INDX+1:INDX+N_GC_ADV) = GC_ADV(1:N_GC_ADV)
            INDX = INDX + N_GC_ADV
            CHK%VNAME(INDX+1:INDX+N_AE_ADV) = AE_ADV(1:N_AE_ADV)
            INDX = INDX + N_AE_ADV
            CHK%VNAME(INDX+1:INDX+N_NR_ADV) = NR_ADV(1:N_NR_ADV)
            INDX = INDX + N_NR_ADV
            CHK%VNAME(INDX+1:INDX+N_TR_ADV) = TR_ADV(1:N_TR_ADV)
            INDX = INDX + N_TR_ADV
            CHK%VNAME(INDX + 1) = 'RHOJ'
            CHK%VTYPE(1:CHK%NVARS) = M3REAL
            CHK%UNITS(1:CHK%NVARS) = 'Coupled CGRID'
            CHK%VDESC(1:CHK%NVARS) = 'Vertical advection checkpoint for ' //
     &                                CHK%VNAME(1:CHK%NVARS)

         ! Added to checkpoint emissions (mdt|slc.1.2014)
         ! Suboptimal method for incorporating emissions sensitivities
         CASE ( ADJ_EMIS_CHK )
            CHK%FDESC(1) = 'Emissions checkpoint file'
            CHK%FDESC(2) = 'From CMAQ model dyn alloc version CTM'
            CHK%NVARS = N_GC_EMIS + N_AE_EMIS + N_NR_EMIS + N_TR_EMIS
            ALLOCATE(CHK%VTYPE(CHK%NVARS))
            ALLOCATE(CHK%VNAME(CHK%NVARS))
            ALLOCATE(CHK%UNITS(CHK%NVARS))
            ALLOCATE(CHK%VDESC(CHK%NVARS))
            ALLOCATE(CHK%VINDX(CHK%NVARS))
            INDX = 0
            DO S = INDX+1, INDX+N_GC_EMIS
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = GC_EMIS(S-INDX)
               CHK%UNITS(S) = 'mol/s'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO
            INDX = INDX + N_GC_EMIS
            DO S = INDX+1, INDX+N_AE_EMIS
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = AE_EMIS(S-INDX)
               IF ( CHK%VNAME(S)(1:3) == 'NUM' ) THEN
                  CHK%UNITS(S) = 'number/m**3'
               ELSE IF ( CHK%VNAME(S)(1:3) == 'SRF' ) THEN
                  CHK%UNITS(S) = 'm**2/m**3'
               ELSE
                  CHK%UNITS(S) = 'g/s'
               END IF
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(CHK%NVARS))
            END DO 
            INDX = INDX + N_AE_EMIS
            DO S = INDX+1, INDX+N_NR_EMIS
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = NR_EMIS(S-INDX)
               CHK%UNITS(S) = 'mol/s'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO 
            INDX = INDX + N_NR_EMIS
            DO S = INDX+1, INDX+N_TR_EMIS
               CHK%VTYPE(S) = M3REAL
               CHK%VNAME(S) = TR_EMIS(S-INDX)
               CHK%UNITS(S) = 'mol/s'
               CHK%VDESC(S) = 'Variable '//TRIM(CHK%VNAME(S))
            END DO 
            INDX = INDX + N_TR_EMIS

         CASE DEFAULT
            XMSG =  '"' // TRIM(CHK%FNAME) // '"' //
     &           ' is not a valid file selection'
            CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT3)

      END SELECT

      IF ( CREATENEW ) THEN
         IF ( CONST_SYNC_STEP .AND. .NOT. PRESENT(TSTEP) ) THEN
            XMSG = 'TSTEP must be a calling argument for CHK_FILE_OPEN'
            CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT1)
         END IF

         CALL CHK_FILE_CREATE_IMPL(CHK_ID, TSTEP)
      ELSE
         CALL CHK_FILE_OPEN_IMPL(CHK_ID)
      ENDIF

      RETURN
      END SUBROUTINE CHK_FILE_OPEN


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_CLOSE
!-----------------------------------------------------------------------
! Function:
!     Subroutine closes checkpoint files

! Subroutines and functions called:
!     CHK_FILE_CLOSE_IMPL

! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      USE CHK_FILES_IMPL
      IMPLICIT NONE

! Local variables
      INTEGER                 :: CHK_ID       ! checkpoint file identification

!-----------------------------------------------------------------------


!     ----Close checkpoint files----
      DO CHK_ID = 1,N_CHKS
         IF ( CHKS(CHK_ID)%ISUSED ) THEN
            CALL CHK_FILE_CLOSE_IMPL (CHK_ID)
         END IF
      END DO

      RETURN

      END SUBROUTINE CHK_FILE_CLOSE


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_REFRESH
!-----------------------------------------------------------------------
! Function:
!     Subroutine refreshes checkpoint files in all processors
!     Fake subroutine in some implementations

! Subroutines and functions called:
!     CHK_FILE_REFRESH_IMPL

! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      USE CHK_FILES_IMPL
      IMPLICIT NONE

! Local variables
      INTEGER                 :: CHK_ID       ! checkpoint file identification

!      ----Refresh checkpoint files----
      DO CHK_ID = 1,N_CHKS
        IF ( CHKS(CHK_ID)%ISUSED ) THEN
           CALL CHK_FILE_REFRESH_IMPL (CHK_ID)
        END IF
      END DO


      END SUBROUTINE CHK_FILE_REFRESH


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_READ ( CHK_ID, JDATE, JTIME, CHK_GRID )
!-----------------------------------------------------------------------
! Function:
!     Subroutine to perform reading from checkpoint file CHK_GRID
!
! Input:
!     CHK_ID, JDATE, JTIME
!
! Output:
!     CHK_GRID
!
! Subroutines and functions called:
!     CHK_FILE_READ_IMPL, M3EXIT
!
! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      USE CHK_FILES_IMPL
      IMPLICIT NONE

! Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

! Input parameters

      INTEGER                 :: CHK_ID       ! chk file id
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS
      REAL                    :: CHK_GRID( :,:,:,: )

!  Local variables:
      CHARACTER(len=NAMLEN_ADJ) :: PNAME = 'CHK_FILE_READ'
      CHARACTER( 250 )         :: XMSG = ' '

!-----------------------------------------------------------------------

      IF ( .NOT. CHKS(CHK_ID)%ISUSED ) THEN
         XMSG = 'The chkfile '//TRIM(CHK_NAMES(CHK_ID))//' is not opened and can not be read.'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF
      CALL CHK_FILE_READ_IMPL ( CHK_ID, JDATE, JTIME, CHK_GRID )

      RETURN

      END SUBROUTINE CHK_FILE_READ



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_WRITE ( CHK_ID, JDATE, JTIME, CHK_GRID )
!-----------------------------------------------------------------------
! Function:
!     Subroutine writes array CHK_GRID to the checkpoint file CHK_ID
!
! Input:
!     CHK_ID, JDATE, JTIME, CHK_GRID
!
! Subroutines and functions called:
!     CHK_FILE_WRITE_IMPL, M3EXIT

! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      USE CHK_FILES_IMPL
      IMPLICIT NONE

! Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

! input parameters:

      INTEGER                 :: CHK_ID       ! chk file id
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS
      REAL                    :: CHK_GRID( :,:,:,: )

!     Local variables:
      CHARACTER(len=NAMLEN_ADJ) :: PNAME = 'CHK_FILE_WRITE'
      CHARACTER( 250 )         :: XMSG = ' '

!-----------------------------------------------------------------------

      IF ( .NOT. CHKS(CHK_ID)%ISUSED ) THEN
         XMSG = 'The file '//TRIM(CHK_NAMES(CHK_ID))//' is not opened and can not be written.'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF
      CALL CHK_FILE_WRITE_IMPL ( CHK_ID, JDATE, JTIME, CHK_GRID )

      RETURN

      END SUBROUTINE CHK_FILE_WRITE


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


      END MODULE CHK_FILES
