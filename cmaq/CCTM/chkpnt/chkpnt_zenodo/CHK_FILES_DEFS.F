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
      MODULE CHK_FILES_DEFS

!-----------------------------------------------------------------------
! Function:
!   Manipulation of checkpoint files
!   This module contains declarations of the common types and variables
!   needed in both abstract module and implementation module
!
! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version from Peter Percell
!   Oct  2011 by Peter Percell at UH-IMAQS: removed unused checkpoint files
!   May  2013 by Peter Percell at UH-IMAQS:
!   - changed ADJ_RHOJ_CHK to ADJ_VA_RHOJ_CHK and added ADJ_HA_RHOJ_CHK
!   - removed the overly restrictive function CHK_FILE_TSTEP which forces
!     the user to pick the synchronization step by setting CTM_MAXSYNC and
!     CTM_MINSYNC to the same value in the run script
!   - added the logical flag CONST_SYNC_STEP to signal whether or not the
!     checkpointing implementation requires a constant synchronization step
!   Aug 2013 by Peter Percell at UH-IMAQS:
!   - reduced vdiff checkpoint file to just layer 1 of cgrid
!   Jan 2014 - S.Capps (US EPA) | M.Turner (UC Boulder)
!   - included emissions checkpointing file
!   Mar 2014 by Peter Percell, UH:
!     - replaced ADJ_XADV_CHK and ADJ_YADV_CHK with just ADJ_HADV_CHK
!     - renamed ADJ_ZADV_CHK to ADJ_VADV_CHK
!-----------------------------------------------------------------------

      IMPLICIT NONE


! Number of checkpoint files
      INTEGER, PARAMETER                   :: N_CHKS = 8

! Checkpoint files ids
      INTEGER, PARAMETER                   :: ADJ_CHEM_CHK = 1
      INTEGER, PARAMETER                   :: ADJ_AERO_CHK = 2
      INTEGER, PARAMETER                   :: ADJ_VDIFF_CHK = 3
      INTEGER, PARAMETER                   :: ADJ_HA_RHOJ_CHK = 4
      INTEGER, PARAMETER                   :: ADJ_VA_RHOJ_CHK = 5
      INTEGER, PARAMETER                   :: ADJ_HADV_CHK = 6
      INTEGER, PARAMETER                   :: ADJ_VADV_CHK = 7
      INTEGER, PARAMETER                   :: ADJ_CLD_CHK = 8

! Checkpoint file names
      CHARACTER(16), PARAMETER :: CHK_NAMES(N_CHKS) =
     & ( / 'ADJ_CHEM_CHK   ','ADJ_AERO_CHK   ','ADJ_VDIFF_CHK  ','ADJ_HA_RHOJ_CHK',
     &     'ADJ_VA_RHOJ_CHK','ADJ_HADV_CHK   ','ADJ_VADV_CHK   ','ADJ_CLD_CHK    ' / )


! We cannot include PARMS3.EXT here, define own parameters
      INTEGER, PARAMETER                     :: MXDLEN_ADJ =  80
      INTEGER, PARAMETER                     :: NAMLEN_ADJ =  16

! True if checkpointing implementation requires a constant synchronization time step
      LOGICAL, PARAMETER                     :: CONST_SYNC_STEP = .TRUE.

! Description of the checkpoint file
      TYPE CHK_FILE
        CHARACTER(NAMLEN_ADJ)                :: FNAME            ! file name pattern
        CHARACTER(MXDLEN_ADJ)                :: FDESC( 2 )       ! file description (MXDLEN_ADJ=80)
        INTEGER                              :: SDATE            ! starting date (YYYYDDD)
        INTEGER                              :: STIME            ! starting time (HHMMSS)
        INTEGER                              :: NLAYS            ! number of layers
        INTEGER                              :: NVARS            ! number of variables
        INTEGER,ALLOCATABLE                  :: VTYPE( : )       ! variable type:  M3(INT|REAL|DBLE)
        CHARACTER(NAMLEN_ADJ),ALLOCATABLE    :: VNAME( : )       ! variable names (length NAMLEN_ADJ=16)
        CHARACTER(NAMLEN_ADJ),ALLOCATABLE    :: UNITS( : )       ! variable units or 'none' (NAMLEN_ADJ=16)
        CHARACTER(MXDLEN_ADJ),ALLOCATABLE    :: VDESC( : )       ! variable descriptions (MXDLEN_ADJ=80)
! Variables initialized in CHK_FILE_CREATE_IMPL
        LOGICAL                              :: ISUSED = .FALSE. ! indicator of used file in particular configuration
        INTEGER                              :: TSTEP            ! time step for checkpoint files (HHMMSS)
                                                                 ! for implementations where the constant tstep is required
        INTEGER                              :: FNUM             ! nc file id
        INTEGER                              :: TINDX            ! netCDF var. ID for timestep-flag
        INTEGER,ALLOCATABLE                  :: VINDX( : )       ! netCDF var. IDs
      END TYPE CHK_FILE

! Array of the checkpoint files descriptions
      TYPE(CHK_FILE),TARGET                  :: CHKS(N_CHKS)

!-----------------------------------------------------------------------

      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_PATH ( CHK_PATH )
!-----------------------------------------------------------------------
! Function:
!     Subroutine reads env. variable CHK_PATH
!
! Output:
!     CHK_PATH
!
! Subroutines and functions called:
!     ENVSTR

! Revision History:
!   Oct  2010 by J.Resler, ICS: Initial version
!-----------------------------------------------------------------------

      IMPLICIT NONE

!  external functions


!  output parameters:
      CHARACTER(*)            :: CHK_PATH

!  local variables
      INTEGER                 :: STATUS

      CALL ENVSTR( 'CHK_PATH', 'The path where the checkpoint files will be stored', '', CHK_PATH, STATUS )

      IF ( STATUS .NE. 0 ) THEN
        CHK_PATH = ''
      ENDIF

      END SUBROUTINE CHK_FILE_PATH


      END MODULE CHK_FILES_DEFS

