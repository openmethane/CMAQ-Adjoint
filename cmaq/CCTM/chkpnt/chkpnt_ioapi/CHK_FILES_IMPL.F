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
      MODULE CHK_FILES_IMPL

!-----------------------------------------------------------------------
! Function:
!   Implementation of the module for manipulation of checkpoint files
!   Serial writes strategy with IOAPI3/PARIO
!
! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!                               inspired by Kumaresh Singh and Adrian Sandu
!                               original adjoint code and IOAPI3 code
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!   Nov  2011 by M. Turner, CU Boulder:
!     - Added output of filename when printing "Reading timestep"
!   Aug 2013 by Peter Percell at UH-IMAQS:
!   -- reduced vdiff checkpoint file to just layer 1 of cgrid
!-----------------------------------------------------------------------

      IMPLICIT NONE

      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_CREATE_IMPL (CHK_ID, TSTEP)
!-----------------------------------------------------------------------
!
!  FUNCTION:  Create netCDF file CHK_ID using info stored in the CHK_FILE
!             description.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     OPEN3, FLUSH3, ENVINT, SEC2TIME
!
!  REVISION  HISTORY:
!   Mar  2008 by J.Resler, ICS: Initial version based on CRTFIL3 from IOAPI-3
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!   Feb 2011 by Peter Percell at UH-IMAQS:
!   -- for this implementation, only need checkpoint files to be open for
!      update on processor 0; opening them for update on other processors
!      can causes problems
!   May 2013 by Peter Percell at UH-IMAQS:
!   -- added TSTEP as a calling argument to pass in a computed sync step
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      USE HGRD_DEFN
      USE VGRD_DEFN
      USE SUBST_MODULES         ! stenex
      IMPLICIT NONE

!...........   INCLUDES:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC

!  EXTERNAL FUNCTIONS and their descriptions

      INTEGER, EXTERNAL        :: ENVINT, SEC2TIME
      LOGICAL, EXTERNAL        :: FLUSH3     !  length after trimming trailing blanks.

!  ARGUMENTS and their descriptions:

      INTEGER             :: CHK_ID   !  identification of the checkpoint file
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                      !   TSTEP(1) = local output step
                                      !   TSTEP(2) = sciproc sync. step (chem)

!  LOCAL VARIABLES and their descriptions:

      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_CREATE_IMPL'
      type(CHK_FILE),pointer   :: CHK
      CHARACTER(255)           :: FNAME         ! checkpoint file name
      INTEGER                  :: V           !  loop counter for file variables
      INTEGER                  :: I,L
      LOGICAL                  :: IERR
      CHARACTER( 96 )          :: XMSG = ' '

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!.............................................................................
!   begin body of subroutine  CHK_FILE_CREATE_IMPL

      IF ( FIRSTIME ) THEN
         LOGDEV = INIT3 ()
         FIRSTIME = .FALSE.
      END IF

!     Checkpoint structure and file name
      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)

!     IOAPI3 implementation of the checkpointing
!     requires constant synchronization timestep
      IF ( .NOT. CONST_SYNC_STEP ) THEN
         XMSG = 'CMAQ adjoint requires fixed synchronization timestep.'
         CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT1 )
      END IF
      CHK%TSTEP = TSTEP(2)

!     Set output file characteristics based on COORD.EXT and open it
      FTYPE3D = GRDDED3
      SDATE3D = CHK%SDATE
      STIME3D = CHK%STIME
      TSTEP3D = CHK%TSTEP
      NLAYS3D = CHK%NLAYS
      DO L = 1, CHK%NLAYS + 1
         VGLVS3D( L ) = VGLVS_GD( L )
      END DO

      NCOLS3D = GL_NCOLS
      NROWS3D = GL_NROWS
      NTHIK3D = NTHIK
      GDTYP3D = GDTYP_GD
      P_ALP3D = P_ALP_GD
      P_BET3D = P_BET_GD
      P_GAM3D = P_GAM_GD
      XORIG3D = XORIG_GD
      YORIG3D = YORIG_GD
      XCENT3D = XCENT_GD
      YCENT3D = YCENT_GD
      XCELL3D = XCELL_GD
      YCELL3D = YCELL_GD
      VGTYP3D = VGTYP_GD
      VGTOP3D = VGTOP_GD
!     VGTPUN3D = VGTPUN_GD ! currently, not defined
      GDNAM3D = GDNAME_GD

!     ----concentration checkpoint files---
      FDESC3D( 1 ) = CHK%FDESC(1)
      FDESC3D( 2 ) = CHK%FDESC(2)
      DO L = 3, MXDESC3
         FDESC3D( L ) = ' '
      END DO

      NVARS3D = MAX( CHK%NVARS, 1 )
      DO V = 1 , CHK%NVARS
         VTYPE3D( V ) = CHK%VTYPE(V)
         VNAME3D( V ) = CHK%VNAME(V)
         UNITS3D( V ) = CHK%UNITS(V)
         VDESC3D( V ) = CHK%VDESC(V)
      END DO

      IF ( MYPE .EQ. 0 ) THEN

!     create checkpoint file
        IF ( .NOT. OPEN3( FNAME, FSNEW3, PNAME ) ) THEN
           XMSG = 'Could not open ' // trim(FNAME) // ' file'
           CALL M3EXIT( PNAME, SDATE3D, STIME3D, XMSG, XSTAT1 )
        END IF
        IERR = FLUSH3(FNAME)
      END IF

c$$$      CALL SUBST_BARRIER
c$$$
c$$$!     Open checkfiles in all other processes
c$$$      IF ( MYPE .GT. 0 ) THEN
c$$$          I = 0
c$$$          IERR = OPEN3( FNAME, FSRDWR3, PNAME )
c$$$!     refreshing sometime takes a little while on NFS disks
c$$$          do while (.not. IERR)
c$$$            call sleep(1)
c$$$            I = I+1
c$$$            write(logdev,*) trim(FNAME), ': ATTEMPT', I
c$$$            IERR = OPEN3( FNAME, FSRDWR3, PNAME )
c$$$          enddo
c$$$      END IF

      RETURN

      END SUBROUTINE CHK_FILE_CREATE_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_OPEN_IMPL (CHK_ID)
!-----------------------------------------------------------------------
!
!  FUNCTION:  Opens netCDF file CHK_ID (IOAPI3/PARIO)
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     OPEN3, FLUSH3, ENVINT, SEC2TIME
!
!  REVISION  HISTORY:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      USE HGRD_DEFN
      USE VGRD_DEFN
      USE SUBST_MODULES         ! stenex
      IMPLICIT NONE

!...........   INCLUDES:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC

!  EXTERNAL FUNCTIONS and their descriptions

      INTEGER, EXTERNAL        :: NAME2FID   !  fname~~> fid lookup

!  ARGUMENTS and their descriptions:

      INTEGER                  :: CHK_ID !  identification of the checkpoint file

!  LOCAL VARIABLES and their descriptions:

      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_OPEN_IMPL'
      type(CHK_FILE),pointer   :: CHK
      CHARACTER(255)           :: FNAME         ! checkpoint file name
      INTEGER                  :: FID           ! subscript for STATE3 arrays
      CHARACTER( 96 )          :: XMSG = ' '

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!.............................................................................
!   begin body of subroutine  CHK_FILE_OPEN_IMPL

!     Checkpoint structure and file name
      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)

!     open checkpoint file
      IF ( .NOT. OPEN3( FNAME, FSREAD3, PNAME ) ) THEN
         XMSG = 'Could not open ' // trim(FNAME) // ' file'
         CALL M3EXIT( PNAME, CHK%SDATE, CHK%STIME, XMSG, XSTAT1 )
      END IF

!  get the characteristics of the file
      IF ( .NOT. DESC3( FNAME ) ) THEN
         XMSG = 'Could not get ' // TRIM(FNAME) // ' file description'
         CALL M3EXIT( PNAME, CHK%SDATE, CHK%STIME, XMSG, XSTAT1 )
      END IF

      CHK%TSTEP = TSTEP3D
      CHK%SDATE = SDATE3D
      CHK%STIME = STIME3D

      RETURN

      END SUBROUTINE CHK_FILE_OPEN_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_CLOSE_IMPL (CHK_ID)
!-----------------------------------------------------------------------
! Function:
!     Subroutine closes checkpoint files
!
! Subroutines and functions called:
!     CLOSE3, M3WARN
!
! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      IMPLICIT NONE

!  INCLUDES:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

!  EXTERNAL FUNCTIONS
      LOGICAL, EXTERNAL       :: CLOSE3

!  ARGUMENTS and their descriptions:

      INTEGER                 :: CHK_ID

!  Local variables:
      CHARACTER( 16 )         :: PNAME = 'CHK_FILE_CLOSE'
      CHARACTER( 96 )         :: XMSG = ' '
      type(CHK_FILE),pointer  :: CHK
      CHARACTER(NAMLEN3)      :: FNAME

!-----------------------------------------------------------------------

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)
      IF ( .NOT. CLOSE3(FNAME ) ) THEN
         XMSG = 'Could not close ' // TRIM(FNAME)
         CALL M3WARN( PNAME, CHK%SDATE, CHK%STIME, XMSG )
      END IF

      RETURN
      END SUBROUTINE CHK_FILE_CLOSE_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_REFRESH_IMPL (CHK_ID)
!-----------------------------------------------------------------------
! Function:
!     Subroutine refreshes checkpoint files in all processes
!     The FLUSH3 doesn't work reliably on NFS, close-open used instead
!
! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      USE SUBST_MODULES         ! stenex
      IMPLICIT NONE

!  INCLUDES:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions

!  EXTERNAL FUNCTIONS
      LOGICAL, EXTERNAL       :: FLUSH3
      LOGICAL, EXTERNAL       :: OPEN3
      LOGICAL, EXTERNAL       :: CLOSE3

!  ARGUMENTS and their descriptions:

      INTEGER                 :: CHK_ID

      CHARACTER( 16 )         :: PNAME = 'CHK_FILE_REFRESH_IMPL'
      type(CHK_FILE),pointer  :: CHK
      CHARACTER( 96 )         :: XMSG = ' '
      CHARACTER(NAMLEN_ADJ)   :: FNAME

!-----------------------------------------------------------------------

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)

!     IF ( .NOT. FLUSH3(FNAME ) ) THEN
      IF ( .NOT. CLOSE3(FNAME ) ) THEN
          XMSG = 'Could not CLOSE ' // TRIM(FNAME)
          CALL M3WARN( PNAME, CHK%SDATE, CHK%STIME, XMSG )
      END IF

      CALL SUBST_BARRIER

!     IF ( .NOT. FLUSH3(FNAME ) ) THEN
      IF ( .NOT. OPEN3(FNAME,  FSRDWR3, PNAME ) ) THEN
          XMSG = 'Could not OPEN ' // TRIM(FNAME)
          CALL M3WARN( PNAME, CHK%SDATE, CHK%STIME, XMSG )
      END IF


      END SUBROUTINE CHK_FILE_REFRESH_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_READ_IMPL ( CHK_ID, JDATE, JTIME, CHK_GRID )
!-----------------------------------------------------------------------
! Function:
!     Function to perform reading from checkpoint file CONC_CHK
!
! OUTPUT:
!     CHK_GRID
!
! Revision History:
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: Rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!   April 2013 by Peter Percell at UH-IMAQS:
!   - As suggested by Jaroslav, XTRACT3 now writes directly to CHK_GRID(:,:,:,V)
!-----------------------------------------------------------------------

      USE HGRD_DEFN
      USE VGRD_DEFN
      USE CHK_FILES_DEFS
      IMPLICIT NONE

! Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
!      INCLUDE SUBST_FILES_ID    ! I/O definitions and declarations
!      INCLUDE SUBST_MPICH       ! MPI definitions and parameters

! EXTERNAL FUNCTIONS and their descriptions:
!      LOGICAL, EXTERNAL       :: DESC3
!      LOGICAL, EXTERNAL       :: XTRACT3
!      INTEGER, EXTERNAL       :: JSTEP3     !  compute time step record numbers

! Parameters
      INTEGER                 :: CHK_ID       ! chk file id
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS
      REAL                    :: CHK_GRID( :,:,:,: )

! Local variables:
      type(CHK_FILE),pointer  :: CHK          ! chk file description struct
      CHARACTER(NAMLEN3)      :: FNAME
      INTEGER                 :: R,C,L,V     ! loop counters

      CHARACTER(NAMLEN_ADJ)   :: PNAME = 'CHK_FILE_READ_IMPL'
      CHARACTER( 96 )         :: XMSG = ' '

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      IF (FIRSTIME) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
      ENDIF

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)

      write(logdev,*) 'Reading timestep: ', jdate, jtime, ' of ', FNAME

      IF ( .NOT. DESC3( FNAME ) ) THEN
         XMSG = 'Could not get ' // TRIM(FNAME) // ' file description'
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

      DO V = 1, CHK%NVARS
         IF ( .NOT. XTRACT3 ( FNAME, CHK%VNAME(V), 1, CHK%NLAYS,
     &                       ROWSX_PE(1,MYPE+1),ROWSX_PE(2,MYPE+1),
     &                       COLSX_PE(1,MYPE+1),COLSX_PE(2,MYPE+1),
     &                       JDATE, JTIME, CHK_GRID(:,:,:,V) ) ) THEN
            XMSG = 'Could not read  ' // CHK%VNAME( V )
     &           // ' from ' // TRIM(FNAME)
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
      END DO

      RETURN

      END SUBROUTINE CHK_FILE_READ_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_WRITE_IMPL ( CHK_ID, JDATE, JTIME, CHK_GRID )
!-----------------------------------------------------------------------
! Function:
!     Function performs writes to the checkpoint file CHK_ID
!
! INPUT:
!     CHK_ID, JDATE, JTIME
!
! OUTPUT:
!     CHK_GRID
!
!   Mar  2008 by J.Resler, ICS: Initial version
!   July 2010 by J.Resler, ICS: rewritten to allow different strategies
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      IMPLICIT NONE

! Include Files:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

!  EXTERNAL FUNCTIONS

#ifdef parallel
      LOGICAL, EXTERNAL       :: WRITE3
#endif

!  ARGUMENTS and their descriptions:

      INTEGER                 :: CHK_ID       ! chk file id
      REAL                    :: CHK_GRID( :,:,:,: )
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS


! Local variables:
      type(CHK_FILE),pointer  :: CHK          ! chk file description struct
      CHARACTER(NAMLEN3)      :: FNAME
      INTEGER                 :: V

      CHARACTER(NAMLEN_ADJ)   :: PNAME = 'CHK_FILE_WRITE_IMPL'
      CHARACTER( 96 )         :: XMSG = ' '

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()
      ENDIF

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)

!      IF ( .NOT. DESC3( FNAME ) ) THEN
!         XMSG = 'Could not get ' // TRIM(FNAME) // ' file description'
!         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!      END IF

      DO V = 1, CHK%NVARS
          IF ( .NOT. WRITE3( FNAME,  CHK%VNAME(V),
     &         JDATE, JTIME, CHK_GRID(:,:,:,V) ) ) THEN
              XMSG = 'Could not write '//CHK%VNAME(V)//' to '// TRIM(FNAME)
              CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF
      END DO

      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &     'Timestep written to', TRIM(FNAME),
     &     'for date and time', JDATE, JTIME


      RETURN

      END SUBROUTINE CHK_FILE_WRITE_IMPL



!-----------------------------------------------------------------------


      END MODULE CHK_FILES_IMPL




