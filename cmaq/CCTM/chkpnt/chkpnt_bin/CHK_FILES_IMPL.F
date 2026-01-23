
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
!   Implementation of the module for manipulation of checkpoint files.
!   Multiple files strategy with direct access fortran binary files.
!   This code should NEVER be used directly!
!   Subroutines from the module CHK_FILES should be called instead.
!   The implementation DOES NOT require fixed timestep.
!   Number of timesteps written is saved in an extra file ADJ_TIMESTEP_CHK.
!   If the unit is compiled with ifort, the "-assume byterecl" option MUST be used.
!   Gfortran and PGF90 require no aditional options
!
! Revision History:
!   June 2011 by J.Resler, ICS: Initial version
!   Aug 2013 by Peter Percell at UH-IMAQS:
!   -- reduced vdiff checkpoint file to just layer 1 of cgrid
!-----------------------------------------------------------------------

      IMPLICIT NONE
      
      INTEGER, PARAMETER                :: ADJ_TIMESTEP_CHK = 0        ! id number of timestep checkpoint file
      CHARACTER(16), PARAMETER          :: TNAME = 'ADJ_TIMESTEP_CHK'  ! name of timestep checkpoint file
      INTEGER, SAVE                      :: TDATE                       ! initial date of timestep checkpoint file
      INTEGER, SAVE                      :: TTIME                       ! initial time of timestep checkpoint file
      INTEGER, SAVE                      :: STEP                        ! current timestep
      INTEGER, SAVE                      :: PREVDATE                    ! date of previous checkpoint file operation
      INTEGER, SAVE                      :: PREVTIME                    ! time of previous checkpoint file operation
      INTEGER, SAVE                      :: STEPNUM = -1                ! logical file unit number of timestep file
      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_CREATE_IMPL ( CHK_ID )
!-----------------------------------------------------------------------
!
!  FUNCTION:  Create file CHK_ID using info stored in the CHK_FILE
!             description.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!   July 2011 by J.Resler, ICS: Initial code
!-----------------------------------------------------------------------

        USE GRID_CONF             ! horizontal & vertical domain specifications
        USE CHK_FILES_DEFS
        IMPLICIT NONE

!...........   INCLUDES:
        INCLUDE SUBST_IOPARMS
        INCLUDE SUBST_STATE
        INCLUDE SUBST_IOFDESC


!  ARGUMENTS and their descriptions:

       INTEGER                    :: CHK_ID !  identification of the checkpoint file


!  LOCAL VARIABLES and their descriptions:
!
        CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_CREATE'
        INTEGER                  :: IERR            !  netCDF error status return

        type(CHK_FILE),pointer   :: CHK
        CHARACTER(255)           :: FULLNAME        ! checkpoint file name incl. path
        INTEGER                  :: RECLEN          ! length of the record
        REAL                      :: REALVAR
        LOGICAL                  :: EX

        LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!.............................................................................
!   begin body of subroutine  CHK_CREATE_FILE

!     Checkpoint structure and file name and path
      CHK => CHKS(CHK_ID)
      print*, 'STEPNUM=',STEPNUM, MYPE
!.......... Check if timestep chkpoint file exists
      IF ( STEPNUM.LT.0 ) THEN
!.......... Create timestep chkpoint file and initialize it
        STEPNUM = CHK_FILE_UNITID()
        TDATE = CHK%SDATE
        TTIME = CHK%STIME
      print*, 'TDATE,TTIME=',TDATE,TTIME
        CALL CHK_FILE_NAME_IMPL (ADJ_TIMESTEP_CHK, FULLNAME)
        RECLEN = (sizeof(STEP) + sizeof(PREVDATE) + sizeof(PREVTIME))
        OPEN(UNIT=STEPNUM,FILE=FULLNAME,FORM='UNFORMATTED',ACCESS='DIRECT',
     &        RECL=RECLEN,STATUS='UNKNOWN',ACTION='READWRITE',IOSTAT=IERR)
        IF ( IERR .NE. 0 ) THEN
           CALL M3EXIT(PNAME, TDATE, TTIME, 'Error creating timestep chkfile', IERR)
        END IF
!........... set the step and time
        STEP = 1
        PREVDATE = CHK%SDATE
        PREVTIME = CHK%STIME
        WRITE(UNIT=STEPNUM,REC=1,IOSTAT=IERR) STEP, PREVDATE, PREVTIME
        IF ( IERR .NE. 0 ) THEN
           CALL M3EXIT(PNAME, PREVDATE, PREVTIME, 'Error writing timestep chkfile', IERR)
        END IF
      ENDIF

!     Checkpoint file name and path
      CALL CHK_FILE_NAME_IMPL (CHK_ID, FULLNAME)

!..........  Get parameters of the new file
      CHK%FNUM = CHK_FILE_UNITID()
      RECLEN = (sizeof(CHK%SDATE)+sizeof(CHK%STIME) + sizeof(REALVAR)*CHK%NVARS*NCOLS*NROWS*CHK%NLAYS) ! suppose full grid checkpointing
      print*, 'RECLEN =', RECLEN
      print*, 'sizeof(REALVAR) =', sizeof(REALVAR)
      print*, 'CHK%NVARS =', CHK%NVARS
      print*, 'NCOLS*NROWS*CHK%NLAYS =', NCOLS*NROWS*CHK%NLAYS

!...........  Create a new file
      OPEN(UNIT=CHK%FNUM,FILE=FULLNAME,FORM='UNFORMATTED',ACCESS='DIRECT',
     &      RECL=RECLEN,STATUS='UNKNOWN',ACTION='READWRITE',IOSTAT=IERR)
      IF ( IERR .NE. 0 ) THEN
          CALL M3EXIT(PNAME, PREVDATE, PREVTIME, 'Error creating CHK file', IERR)
      END IF


      END SUBROUTINE CHK_FILE_CREATE_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_OPEN_IMPL (CHK_ID)
!-----------------------------------------------------------------------
!
!  FUNCTION:  Opens binary checkpoint file CHK_ID
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!
!  REVISION  HISTORY:
!   June 2011 by J.Resler, ICS: Initial version
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      USE SUBST_MODULES         ! stenex
      IMPLICIT NONE

!...........   INCLUDES:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC

!  ARGUMENTS and their descriptions:

      INTEGER                  :: CHK_ID !  identification of the checkpoint file

!  LOCAL VARIABLES and their descriptions:

      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_OPEN'
      type(CHK_FILE),pointer   :: CHK
      INTEGER                  :: IERR
      CHARACTER(255)           :: FULLNAME      ! checkpoint file name incl. path
      INTEGER                  :: RECLEN          ! length of the record
      REAL                      :: REALVAR

      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!.............................................................................
!   begin body of subroutine  CHK_FILE_OPEN_IMPL

!.......... Check if timestep chkpoint file exists
      IF ( STEPNUM.LT.0 ) THEN
!------ open timestep chkpoint file ----
        STEPNUM = CHK_FILE_UNITID()
        TDATE = CHK%SDATE
        TTIME = CHK%STIME
        CALL CHK_FILE_NAME_IMPL (ADJ_TIMESTEP_CHK, FULLNAME)
        RECLEN = (sizeof(STEP) + sizeof(PREVDATE) + sizeof(PREVTIME))
        OPEN(STEPNUM,FILE=FULLNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN,STATUS='OLD',ACTION='READ',IOSTAT=IERR)
        IF ( IERR .NE. 0 ) THEN
          CALL M3EXIT(PNAME, TDATE, TTIME, 'Error opening timestep chkfile', IERR)
        END IF
!------ read timestep chkpoint file ----
        READ(STEPNUM,REC=1,IOSTAT=IERR) STEP, PREVDATE, PREVTIME
        IF ( IERR .NE. 0 ) THEN
          CALL M3EXIT(PNAME, PREVDATE, PREVTIME, 'Error reading timestep chkfile', IERR)
        END IF
      END IF

!     Checkpoint structure and file name and path
      CHK => CHKS(CHK_ID)
      CALL CHK_FILE_NAME_IMPL (CHK_ID, FULLNAME)

!..........  Get parameters of the new file
      CHK%FNUM = CHK_FILE_UNITID()
      RECLEN = (sizeof(CHK%SDATE)+sizeof(CHK%STIME) + sizeof(REALVAR)*CHK%NVARS*NCOLS*NROWS*CHK%NLAYS) ! suppose full grid checkpointing

!...........  Open an old file
      OPEN(UNIT=CHK%FNUM,FILE=FULLNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN,STATUS='OLD',ACTION='READ',IOSTAT=IERR)
      IF ( IERR .NE. 0 ) THEN
          CALL M3EXIT(PNAME, PREVDATE, PREVTIME, 'Error opening CHK file', IERR)
      END IF

      RETURN

      END SUBROUTINE CHK_FILE_OPEN_IMPL


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_CLOSE_IMPL (CHK_ID)
!-----------------------------------------------------------------------
! Function:
!     Function closess checkpoint files

! Subroutines and functions called:
!     CLOSE3, M3WARN

!  REVISION  HISTORY:
!   July 2010 by J.Resler, ICS: Initial version
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      IMPLICIT NONE

!      INCLUDE SUBST_NETCDF

!     Parameters
      INTEGER                  :: CHK_ID

!     Local variables:
      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_CLOSE'
      CHARACTER( 96 )          :: XMSG = ' '
      type(CHK_FILE),pointer   :: CHK
      LOGICAL                  :: OPENED
      INTEGER                  :: IERR


      INTEGER, SAVE           :: LOGDEV
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      CHK => CHKS(CHK_ID)

      CLOSE(CHK%FNUM, IOSTAT=IERR)
      IF ( IERR .NE. 0 ) THEN
         XMSG = 'Could not close ' // TRIM(CHK%FNAME)
         CALL M3WARN( PNAME, CHK%SDATE, CHK%STIME, XMSG )
      END IF

      INQUIRE ( UNIT=STEPNUM, OPENED=OPENED, IOSTAT=IERR )
      IF ( IERR .NE. 0 ) THEN
         CALL M3EXIT(PNAME, CHK%SDATE, CHK%STIME, 'Error testing the timestep chk file.', IERR)
         RETURN
      END IF

      IF ( OPENED ) THEN
         CLOSE(STEPNUM, IOSTAT=IERR)
         IF ( IERR .NE. 0 ) THEN
            CALL M3WARN( PNAME, CHK%SDATE, CHK%STIME, 'Could not close the timestep chk file.' )
         END IF
      END IF

      RETURN

      END SUBROUTINE CHK_FILE_CLOSE_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_REFRESH_IMPL (CHK_ID)
!-----------------------------------------------------------------------
! Function:
!     Subroutine refreshes checkpoint files in all processors
!     Fake subroutine in this implementation
!
! Revision History:
!   July 2010 by J.Resler, ICS: Initial version
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

!     Parameters
      INTEGER                  :: CHK_ID

      RETURN

      END SUBROUTINE CHK_FILE_REFRESH_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_READ_IMPL ( CHK_ID, JDATE, JTIME, CHK_GRID )
!-----------------------------------------------------------------------
! Function:
!     Function implements reading from checkpoint file CONC_CHK

! OUTPUT:
!     CHK_GRID

! Revision History:
!     June 2011 by J.Resler, ICS: Initial version
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      IMPLICIT NONE

! Include Files:
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

! Parameters
      INTEGER                 :: CHK_ID       ! chk file id
      REAL                    :: CHK_GRID( :,:,:,: )
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS

! Local variables:
      type(CHK_FILE),pointer  :: CHK          ! chk file description struct
      INTEGER                 :: IERR

      INTEGER                 :: CDATE        ! checkpoint record date, coded YYYYDDD
      INTEGER                 :: CTIME        ! checkpoint record time, coded HHMMSS

      CHARACTER(NAMLEN_ADJ)   :: PNAME = 'CHK_FILE_READ'
      CHARACTER( 250 )        :: XMSG = ' '

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      IF (FIRSTIME) THEN
         LOGDEV = INIT3()
         FIRSTIME = .FALSE.
      ENDIF

      CHK => CHKS(CHK_ID)

      IF (JDATE.NE.PREVDATE .OR. JTIME.NE.PREVTIME ) THEN
!-------- new timestep arises ------
         IF ( JDATE.GT.PREVDATE .OR. (JDATE.EQ.PREVDATE .AND. JTIME.GT.PREVTIME) ) THEN
            XMSG = 'Inconsistent time step for file:  '//CHK%FNAME
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, -1)
            RETURN
         END IF

         STEP = STEP - 1
         IF ( STEP .LT. 0 ) THEN
            XMSG = 'Time step error for file:  '//CHK%FNAME
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, -1)
            RETURN
         END IF
      END IF

      write(logdev,*) 'Reading chk_file: ', CHK%FNUM, CHK%FNAME, JDATE, JTIME, STEP

      READ(CHK%FNUM, REC=STEP, IOSTAT=IERR) CDATE, CTIME, CHK_GRID
      IF ( IERR .NE. 0 ) THEN
          WRITE(XMSG,*) 'Error reading from CHK file '//TRIM(CHK%FNAME)//', STEP=', STEP
          CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, IERR)
      END IF

!------- Check consistency of the read record with requisted timestep
      IF ( CDATE .NE. JDATE .OR. CTIME .NE. JTIME ) THEN
          WRITE(XMSG,*) 'Consistency error in file '//TRIM(CHK%FNAME)//', STEP=', STEP, CDATE, CTIME
          CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, IERR)
      END IF

!------- Save the last read timestep
      PREVDATE = JDATE
      PREVTIME = JTIME

      RETURN

      END SUBROUTINE CHK_FILE_READ_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_WRITE_IMPL ( CHK_ID, JDATE, JTIME, CHK_GRID )
!-----------------------------------------------------------------------
! Function:
!     Function to perform writing to the checkpoint file CONC_CHK

! INPUT:
!     CHK_GRID

! Revision History:
!     June 2011 by J.Resler, ICS: Initial version
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      IMPLICIT NONE

! Include Files:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

      INTEGER                 :: CHK_ID       ! chk file id
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS
      REAL                    :: CHK_GRID( :,:,:,: )

! Local variables:
      type(CHK_FILE),pointer  :: CHK          ! chk file description struct
      INTEGER                 :: IERR
      integer                 :: reclen

      CHARACTER(NAMLEN_ADJ)   :: PNAME = 'CHK_FILE_WRITE'
      CHARACTER( 250 )        :: XMSG = ' '

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------


      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()
      ENDIF

      CHK => CHKS(CHK_ID)

      IF (JDATE.NE.PREVDATE .OR. JTIME.NE.PREVTIME ) THEN
!-------- new timestep arises ------
         IF ( JDATE.LT.PREVDATE .OR. (JDATE.EQ.PREVDATE .AND. JTIME.LT.PREVTIME) ) THEN
            XMSG = 'Inconsistent time step for file:  '//CHK%FNAME
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, -1)
            RETURN
         END IF
         PREVDATE = JDATE
         PREVTIME = JTIME
         STEP = STEP + 1
         WRITE(STEPNUM,REC=1,IOSTAT=IERR) STEP, PREVDATE, PREVTIME
         IF ( IERR .NE. 0 ) THEN
            CALL M3EXIT(PNAME, PREVDATE, PREVTIME, 'Error writing timestep chkfile', IERR)
         END IF
      END IF

      inquire(CHK%FNUM, RECL=reclen)
      write(logdev,*) 'Writing chk_file: ', CHK%FNUM, CHK%FNAME, JDATE, JTIME, STEP, reclen, size(CHK_GRID)
!-------- perform the write ----------
      WRITE(CHK%FNUM, REC=STEP, IOSTAT=IERR) JDATE, JTIME, CHK_GRID
      IF ( IERR .NE. 0 ) THEN
          WRITE(XMSG,*) 'Error writing CHK file '//TRIM(CHK%FNAME)//', STEP=', STEP
          CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, IERR)
      END IF

      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &         'Timestep written to', CHK%FNAME,
     &         'for date and time', JDATE, JTIME

      RETURN

      END SUBROUTINE CHK_FILE_WRITE_IMPL


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_NAME_IMPL (CHK_ID, FULLNAME)
!-----------------------------------------------------------------------
! Function:
!     Function returns the full name including path of the checkpoint file

! Subroutines and functions called:
!     M3WARN

!  REVISION  HISTORY:
!   July 2010 by J.Resler, ICS: Initial version
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      IMPLICIT NONE

! Include Files:
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

!     Parameters
      INTEGER, INTENT(IN)        :: CHK_ID       ! chk file id
      CHARACTER(*), INTENT(OUT)  :: FULLNAME


!     Local variables:
      type(CHK_FILE),pointer   :: CHK          ! chk file description struct
      CHARACTER(16)            :: FNAME
      INTEGER                  :: FDATE
      INTEGER                  :: FTIME
      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_NAME'
      CHARACTER( 96 )          :: XMSG = ' '
      INTEGER                  :: IERR

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()
      ENDIF

      IF (CHK_ID .EQ. ADJ_TIMESTEP_CHK) THEN
         FNAME =  TNAME
         FDATE = TDATE
         FTIME = TTIME
      ELSE
         CHK => CHKS(CHK_ID)
         FNAME = CHK_NAMES(CHK_ID)
         FDATE = CHK%SDATE
         FTIME = CHK%STIME
      END IF

      CALL CHK_FILE_PATH ( FULLNAME )

!   Create the real file name
      IF ( TRIM(FULLNAME) .EQ. '' ) THEN
         WRITE(LOGDEV,92001) TRIM(FNAME),FDATE, MYPE
         WRITE(FULLNAME,92001) TRIM(FNAME),FDATE, MYPE
92001    FORMAT (  A , "_", I7, "_", I3.3 )
      ELSE
         WRITE(LOGDEV,92002) TRIM(FULLNAME),TRIM(FNAME),FDATE, MYPE
         WRITE(FULLNAME,92002) TRIM(FULLNAME),TRIM(FNAME),FDATE, MYPE
92002    FORMAT (  A , "/", A , "_", I7, "_", I3.3 )
      END IF


      RETURN

      END SUBROUTINE CHK_FILE_NAME_IMPL


!-----------------------------------------------------------------------
      INTEGER FUNCTION CHK_FILE_UNITID ()
!-----------------------------------------------------------------------
! Function:
!     Function returns next unused number of logical unit

! Subroutines and functions called:
!

!  REVISION  HISTORY:
!   June 2011 by J.Resler, ICS: Initial version
!-----------------------------------------------------------------------

      USE CHK_FILES_DEFS
      IMPLICIT NONE

!     Local variables:

      INTEGER, SAVE             :: LASTUNITUSED = 50
      INTEGER                   :: MAXLUN = 1000
      LOGICAL                   :: OPENED
      INTEGER                   :: IERR
      INTEGER                   :: FNUM

      CHARACTER(NAMLEN_ADJ)     :: PNAME = 'CHK_FILE_UNITID'
      CHARACTER( 350 )          :: XMSG = ' '

!-----------------------------------------------------------------------

      CHK_FILE_UNITID = -1
!     Test next unit number and skip the occupied unit numbers.
      DO FNUM = LASTUNITUSED+1, MAXLUN
          INQUIRE ( UNIT=FNUM, OPENED=OPENED, IOSTAT=IERR )
          IF ( IERR .NE. 0 ) THEN
            CALL M3EXIT(PNAME, TDATE, TTIME, XMSG, IERR)
            RETURN
         END IF
 
         IF ( .NOT. OPENED ) THEN
            LASTUNITUSED = FNUM
            CHK_FILE_UNITID = FNUM
            RETURN
         END IF
      END DO

      XMSG = 'Max. number of opened files exceeded.'
      CALL M3EXIT(PNAME, PREVDATE, PREVTIME, XMSG, -1)
      RETURN

      END FUNCTION CHK_FILE_UNITID



      END MODULE CHK_FILES_IMPL

