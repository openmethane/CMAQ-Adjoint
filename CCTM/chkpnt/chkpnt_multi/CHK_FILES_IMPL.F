
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
!   Multiple files strategy with netcdf
!
! Revision History:
!   July 2010 by J.Resler, ICS: Initial version
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!   Aug 2013 by Peter Percell at UH-IMAQS:
!   -- reduced vdiff checkpoint file to just layer 1 of cgrid
!-----------------------------------------------------------------------

      IMPLICIT NONE
      
      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_CREATE_IMPL ( CHK_ID, TSTEP )
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
!   July 2010 by J.Resler, ICS: Initial code
!                               partialy inspired by CRTFIL3 from IOAPI3
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!   May 2013 by Peter Percell at UH-IMAQS:
!   -- added TSTEP as a calling argument to pass in a computed sync step
!-----------------------------------------------------------------------

        USE GRID_CONF             ! horizontal & vertical domain specifications
        USE CHK_FILES_DEFS
        USE NETCDF
        IMPLICIT NONE

!...........   INCLUDES:
        INCLUDE SUBST_IOPARMS
        INCLUDE SUBST_STATE
        INCLUDE SUBST_IOFDESC


!  ARGUMENTS and their descriptions:

        INTEGER             :: CHK_ID   !  identification of the checkpoint file
        INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                        !   TSTEP(1) = local output step
                                        !   TSTEP(2) = sciproc sync. step (chem)


!  LOCAL VARIABLES and their descriptions:
!
        CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_CREATE_IMPL'
        CHARACTER( 96 )          :: XMSG = ' '
        INTEGER                  :: ISCR
        INTEGER                  :: IERR            !  netCDF error status return
        INTEGER                  :: VAR             !  loop counter for file variables
        INTEGER                  :: TSIZE           !  time-step mode
        INTEGER                  :: DSIZE           !  D/T width
        INTEGER                  :: CDIM            !  col/level     dim #
        INTEGER                  :: RDIM            !  row/site      dim #
        INTEGER                  :: LDIM            !  layer         dim #
        INTEGER                  :: DDIM            !  timestamp date-and-time dimension #
        INTEGER                  :: TDIM            !  timestamp date-and-time dimension #
        INTEGER                  :: VDIM            !  variable-flag dim #
        INTEGER                  :: NDIMS           !  number of dims for NCVDEF()
        INTEGER                  :: DIMS( 4 )       !  array of dims for NCVDEF()
        INTEGER                  :: DELS( 4 )       !  array of dims for NCVPT()
        CHARACTER(NAMLEN_ADJ)    :: DSCBUF          !  scratch buffer for descriptions
        CHARACTER(MXDLEN_ADJ)    :: DSCBU2          !  scratch buffer for descriptions

        type(CHK_FILE),pointer   :: CHK
        CHARACTER(16)            :: FNAME           ! checkpoint file name
        CHARACTER(255)           :: FULLNAME        ! checkpoint file name incl. path
        INTEGER                  :: FNUM            ! netcdf file number
        INTEGER                  :: FMODE           ! netcdf creation mode

        INTEGER                   :: VLEN
        INTEGER,parameter       :: ONELEN=1
        INTEGER                   :: IVAL(1)
        REAL                      :: RVAL(1)
        DOUBLE PRECISION         :: DVAL(1)

        REAL*8                     :: XORIG         ! xorig of the local patch (REAL*8 in accordance with IOAPI-3)
        REAL*8                     :: YORIG         ! yorig of the local patch

        CHARACTER(MXVARS3*NAMLEN3):: CHARVARS
        CHARACTER(2*MXDLEN3)      :: CHARDESC
        CHARACTER(MXDLEN3*MXDESC3):: CHARSCNDSC
        INTEGER                  :: I

!.............................................................................
!   begin body of subroutine  CHK_CREATE_FILE

!     Checkpoint structure and file name and path
      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)
      CALL CHK_FILE_NAME_IMPL (CHK_ID, FULLNAME)

!   the present implementation of the checkpointing
!   with multiple netcdf files requires constant
!   synchronization timestep if we wish to produce
!   IOAPI3 compatible files
!   TODO: rewrite without this requirement
      IF ( .NOT. CONST_SYNC_STEP ) THEN
         XMSG = 'CMAQ adjoint requires fixed synchronization timestep.'
         CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT1 )
      END IF
      CHK%TSTEP = TSTEP(2)

!...........  Create the netCDF file, and put it into definition mode.

      IF ( CHK%NVARS .GT. MXVARS3 ) THEN
          CALL M3ABORT( FNAME, 0, IERR, 'Max NVARS for this build exceeded' )
          GO TO 999
      END IF


!...........  Create netCDF file and proceed...
      FMODE = NF90_CLOBBER
      IERR = NF90_CREATE(FULLNAME, FMODE, FNUM)
      IF ( IERR .NE. NF90_NOERR ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR, 'Error creating netCDF file' )
          CALL M3MSG2( FNAME )
          CALL NCCLOS( FNUM, IERR )
          GO TO 999
      END IF          !  ierr nonzero:  NCCRE() failed
      CHK%FNUM = FNUM
      write(logdev,*) 'po create'
!
!...........   Set attributes valid for all file types:
!.......   IOAPI_VERSION:  I/O API version and date
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'IOAPI_VERSION', TRIM( VERSN3 ))
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute IOAPI_VERSION' )
          GO TO 999
      END IF

!.......   EXEC_ID:  execution ID
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'EXEC_ID', EXECN3)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute EXEC_ID ' )
          GO TO 999
      END IF

!.......   FTYPE:  file type ID
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'FTYPE', GRDDED3)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute FTYPE' )
          GO TO 999
      END IF

!.......   CDATE:  creation date
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'CDATE', CURDATE)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute CDATE' )
          GO TO 999
      END IF

!.......   CTIME:  creation time
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'CTIME', CURTIME)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute CTIME' )
          GO TO 999
      END IF

!.......   WDATE:  date of last update
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'WDATE', CURDATE)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute WDATE' )
          GO TO 999
      END IF

!.......   WTIME:  time of last update
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL,'WTIME', CURTIME)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute WTIME' )
          GO TO 999
      END IF

!...............   Set up remaining data file attributes:
!.......   SDATE:STIME  starting date&time (normalized by NEXTIME(:,:,0) )
      CALL NEXTIME( CHK%SDATE, CHK%STIME, 0 )

      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'SDATE', CHK%SDATE)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute SDATE' )
          GO TO 999
      END IF

      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'STIME', CHK%STIME)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &    'Error creating netCDF file attribute STIME' )
          GO TO 999
      END IF

!.......   TSTEP:  time step
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'TSTEP', CHK%TSTEP)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &    'Error creating netCDF file attribute TSTEP' )
          GO TO 999
      END IF

!.......   NCOLS:  number of grid columns of the local processor patch
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'NCOLS', NCOLS)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &    'Error creating netCDF file attribute NCOLS' )
          GO TO 999
      END IF

!.......   NROWS:  number of grid rows of the local processor patch
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'NROWS', NROWS)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &    'Error creating netCDF file attribute NROWS' )
          GO TO 999
      END IF

!.......   NLAYS:  number of layers
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'NLAYS', CHK%NLAYS)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute NLAYS' )
          GO TO 999
      END IF

!.......   NVARS:  number of variables
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'NVARS', CHK%NVARS)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute NVARS' )
          GO TO 999
      END IF

!.......   GDTYP:  grid type ID (lat-lon, UTM, RADM, etc...)
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'GDTYP', GDTYP_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute GDTYP' )
          GO TO 999
      END IF

!.......   P_ALP:  first map-projection descriptive angle
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'P_ALP', P_ALP_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute P_ALP ' )
          GO TO 999
      END IF

!.......   P_BET:  second map-projection descriptive angle
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'P_BET', P_BET_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute P_BET ' )
          GO TO 999
      END IF

!.......   P_GAM:  third map-projection descriptive angle
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'P_GAM', P_GAM_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute P_GAM ' )
          GO TO 999
      END IF

!.......   XCENT:  lon of coordinate-system origin
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'XCENT', XCENT_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute XCENT ' )
          GO TO 999
      END IF

!.......   YCENT:  lat of coordinate-system origin
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'YCENT', YCENT_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute YCENT ' )
          GO TO 999
      END IF

!.......   XORIG:  X-coord of grid origin of the local processor patch
!.......  Find origins of local processor grid patch
      XORIG = XORIG_GD+(COLSX_PE(1,MYPE+1)-1)*XCELL_GD
      YORIG = YORIG_GD+(ROWSX_PE(1,MYPE+1)-1)*YCELL_GD

      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'XORIG', XORIG)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute XORIG ' )
          GO TO 999
      END IF

!.......   YORIG:  Y-coord of grid origin of the local processor patch
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'YORIG', YORIG)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute YORIG ' )
          GO TO 999
      END IF

!.......   XCELL:  cell width (X direction)
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'XCELL', XCELL_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute XCELL ' )
          GO TO 999
      END IF

!.......   YCELL:  cell width (Y direction)
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'YCELL', YCELL_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute YCELL ' )
          GO TO 999
      END IF

!.......   VGTYP:  vertical coordinate type (VGSIGP3, ... ) or IMISS3
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'VGTYP', VGTYP_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute VGTYP ' )
          GO TO 999
      END IF

!.......   VGTOP:
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'VGTOP', VGTOP_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute VGTOP ' )
          GO TO 999
      END IF

!.......   VGLVS( 1, ..., NLAYS3D+1 ):  vertical coordinate values
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'VGLVLS', VGLVS_GD )
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute VGLVLS ' )
          GO TO 999
      END IF

!.......   GDNAM:  grid name
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'GDNAM', GDNAME_GD)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute GDNAM ' )
          GO TO 999
      END IF

!.......   UPNAM:  name of last program to update the file
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'UPNAM', PNAME)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &        'Error creating netCDF file attribute UPNAM ' )
          GO TO 999
      END IF          !  ierr nonzero:  NCAPTC() failed

!.......   Variables-list for the file:
      IF( CHK%NVARS .GT. 0 ) THEN
          DO I = 1, CHK%NVARS
              CHARVARS(1+(I-1)*NAMLEN3:I*NAMLEN3) = CHK%VNAME(I)
          ENDDO
          IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'VAR-LIST', CHARVARS(1:NAMLEN3*CHK%NVARS))
          IF ( IERR .NE. 0 ) THEN
             CALL M3ABORT( FNAME, FNUM, IERR,
     &                'Error creating attribute VAR-LIST ' )
             GO TO 999
          END IF
      END IF

!.......   FILEDESC:  file description
      CHARDESC(1:MXDLEN3) = CHK%FDESC(1)
      CHARDESC(MXDLEN3+1:2*MXDLEN3) = CHK%FDESC(2)
      IERR = NF90_PUT_ATT(FNUM, NF90_GLOBAL, 'FILEDESC', CHARDESC)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating netCDF file attribute FILEDESC ' )
          GO TO 999
      END IF

!...............   Now define dimensions for the file:
!...........   First:  time-step and layer-number dimensions

      TSIZE = NF90_UNLIMITED
      DSIZE = 2               !  for ddim

      IERR = NF90_DEF_DIM(FNUM, 'TSTEP', TSIZE, TDIM)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating netCDF dimension TSTEP ' )
          GO TO 999
      END IF              !  ierr nonzero:  NCDDEF() failed

      IERR = NF90_DEF_DIM(FNUM, 'DATE-TIME', DSIZE, DDIM)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating netCDF dimension DATE-TIME' )
          GO TO 999
      END IF              !  ierr nonzero:  NCDDEF() failed

      IERR = NF90_DEF_DIM(FNUM, 'LAY', CHK%NLAYS, LDIM)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating netCDF dimension LAY' )
          GO TO 999
      END IF              !  ierr nonzero:  NCDDEF() failed

      IERR = NF90_DEF_DIM(FNUM, 'VAR', MAX( CHK%NVARS, 1 ), VDIM)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating netCDF dimension VAR ' )
          GO TO 999
      END IF              !  ierr nonzero:  NCDDEF() failed

      NDIMS = 3               !  for tindx(fid)
      DIMS( 1 ) = DDIM        !  date field or time field
      DIMS( 2 ) = VDIM        !  variable number
      DIMS( 3 ) = TDIM        !  time step dimension
      DSCBUF = '<YYYYDDD,HHMMSS>'
      DSCBU2 = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'


!...........   Create timestep-available-flag variable:
      IERR = NF90_DEF_VAR(FNUM, 'TFLAG', NF90_INT, DIMS(1:NDIMS), CHK%TINDX)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating netCDF variable TFLAG ' )
          GO TO 999
      END IF              !  ierr nonzero:  NCVDEF() failed

      IERR = NF90_PUT_ATT(FNUM, CHK%TINDX, 'units', DSCBUF)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating FLAG attribute "units"' )
          GO TO 999
      END IF

      DSCBUF = 'TFLAG'
      IERR = NF90_PUT_ATT(FNUM, CHK%TINDX, 'long_name', DSCBUF)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating FLAG attribute "long_name" ' )
          GO TO 999
      END IF

      IERR = NF90_PUT_ATT(FNUM, CHK%TINDX, 'var_desc', DSCBU2)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &          'Error creating FLAG attribute VAR_DESC ')
          GO TO 999
      END IF

!...........   Now, other dimensions (depend upon file type):
      IERR = NF90_DEF_DIM(FNUM, 'ROW', NROWS, RDIM)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &              'Error creating netCDF dimension ROW ' )
          GO TO 999
      END IF              !  ierr nonzero:  NCDDEF() failed

      IERR = NF90_DEF_DIM(FNUM, 'COL', NCOLS, CDIM)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &              'Error creating netCDF dimension COL' )
          GO TO 999
      END IF              !  ierr nonzero:  NCDDEF() failed

      DIMS( 1 ) = CDIM
      DIMS( 2 ) = RDIM
      DIMS( 3 ) = LDIM
      DIMS( 4 ) = TDIM

      NDIMS = 4

!.......   Define all the Models-3 variables for this file:

      DO VAR = 1 , CHK%NVARS

            ISCR = CHK%VTYPE( VAR )

            IF ( ISCR .NE. M3INT   .AND.
     &           ISCR .NE. M3REAL  .AND.
     &           ISCR .NE. M3DBLE ) THEN
                WRITE( LOGDEV,91010 )
     &          'Error creating netCDF variable ' // CHK%VNAME( VAR ) //
     &          'for file ' // FNAME,
     &          'Illegal data type', ISCR
                CALL NCABOR( FNUM, IERR )
                GO TO 999
            END IF

            IERR = NF90_DEF_VAR(FNUM, CHK%VNAME(VAR), CHK%VTYPE(VAR), DIMS(1:NDIMS), CHK%VINDX(VAR))
            IF ( IERR .NE. 0 ) THEN
                DSCBU2 =
     &           'Error creating netCDF variable ' // CHK%VNAME( VAR )
                CALL M3ABORT( FNAME, FNUM, IERR, DSCBU2 )
                GO TO 999
            END IF              !  ierr nonzero:  NCVDEF() failed

!...............   Now define attributes for the variables:
!...........   Note:  "long_name" and "units" are standard netCDF attributes.
            IERR = NF90_PUT_ATT(FNUM, CHK%VINDX(VAR), 'long_name', CHK%VNAME(VAR))
            IF ( IERR .NE. 0 ) THEN
                DSCBU2 =
     &          'Error creating  attribute LONG_NAME for variable ' //
     &          CHK%VNAME( VAR )
                CALL M3ABORT( FNAME, FNUM, IERR, DSCBU2 )
                GO TO 999
            END IF

            IERR = NF90_PUT_ATT(FNUM, CHK%VINDX(VAR), 'units', CHK%UNITS(VAR))
            IF ( IERR .NE. 0 ) THEN
                WRITE( LOGDEV,91010 )
     &              'Error creating ' // CHK%VNAME( VAR ) //
     &                  ' attribute UNITS for file ' // FNAME,
     &              'netCDF error number', IERR
                CALL NCABOR( FNUM, IERR )
                GO TO 999
            END IF

            IERR = NF90_PUT_ATT(FNUM, CHK%VINDX(VAR), 'var_desc', CHK%VDESC(VAR))
            IF ( IERR .NE. 0 ) THEN
                DSCBU2 =
     &          'Error creating  attribute UNITS for variable ' //
     &          CHK%VNAME( VAR )
                CALL M3ABORT( FNAME, FNUM, IERR, DSCBU2 )
                GO TO 999
            END IF

      END DO        ! end loop on VAR, defining variables & attributes

!...........   Put FNUM back into data mode:  attributes and variables now defined.
      IERR = NF90_ENDDEF(FNUM)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, FNUM, IERR,
     &              'Error putting netCDF file into data mode.' )
          GO TO 999
      END IF          !  ierr nonzero:  NCENDF() failed


999   CONTINUE


      RETURN

!...........   Error and warning message formats..... 91xxx

91010 FORMAT ( //5X , '>>> WARNING in subroutine CHK_CREATE_FILE <<<',
     &            2 ( /5X , A , : ) , I5, // )

      END SUBROUTINE CHK_FILE_CREATE_IMPL



!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_OPEN_IMPL (CHK_ID)
!-----------------------------------------------------------------------
!
!  FUNCTION:  Opens netCDF checkpoint file CHK_ID
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!
!  REVISION  HISTORY:
!   July 2010 by J.Resler, ICS: Initial version
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      USE SUBST_MODULES         ! stenex
      USE NETCDF
      IMPLICIT NONE

!...........   INCLUDES:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
!      INCLUDE SUBST_STATE
      INCLUDE SUBST_IOFDESC
!      INCLUDE SUBST_ATDSC

!  EXTERNAL FUNCTIONS and their descriptions

      INTEGER, EXTERNAL        :: ENVINT, SEC2TIME
      LOGICAL, EXTERNAL        :: FLUSH3     !  length after trimming trailing blanks.

!  ARGUMENTS and their descriptions:

      INTEGER                  :: CHK_ID !  identification of the checkpoint file

!  LOCAL VARIABLES and their descriptions:

      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_OPEN_IMPL'
      type(CHK_FILE),pointer   :: CHK
      CHARACTER(16)            :: FNAME         ! checkpoint file name
      CHARACTER(255)           :: FULLNAME      ! checkpoint file name incl. path
      INTEGER                  :: FMODE           ! netcdf creation mode
      INTEGER                  :: IERR
      CHARACTER( 96 )          :: XMSG = ' '
      INTEGER                   :: VAR

      INTEGER, SAVE           :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!.............................................................................
!   begin body of subroutine  CHK_FILE_OPEN_IMPL

      IF ( FIRSTIME ) THEN
         LOGDEV = INIT3 ()
         FIRSTIME = .FALSE.
      END IF

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)
      CALL CHK_FILE_NAME_IMPL (CHK_ID, FULLNAME)

!     open checkpoint file
      FMODE = NF90_NOWRITE
      IERR = NF90_OPEN(FULLNAME, FMODE, CHK%FNUM)

      IF ( IERR .NE. NF90_NOERR ) THEN
         XMSG = 'Could not open ' // TRIM(FNAME) // ' file'
         CALL M3EXIT( PNAME, SDATE3D, STIME3D, XMSG, XSTAT1 )
      END IF

!.......   TSTEP:  time step
      IERR = NF90_GET_ATT(CHK%FNUM, NF90_GLOBAL, 'TSTEP', CHK%TSTEP)
      IF ( IERR .NE. NF90_NOERR ) THEN
          CALL M3ABORT( FNAME, CHK%FNUM, IERR,
     &    'Error reading netCDF file attribute TSTEP' )
          GO TO 998
      END IF

      IERR = NF90_GET_ATT(CHK%FNUM, NF90_GLOBAL, 'SDATE', CHK%SDATE)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, CHK%FNUM, IERR,
     &        'Error reading netCDF file attribute SDATE' )
          GO TO 998
      END IF

      IERR = NF90_GET_ATT(CHK%FNUM, NF90_GLOBAL, 'STIME', CHK%STIME)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, CHK%FNUM, IERR,
     &    'Error creating netCDF file attribute STIME' )
          GO TO 998
      END IF

!......   TINDX:  index of time variable
      IERR = NF90_INQ_VARID(CHK%FNUM, 'TFLAG', CHK%TINDX)
      IF ( IERR .NE. 0 ) THEN
          CALL M3ABORT( FNAME, CHK%FNUM, IERR,
     &        'Error inquiring variable TFLAG ' )
          GO TO 998
      END IF

!......   VINDX:  index of variables
      DO  VAR = 1 , CHK%NVARS
         IERR = NF90_INQ_VARID(CHK%FNUM, CHK%VNAME(VAR), CHK%VINDX(VAR))
         IF ( IERR .NE. 0 ) THEN
             CALL M3ABORT( FNAME, CHK%FNUM, IERR,
     &           'Error inquiring variable '//TRIM(CHK%VNAME(VAR)) )
             GO TO 998
         END IF

      END DO

      write(logdev,*) 'end of open_impl'

998   CONTINUE

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
      USE NETCDF
      IMPLICIT NONE

!      INCLUDE SUBST_NETCDF

!     Parameters
      INTEGER                  :: CHK_ID

!     Local variables:
      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_CLOSE_IMPL'
      CHARACTER( 96 )          :: XMSG = ' '
      type(CHK_FILE),pointer   :: CHK
      CHARACTER(16)            :: FNAME           ! checkpoint file name
      INTEGER                  :: IERR


      INTEGER, SAVE           :: LOGDEV
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)

      IERR = NF90_CLOSE(CHK%FNUM)

      IF ( IERR .NE. NF90_NOERR ) THEN
         XMSG = 'Could not close ' // TRIM(FNAME)
         CALL M3WARN( PNAME, CHK%SDATE, CHK%STIME, XMSG )
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
!   July 2010 by J.Resler, ICS: Initial version
!   Oct  2010 by J.Resler, ICS: Adapted for integration into the new
!                               CMAQ adjoint version
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE NETCDF
      USE CHK_FILES_DEFS
      IMPLICIT NONE

! Include Files:
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

! EXTERNAL FUNCTIONS and their descriptions:
      INTEGER                 :: JSTEP3     !  compute time step record numbers

! Parameters
!      CHARACTER(len=NAMLEN_ADJ) :: FNAME
!     REAL         CGRID( NCOLS,NROWS,NLAYS,* )
      INTEGER                 :: CHK_ID       ! chk file id
      REAL                     :: CHK_GRID( :,:,:,: )
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS

! Local variables:
      type(CHK_FILE),pointer  :: CHK          ! chk file description struct
      CHARACTER(16)           :: FNAME           ! checkpoint file name
      INTEGER                 :: STEP         !
      INTEGER                 :: VAR
      INTEGER                 :: IERR
      INTEGER                 :: INDX

!      INTEGER                :: NDIMS           !  number of dims for NCVDEF()
      INTEGER                 :: DIMS( 4 )       !  array of dims for NCVDEF()
      INTEGER                 :: DELS( 4 )       !  array of dims for NCVPT()
      INTEGER                 :: DELTA
      REAL, POINTER           :: BUFFER( : )
      INTEGER                 :: C,R,L           ! loop variables

      CHARACTER(NAMLEN_ADJ)   :: PNAME = 'CHK_FILE_READ_IMPL'
      CHARACTER( 96 )         :: XMSG = ' '

      INTEGER, SAVE         :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE         :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      IF (FIRSTIME) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
      ENDIF

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)

      STEP = JSTEP3( JDATE, JTIME, CHK%SDATE, CHK%STIME, CHK%TSTEP )
      IF ( STEP .LT. 0 ) THEN
          XMSG = 'Time step error for file:  '//FNAME
          CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
          RETURN
      END IF

      DIMS( 1 ) = 1
      DELS( 1 ) = NCOLS
      DIMS( 2 ) = 1
      DELS( 2 ) = NROWS
      DIMS( 3 ) = 1
      DELS( 3 ) = CHK%NLAYS
      DIMS( 4 ) = STEP
      DELS ( 4 ) = 1
      ALLOCATE(BUFFER(NCOLS*NROWS*CHK%NLAYS))

      write(logdev,*) 'Reading chk_file: ', CHK%FNUM, FNAME, JDATE, JTIME, STEP

      DO VAR = 1, CHK%NVARS

!          CALL NCVGT( CHK%FNUM, CHK%VINDX( VAR ),DIMS, DELTS, BUFFER, IERR )
          IERR = NF90_GET_VAR(CHK%FNUM, CHK%VINDX( VAR ), BUFFER, DIMS, DELS)
!          IERR = NF_GET_VARA_REAL(CHK%FNUM, CHK%VINDX( VAR ), DIMS, DELS, BUFFER)

          IF ( IERR .NE. NF90_NOERR ) THEN
              WRITE( XMSG,'( A, I6, 2X, 5 A )' )
     &             'netCDF error number', IERR,
     &             'Error reading variable "', CHK%VNAME( VAR ),
     &             '" from file "', FNAME, '"'
              CALL M3MSG2( XMSG )
              DEALLOCATE(BUFFER)
              RETURN
          END IF

          INDX = 0
          DO L = 1, CHK%NLAYS
            DO R = 1, NROWS
              DO C = 1, NCOLS
                INDX = INDX + 1
                CHK_GRID(C,R,L,VAR) = BUFFER(INDX)
              END DO
            END DO
          END DO

      END DO    !  end loop on variables VAR

      DEALLOCATE(BUFFER)

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
!     Resler Jaroslav, ICS AS Prague, Aug 2001 - Initial code
!-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      USE NETCDF

      IMPLICIT NONE

! Include Files:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
!      INCLUDE SUBST_NETCDF

      INTEGER                 :: CHK_ID       ! chk file id
      REAL                     :: CHK_GRID( :,:,:,: )
      INTEGER                 :: JDATE        ! current model date, coded YYYYDDD
      INTEGER                 :: JTIME        ! current model time, coded HHMMSS

C...........   EXTERNAL FUNCTIONS and their descriptions:
      INTEGER                 :: JSTEP3     !  compute time step record numbers

! Local variables:
      type(CHK_FILE),pointer  :: CHK          ! chk file description struct
      CHARACTER(16)           :: FNAME           ! checkpoint file name
      INTEGER                 :: STEP         !
      INTEGER                 :: VAR
      INTEGER                 :: IERR
      INTEGER                 :: INDX

      INTEGER                 :: DIMS( 4 )       !  array of dims for NCVDEF()
      INTEGER                 :: DELS( 4 )       !  array of dims for NCVPT()
      INTEGER                 :: DIMTS( 3 )      !  array of dims for NCVDEF()
      INTEGER                 :: DELTS( 3 )      !  array of dims for NCVPT()
      INTEGER                 :: FLAGS( 2 )      !  ( jdate yyyyddd, jtime hhmmss )
      REAL, POINTER           :: BUFFER( : )
      INTEGER                 :: C,R,L           ! loop variables

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

      STEP = JSTEP3( JDATE, JTIME, CHK%SDATE, CHK%STIME, CHK%TSTEP )
      IF ( STEP .LT. 0 ) THEN
          XMSG = 'Time step error for file:  '//CHK%FNAME
          CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
          RETURN
      END IF

      DIMS( 1 ) = 1
      DELS( 1 ) = NCOLS
      DIMS( 2 ) = 1
      DELS( 2 ) = NROWS
      DIMS( 3 ) = 1
      DELS( 3 ) = CHK%NLAYS
      DIMS( 4 ) = STEP
      DELS( 4 ) = 1

      FLAGS( 1 ) = JDATE
      FLAGS( 2 ) = JTIME
      DIMTS( 1 ) = 1          !  field:  date or time
      DELTS( 1 ) = 2          !  extent:  both date and time
      DIMTS( 2 ) = 1          !  variable-number
      DELTS( 2 ) = 1          !  extent:  1
      DIMTS( 3 ) = STEP       !  time start
      DELTS( 3 ) = 1          !  time extent:  1

      write(logdev,*) 'Writing chk_file: ', CHK%FNUM, FNAME, JDATE, JTIME, STEP
!...............   Write contents of the buffer for each variable

      ALLOCATE(BUFFER(NCOLS*NROWS*CHK%NLAYS))
      DO  VAR = 1, CHK%NVARS
          INDX = 0
          DO L = 1, CHK%NLAYS
            DO R = 1, NROWS
              DO C = 1, NCOLS
                INDX = INDX + 1
                BUFFER(INDX) = CHK_GRID(C,R,L,VAR)
              END DO
            END DO
          END DO

          IERR = NF90_PUT_VAR(CHK%FNUM, CHK%VINDX( VAR ), BUFFER, DIMS, DELS)

          IF ( IERR .NE. NF90_NOERR ) THEN
!$OMP CRITICAL( S_LOGOUT )
              XMSG = 'Error writing variable '
     &               // CHK%VNAME( VAR ) //
     &               'to file ' // FNAME
              WRITE( LOGDEV,* ) TRIM( XMSG )
              WRITE( LOGDEV,* ) 'netCDF error number', IERR
              WRITE( LOGDEV,* ) 'CHK ID      ', CHK_ID
              WRITE( LOGDEV,* ) 'netCDF ID   ', CHK%FNUM
              WRITE( LOGDEV,* ) 'vble        ', CHK%VINDX( VAR )
              WRITE( LOGDEV,* ) 'dims array  ', DIMS
              WRITE( LOGDEV,* ) 'delts array ', DELS
              WRITE( LOGDEV,* ) 'offset      ', INDX
              WRITE( LOGDEV,* )
!$OMP END CRITICAL( S_LOGOUT )

               DEALLOCATE(BUFFER)
               RETURN

          END IF

          IERR = NF90_PUT_VAR(CHK%FNUM, CHK%TINDX, FLAGS, DIMTS, DELTS)
          IF ( IERR .NE. 0 ) THEN
!$OMP CRITICAL( S_LOGOUT )
              WRITE( LOGDEV,* )
     &          'FDV_WR_CHK:  Error writing time-flag to file ' // FNAME
	          WRITE( LOGDEV,* ) 'netCDF error number', IERR
	          WRITE( LOGDEV,* ) 'chk file ID ', CHK_ID
	          WRITE( LOGDEV,* ) 'netCDF ID     ', CHK%FNUM
	          WRITE( LOGDEV,* ) 'dims array    ', DIMTS
	          WRITE( LOGDEV,* ) 'delts array   ', DELTS
	          WRITE( LOGDEV,* )
!$OMP END CRITICAL( S_LOGOUT )
              RETURN
          END IF

      END DO

      DEALLOCATE(BUFFER)
      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &         'Timestep written to', FNAME,
     &         'for date and time', JDATE, JTIME

      RETURN

      END SUBROUTINE CHK_FILE_WRITE_IMPL


!-----------------------------------------------------------------------
      SUBROUTINE CHK_FILE_NAME_IMPL (CHK_ID, FULLNAME)
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

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS
      IMPLICIT NONE

!     Parameters
      INTEGER, INTENT(IN)        :: CHK_ID       ! chk file id
      CHARACTER(*), INTENT(OUT)  :: FULLNAME


!     Local variables:
      type(CHK_FILE),pointer   :: CHK          ! chk file description struct
      CHARACTER(16)            :: FNAME
      CHARACTER(NAMLEN_ADJ)    :: PNAME = 'CHK_FILE_NAME_IMPL'
      CHARACTER( 96 )          :: XMSG = ' '
      INTEGER                  :: IERR

      INTEGER, SAVE           :: LOGDEV
      LOGICAL, SAVE           :: FIRSTIME = .TRUE.

!-----------------------------------------------------------------------

      CHK => CHKS(CHK_ID)
      FNAME = CHK_NAMES(CHK_ID)
      CALL CHK_FILE_PATH ( FULLNAME )

!   Create the real file name
      IF ( TRIM(FULLNAME) .EQ. '' ) THEN
         WRITE(LOGDEV,92001) TRIM(FNAME),CHK%SDATE, MYPE
         WRITE(FULLNAME,92001) TRIM(FNAME),CHK%SDATE, MYPE
92001    FORMAT (  A , "_", I7, "_", I3.3 )
      ELSE
         WRITE(LOGDEV,92002) TRIM(FULLNAME),TRIM(FNAME),CHK%SDATE, MYPE
         WRITE(FULLNAME,92002) TRIM(FULLNAME),TRIM(FNAME),CHK%SDATE, MYPE
92002    FORMAT (  A , "/", A , "_", I7, "_", I3.3 )
      END IF


      RETURN

      END SUBROUTINE CHK_FILE_NAME_IMPL




      END MODULE CHK_FILES_IMPL

