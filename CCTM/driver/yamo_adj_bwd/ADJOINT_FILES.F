C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      MODULE ADJOINT_FILES

C-----------------------------------------------------------------------
C Function:
C   Define adjoint file functions

C Revision History:
C   Aug 2007 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C   Oct 2010 J.Resler ICS: checkpoint files moved to CHK_FILES module
C   Feb 2011 by Peter Percell at UH-IMAQS: removed unneeded code
C   May 2013 by Peter Percell at UH-IMAQS:
C   - the adjoint file implementation here uses IOAPI3 tools, so a constant
C     output step is required
C   Aug 2013 by Peter Percell at UH-IMAQS:
C   - made the LGRID output frequency selectable by user
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

C-----------------------------------
C Data
C-----------------------------------

C Parameters
      CHARACTER(16), PARAMETER :: CTM_XFIRST_IN   = 'CTM_XFIRST_IN'
      CHARACTER(16), PARAMETER :: CTM_XFIRST_OUT  = 'CTM_XFIRST_OUT'

!      CHARACTER(16), PARAMETER :: CTM_CONC_FWD    = 'CTM_CONC_FWD'
      CHARACTER(16) :: CTM_CONC_FWD
      CHARACTER(16), PARAMETER :: INIT_LGRID_1    = 'INIT_LGRID_1'
      CHARACTER(16), PARAMETER :: INIT_EM_1       = 'INIT_EM_1'
      CHARACTER(16), PARAMETER :: INIT_EM_SF_1    = 'INIT_EM_SF_1'
      CHARACTER(16), PARAMETER :: ADJ_LGRID       = 'ADJ_LGRID'
      CHARACTER(16), PARAMETER :: ADJ_FORCE_FILE  = 'ADJ_FORCE'
      CHARACTER(16), PARAMETER :: ADJ_LGRID_EM_SF = 'ADJ_LGRID_EM_SF'
      CHARACTER(16), PARAMETER :: ADJ_LGRID_EM    = 'ADJ_LGRID_EM'
      CHARACTER(16), PARAMETER :: ADJ_LGRID_EM_NRM= 'ADJ_LGRID_EM_NRM'

      CHARACTER(16) :: LGRID_OUTPUT_FREQ

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      CONTAINS

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Create new adjoint file
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE OP_ADJ_FILE(JDATE, JTIME, TSTEP, ADJ_FILE)

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE EMIS_DEFN, ONLY : EMSENSL

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_DIFF     ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF     ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF     ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF     ! tracer diffusion species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Parameters:

C Arguments:
      INTEGER, INTENT(IN) :: JDATE ! starting date (YYYYDDD)
      INTEGER, INTENT(IN) :: JTIME ! starting time (HHMMSS)
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                !   TSTEP(1) = local output step
                                !   TSTEP(2) = sciproc sync. step (chem)
      CHARACTER(16), INTENT(IN) :: ADJ_FILE

C Local Variables:
      INTEGER, SAVE :: LOGDEV   ! FORTRAN unit number for log file

      CHARACTER(16) :: PNAME = 'OP_ADJ_FILE'
      CHARACTER(96) :: XMSG = ' '

      CHARACTER(80) :: VARDESC  ! environment variable description

      INTEGER :: STATUS

      INTEGER :: K, S           ! loop counters

C-----------------------------------------------------------------------

      LOGDEV = INIT3()

      VARDESC =
     &     'LGRID output frequency: SYNC_STEP, OUTPUT_STEP or ONCE (at start time)'
      LGRID_OUTPUT_FREQ = 'SYNC_STEP'
      CALL ENVSTR('ADJ_LGRID_FREQ', VARDESC, LGRID_OUTPUT_FREQ,
     &     LGRID_OUTPUT_FREQ, STATUS)
      IF ( STATUS /= 0 ) WRITE(LOGDEV, '(5X, A)') VARDESC
      IF ( STATUS == 1 ) THEN
         XMSG = 'Environment variable improperly formatted'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
      ELSE IF ( STATUS == -1 ) THEN
         XMSG = 'Environment variable set, but empty ... Using default'
         WRITE(LOGDEV, '(5X, A)') XMSG
      ELSE IF ( STATUS == -2 ) THEN
         XMSG = 'Environment variable not set ... Using default'
         WRITE(LOGDEV, '(5X, A)') XMSG
      END IF

      IF ( TRIM(LGRID_OUTPUT_FREQ) /= 'SYNC_STEP' .AND.
     &     TRIM(LGRID_OUTPUT_FREQ) /= 'OUTPUT_STEP' .AND.
     &     TRIM(LGRID_OUTPUT_FREQ) /= 'ONCE' ) THEN
         WRITE(LOGDEV, '(5X, A)')
     &        'ADJ_LGRID_FREQ was set to ' // LGRID_OUTPUT_FREQ
         XMSG = 'It must be SYNC_STEP, OUTPUT_STEP or ONCE'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)
      END IF

      IF ( .NOT. CONST_SYNC_STEP .AND.
     &     TRIM(LGRID_OUTPUT_FREQ) == 'SYNC_STEP' ) THEN
         XMSG = 'ADJ_LGRID_FREQ cannot be SYNC_STEP ... Using OUTPUT_STEP'
         CALL M3WARN(PNAME, JDATE, JTIME, XMSG)
         LGRID_OUTPUT_FREQ = 'OUTPUT_STEP'
      END IF

      IF ( MYPE == 0 ) THEN
!
!         IF ( ADJ_FILE == ADJ_LGRID ) THEN
!
!C Get file characteristics of CTM_CONC_FWD
!            IF ( .NOT. OPEN3(CTM_CONC_FWD, FSREAD3, PNAME) ) THEN
!               XMSG = 'Could not open ' // TRIM(CTM_CONC_FWD) //
!     &              ' file as OLD:READ-ONLY'
!               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
!            END IF
!
!            IF ( .NOT. DESC3(CTM_CONC_FWD) ) THEN
!               XMSG = 'Could not get ' // TRIM(CTM_CONC_FWD) // ' file description'
!               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
!            END IF
!
!            IF ( .NOT. CLOSE3(CTM_CONC_FWD) ) THEN
!               XMSG = 'Could not close ' // TRIM(CTM_CONC_FWD)
!               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
!            END IF
!
!         ELSE

C Get file characteristics of CTM_CONC_1
            IF ( .NOT. OPEN3(CTM_CONC_1, FSREAD3, PNAME) ) THEN
               XMSG = 'Could not open ' // TRIM(CTM_CONC_1) //
     &              ' file as OLD:READ-ONLY'
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
            END IF

            IF ( .NOT. DESC3(CTM_CONC_1) ) THEN
               XMSG = 'Could not get ' // TRIM(CTM_CONC_1) // ' file description'
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
            END IF

            IF ( .NOT. CLOSE3(CTM_CONC_1) ) THEN
               XMSG = 'Could not close ' // TRIM(CTM_CONC_1)
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
            END IF

!         END IF

C Set file characteristics that differ from the base file, CTM_CONC_FWD or
C CTM_CONC_1
         SELECT CASE ( ADJ_FILE )

         CASE ( ADJ_LGRID )

            IF ( TRIM(LGRID_OUTPUT_FREQ) == 'SYNC_STEP' ) THEN
               TSTEP3D = TSTEP(2) ! sync step
            ELSE IF ( TRIM(LGRID_OUTPUT_FREQ) == 'OUTPUT_STEP' ) THEN
               TSTEP3D = TSTEP(1) ! output step
            ELSE
               TSTEP3D = 0      ! make it a time-independent file
            END IF
            
            FDESC3D(1) = 'Sensitivity output'
            FDESC3D(2) = 'From CMAQ model dyn alloc version CTM'
            DO K = 3, MXDESC3
               FDESC3D(K) = ' '
            END DO

            NLAYS3D = NLAYS
            NVARS3D = 0

C Note: "CF" stands for "cost function units"
            DO S = 1, N_GC_SPC
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = GC_SPC(S)
               UNITS3D(NVARS3D) = 'CF/ppmV'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

            NVARS3D = NVARS3D + 1
            VTYPE3D(NVARS3D) = M3REAL
            VNAME3D(NVARS3D) = 'RHOJ'
            UNITS3D(NVARS3D) = 'CF/m*kg/m**3'
            VDESC3D(NVARS3D) = 'Sensitivity of cost function to ' //
     &           TRIM(VNAME3D(NVARS3D))

!     Aerosol species - slc.4.2014
            DO S = 1, N_AE_SPC
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = AE_SPC(S)
               IF ( VNAME3D(NVARS3D)(1:3) == 'NUM' ) THEN
                  UNITS3D(NVARS3D) = 'CF/number/m**3'
               ELSE IF ( VNAME3D(NVARS3D)(1:3) == 'SRF' ) THEN
                  UNITS3D(NVARS3D) = 'CF/m**2/m**3'
               ELSE
                  UNITS3D(NVARS3D) = 'CF/microg/m**3'
               END IF
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Ammonia and semi-volatile compounds
            DO S = 1, N_NR_SPC
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = NR_SPC(S)
               UNITS3D(NVARS3D) = 'CF/ppmV'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Tracer species - slc.4.2014
            DO S = 1, N_TR_SPC
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = TR_SPC(S)
               UNITS3D(NVARS3D) = 'CF/ppmV'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

         CASE ( ADJ_LGRID_EM )

            TSTEP3D = TSTEP(1)  ! output step

            FDESC3D(1) = 'Sensitivity output wrt emissions'
            FDESC3D(2) = 'From CMAQ model dyn alloc version CTM'
            DO K = 3, MXDESC3
               FDESC3D(K) = ' '
            END DO

            NLAYS3D = EMSENSL
            NVARS3D = 0

            IF ( TRIM(LGRID_OUTPUT_FREQ) == 'SYNC_STEP' ) THEN
               TSTEP3D = TSTEP(2) ! sync step
            ELSE IF ( TRIM(LGRID_OUTPUT_FREQ) == 'OUTPUT_STEP' ) THEN
               TSTEP3D = TSTEP(1) ! output step
            ELSE
               TSTEP3D = 0      ! make it a time-independent file
            END IF

C Note: "CF" stands for "cost function units"
!       Only diffused species are calculated as emissions currently. so
!       the loop is over the diffused species. 

            ! dbg.1.2015.slc
            !write(*,*) 'N_*_DIFF, GC: ',N_GC_DIFF,' AE: ',N_AE_DIFF
            !write(*,*) 'NR: ',N_NR_DIFF,' TR: ',N_TR_DIFF 
            DO S = 1, N_GC_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = GC_DIFF(S)
               UNITS3D(NVARS3D) = 'CF/(ppmv/s)'
               VDESC3D(NVARS3D) = 
     &              'Sensitivity of cost function to emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!            NVARS3D = NVARS3D + 1
!            VTYPE3D(NVARS3D) = M3REAL
!            VNAME3D(NVARS3D) = 'RHOJ'
!            UNITS3D(NVARS3D) = 'CF/m*kg/m**3'
!            VDESC3D(NVARS3D) = 
!     &           'Sensitivity of cost function to emissions of  ' //
!     &           TRIM(VNAME3D(NVARS3D))

!     Aerosol species - slc.7.2014
            DO S = 1, N_AE_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = AE_DIFF(S)
               IF ( VNAME3D(NVARS3D)(1:3) == 'NUM' ) THEN
                  UNITS3D(NVARS3D) = 'CF/(g/s)'
               ELSE IF ( VNAME3D(NVARS3D)(1:3) == 'SRF' ) THEN
                  UNITS3D(NVARS3D) = 'CF/(g/s)'
               ELSE
                  UNITS3D(NVARS3D) = 'CF/(g/s)'
               END IF
               VDESC3D(NVARS3D) = 
     &              'Sensitivity of cost function to emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Ammonia and semi-volatile compounds
            DO S = 1, N_NR_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = NR_DIFF(S)
               UNITS3D(NVARS3D) = 'CF/(mol/s)'
               VDESC3D(NVARS3D) = 
     &              'Sensitivity of cost function to emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Tracer species - slc.7.2014
            DO S = 1, N_TR_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = TR_DIFF(S)
               UNITS3D(NVARS3D) = 'CF/(mol/s)'
               VDESC3D(NVARS3D) = 
     &              'Sensitivity of cost function to emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO


         CASE ( ADJ_LGRID_EM_SF )

            TSTEP3D = TSTEP(1)  ! output step

            FDESC3D(1) = 'Sensitivity output wrt emissions scaling  ' //
     &                    ' factor'
            FDESC3D(2) = 'From CMAQ model dyn alloc version CTM'
            DO K = 3, MXDESC3
               FDESC3D(K) = ' '
            END DO

            NLAYS3D = EMSENSL ! only surface value
            NVARS3D = 0

            IF ( TRIM(LGRID_OUTPUT_FREQ) == 'SYNC_STEP' ) THEN
               TSTEP3D = TSTEP(2) ! sync step
            ELSE IF ( TRIM(LGRID_OUTPUT_FREQ) == 'OUTPUT_STEP' ) THEN
               TSTEP3D = TSTEP(1) ! output step
            ELSE
               TSTEP3D = 0      ! make it a time-independent file
            END IF

C Note: "CF" stands for "cost function units"
            DO S = 1, N_GC_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = GC_DIFF(S)
               UNITS3D(NVARS3D) = 'CF'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to
     %              emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!            NVARS3D = NVARS3D + 1
!            VTYPE3D(NVARS3D) = M3REAL
!            VNAME3D(NVARS3D) = 'RHOJ'
!            UNITS3D(NVARS3D) = 'CF'
!            VDESC3D(NVARS3D) = 'Sensitivity of cost function to
!     &           emissions of  ' //
!     &           TRIM(VNAME3D(NVARS3D))

!     Aerosol species - slc.7.2014
            DO S = 1, N_AE_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = AE_DIFF(S)
               IF ( VNAME3D(NVARS3D)(1:3) == 'NUM' ) THEN
                  UNITS3D(NVARS3D) = 'CF'
               ELSE IF ( VNAME3D(NVARS3D)(1:3) == 'SRF' ) THEN
                  UNITS3D(NVARS3D) = 'CF'
               ELSE
                  UNITS3D(NVARS3D) = 'CF'
               END IF
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to
     &              emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Ammonia and semi-volatile compounds
            DO S = 1, N_NR_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = NR_DIFF(S)
               UNITS3D(NVARS3D) = 'CF'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to
     &              emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Tracer species - slc.7.2014
            DO S = 1, N_TR_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = TR_DIFF(S)
               UNITS3D(NVARS3D) = 'CF'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to
     &              emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO


         CASE ( ADJ_LGRID_EM_NRM )

            TSTEP3D = TSTEP(1)  ! output step

            FDESC3D(1) = 'Sensitivity output wrt emissions - fully' //
     &                    ' normalized'
            FDESC3D(2) = 'From CMAQ model dyn alloc version CTM'
            DO K = 3, MXDESC3
               FDESC3D(K) = ' '
            END DO

            NLAYS3D = EMSENSL ! NLAYS
            NVARS3D = 0

            IF ( TRIM(LGRID_OUTPUT_FREQ) == 'SYNC_STEP' ) THEN
               TSTEP3D = TSTEP(2) ! sync step
            ELSE IF ( TRIM(LGRID_OUTPUT_FREQ) == 'OUTPUT_STEP' ) THEN
               TSTEP3D = TSTEP(1) ! output step
            ELSE
               TSTEP3D = 0      ! make it a time-independent file
            END IF

C Note: "CF" stands for "cost function units"
            DO S = 1, N_GC_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = GC_SPC(S)
               UNITS3D(NVARS3D) = 'CF/(mol/s)'
               VDESC3D(NVARS3D) = 'Normalized ensitivity of cost function to emissions of ' // 
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!            NVARS3D = NVARS3D + 1
!            VTYPE3D(NVARS3D) = M3REAL
!            VNAME3D(NVARS3D) = 'RHOJ'
!            UNITS3D(NVARS3D) = 'CF/m*kg/m**3'
!            VDESC3D(NVARS3D) = 'Normalized ensitivity of cost function to emissions of ' //
!     &           TRIM(VNAME3D(NVARS3D))

!     Aerosol species - slc.7.2014
            DO S = 1, N_AE_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = AE_DIFF(S)
               IF ( VNAME3D(NVARS3D)(1:3) == 'NUM' ) THEN
                  UNITS3D(NVARS3D) = 'CF/(g/s)'
               ELSE IF ( VNAME3D(NVARS3D)(1:3) == 'SRF' ) THEN
                  UNITS3D(NVARS3D) = 'CF/(g/s)'
               ELSE
                  UNITS3D(NVARS3D) = 'CF/(g/s)'
               END IF
               VDESC3D(NVARS3D) = 'Normalized ensitivity of cost function to emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Ammonia and semi-volatile compounds
            DO S = 1, N_NR_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = NR_DIFF(S)
               UNITS3D(NVARS3D) = 'CF/(mol/s)'
               VDESC3D(NVARS3D) = 'Normalized ensitivity of cost function to emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

!        Tracer species - slc.7.2014
            DO S = 1, N_TR_DIFF
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = TR_DIFF(S)
               UNITS3D(NVARS3D) = 'CF/(mol/s)'
               VDESC3D(NVARS3D) = 'Normalized ensitivity of cost function to emissions of ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

         CASE DEFAULT

            XMSG =  '"' // TRIM(ADJ_FILE) // '"' //
     &           ' is not a valid file selection'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT3)

         END SELECT

C Create new file
         IF ( .NOT. OPEN3(ADJ_FILE, FSCREA3, PNAME) ) THEN
            XMSG = 'Could not create ' // TRIM(ADJ_FILE) // ' file'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

      END IF                    ! MYPE == 0

      END SUBROUTINE OP_ADJ_FILE

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Write to an adjoint file
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE WR_ADJ_FILE(ARRAY, JDATE, JTIME, ADJ_FILE)

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
#include      SUBST_IODECL      # I/O definitions and declarations

C Arguments:
      REAL, POINTER :: ARRAY(:,:,:,:)
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      CHARACTER(16), INTENT(IN) :: ADJ_FILE

C Local variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      CHARACTER(16) :: PNAME = 'WR_ADJ_FILE'
      CHARACTER(96) :: XMSG = ' '

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
      END IF

      !WRITE(LOGDEV,*) 'WR_ADJ_FILE, sum(ARRAY): ',sum(ARRAY)
      !WRITE(LOGDEV,*) 'WR_ADJ_FILE, shape(ARRAY): ',shape(ARRAY)
      !WRITE(LOGDEV,*) 'JDATE.TIME: ',JDATE, JTIME

      IF ( .NOT. WRITE3(ADJ_FILE, ALLVAR3, JDATE, JTIME, ARRAY) ) THEN
         XMSG = 'Could not write to ' // TRIM(ADJ_FILE)
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      WRITE(LOGDEV, '(/5X, 3(A, :, 1X), I8, ":", I6.6)')
     &     'Timestep written to', ADJ_FILE,
     &     'for date and time', JDATE, JTIME

      END SUBROUTINE WR_ADJ_FILE

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Read from an adjoint file
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_ADJ_FILE(ARRAY, JDATE, JTIME, ADJ_FILE)

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE VDIFF_MAP

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC     ! file header data structure
      REAL, INTENT(OUT) :: ARRAY(*)
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      CHARACTER(16), INTENT(IN) :: ADJ_FILE

C Local variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: PNAME = 'RD_ADJ_FILE'
      CHARACTER(96) :: XMSG = ' '

      REAL :: BUFFER(1:NCOLS, 1:NROWS, 1:NLAYS)

      INTEGER, SAVE :: STRTCOL, ENDCOL, STRTROW, ENDROW ! local window in global domain

      INTEGER :: K, C, R, L, V  ! loop counters

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

C Get local window in global domain
         STRTCOL = COLSX_PE(1, MYPE+1)
         ENDCOL = STRTCOL + NCOLS - 1
         STRTROW = ROWSX_PE(1, MYPE+1)
         ENDROW = STRTROW + NROWS - 1
      END IF

C Get file characteristics of ADJ_FILE
      IF ( .NOT. DESC3(ADJ_FILE) ) THEN
         XMSG = 'Could not get ' // TRIM(ADJ_FILE) // ' file description'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      K = 0
      DO V = 1, N_SPC_DIFF
         IF ( .NOT. XTRACT3(ADJ_FILE, VNAME3D(V),
     &        1, NLAYS, STRTROW, ENDROW, STRTCOL, ENDCOL,
     &        JDATE, JTIME, BUFFER) ) THEN
            XMSG = 'Could not read  ' // VNAME3D(V)
     &           // ' from ' // TRIM(ADJ_FILE)
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

         DO L = 1, NLAYS
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  K = K + 1
                  ARRAY(K) = BUFFER(C, R, L)
               END DO
            END DO
         END DO
      END DO

      END SUBROUTINE RD_ADJ_FILE



C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Read from an adjoint file
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_ADJ_EM_FILE(ARRAY, JDATE, JTIME, ADJ_FILE)

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE VDIFF_MAP
      USE EMIS_DEFN, ONLY : EMSENSL ! layers of emissions sensitivity 

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC     ! file header data structure
      REAL, INTENT(OUT) :: ARRAY(*)
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      CHARACTER(16), INTENT(IN) :: ADJ_FILE

C Local variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER(16) :: PNAME = 'RD_ADJ_EM_FILE'
      CHARACTER(96) :: XMSG = ' '

      REAL :: BUFFER(1:NCOLS, 1:NROWS, 1:EMSENSL)

      INTEGER, SAVE :: STRTCOL, ENDCOL, STRTROW, ENDROW ! local window in global domain

      INTEGER :: K, C, R, L, V  ! loop counters

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

C Get local window in global domain
         STRTCOL = COLSX_PE(1, MYPE+1)
         ENDCOL = STRTCOL + NCOLS - 1
         STRTROW = ROWSX_PE(1, MYPE+1)
         ENDROW = STRTROW + NROWS - 1
      END IF

C Get file characteristics of ADJ_FILE
      IF ( .NOT. DESC3(ADJ_FILE) ) THEN
         XMSG = 'Could not get ' // TRIM(ADJ_FILE) // ' file description'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

      K = 0
      DO V = 1, N_SPC_DIFF
         IF ( .NOT. XTRACT3(ADJ_FILE, VNAME3D(V),
     &        1, EMSENSL, STRTROW, ENDROW, STRTCOL, ENDCOL,
     &        JDATE, JTIME, BUFFER) ) THEN
            XMSG = 'Could not read  ' // VNAME3D(V)
     &           // ' from ' // TRIM(ADJ_FILE)
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

         DO L = 1, EMSENSL
            DO R = 1, NROWS
               DO C = 1, NCOLS
                  K = K + 1
                  ARRAY(K) = BUFFER(C, R, L)
               END DO
            END DO
         END DO
      END DO

      END SUBROUTINE RD_ADJ_EM_FILE



C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Read from an offline forcing file
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RD_FORCE_FILE(ARRAY, JDATE, JTIME) !, ADJ_FORCE)

      USE CGRID_SPCS            ! CGRID species number and offsets
      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_FILES_ID    ! file name parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table

C Arguments:

      INTEGER, EXTERNAL :: ENVINT
      INTEGER, EXTERNAL :: INDEX1
      INTEGER           :: IOS
      INTEGER           :: NLAYS_FRC
      INTEGER           :: ALLOCSTAT

      REAL, INTENT(OUT) :: ARRAY(:,:,:,:)
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
!      CHARACTER(16), INTENT(IN) :: ADJ_FORCE_FILE

C Local variables:
      LOGICAL, SAVE :: FIRSTIME_FRC = .TRUE.

      CHARACTER(16) :: PNAME = 'RD_FORCE_FILE'
      CHARACTER(96) :: XMSG = ' '
      
      INTEGER, SAVE :: NVARS_FRC
      INTEGER, ALLOCATABLE, SAVE :: VMAP_FRC(:)

!      INTEGER, PARAMETER :: NLAYS_FRC = 1

!      REAL :: BUFFER(1:NCOLS, 1:NROWS, 1:NLAYS_FRC)
!      REAL :: BUFFER(1:NCOLS, 1:NROWS)
      REAL, ALLOCATABLE :: BUFFER(:, :, :)

      INTEGER :: STRTCOL, ENDCOL, STRTROW, ENDROW ! local window in global domain
      INTEGER :: CARRAY, RARRAY, LENROW, LENCOL


      INTEGER :: C, R, L, V, S  ! loop counters

C-----------------------------------------------------------------------

      IF ( FIRSTIME_FRC ) THEN
         FIRSTIME_FRC = .FALSE.
         IF ( .NOT. DESC3(ADJ_FORCE_FILE) ) THEN
            XMSG = 'Could not get ' // TRIM(ADJ_FORCE_FILE) // ' file description'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF
         NVARS_FRC = NVARS3D
         ALLOCATE( VMAP_FRC(NVARS_FRC), STAT=ALLOCSTAT )
         DO V = 1, NVARS_FRC
            S = INDEX1( VNAME3D(V), N_GC_SPC, GC_SPC )
            IF ( S .EQ. 0 ) THEN
               S = INDEX1( VNAME3D(V), N_AE_SPC, AE_SPC )
               IF ( S .EQ. 0 ) THEN
                  S = INDEX1( VNAME3D(V), N_NR_SPC, NR_SPC )
                  IF ( S .EQ. 0 ) THEN
                     S = INDEX1( VNAME3D(V), N_TR_SPC, TR_SPC )
                     IF ( S .EQ. 0 ) THEN
                        XMSG = 'Could not find ' // VNAME3D(V) // ' in LGRID'
                        CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
                     ELSE
                        S = TR_STRT - 1 + S
                     END IF
                  ELSE
                     S = NR_STRT - 1 + S
                  END IF
               ELSE
                  S = AE_STRT - 1 + S
               END IF
            ELSE
               S = GC_STRT - 1 + S
            END IF
            VMAP_FRC(V) = S
         END DO
      END IF !firstime_frc

         NLAYS_FRC = ENVINT('NLAYS_FRC', 'No. of layers in forcing file', 1, IOS)
         ALLOCATE( BUFFER(1:NCOLS, 1:NROWS, 1:NLAYS_FRC), STAT=ALLOCSTAT )

C Get local window in global domain
         STRTCOL = COLSX_PE(1, MYPE+1)
         ENDCOL = STRTCOL + NCOLS - 1
         STRTROW = ROWSX_PE(1, MYPE+1)
         ENDROW = STRTROW + NROWS - 1
         
         BUFFER = 0.d0
         ARRAY  = 0.d0

C Get file characteristics of ADJ_FILE_FILE
      IF ( .NOT. DESC3(ADJ_FORCE_FILE) ) THEN
         XMSG = 'Could not get ' // TRIM(ADJ_FORCE_FILE) // ' file description'
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
      END IF

!     Debug (6.2014.slc)
      !WRITE(*,*) 'NVARS3D: ',NVARS3D
      !WRITE(*,*) 'VNAME3D(1): ',VNAME3D(1)
      !WRITE(*,*) 'FTYPE3D: ',FTYPE3D
      !WRITE(*,*) 'SDATE3D: ',SDATE3D
      !WRITE(*,*) 'ADJOINT files'
      !WRITE(*,*) 'NCOLS: ',NCOLS,' NROWS: ',NROWS
      !WRITE(*,*) 'NLAYS: ',NLAYS,'NSPCSD: ',NSPCSD
      !WRITE(*,*) 'STRTCOL: ',STRTCOL
      !WRITE(*,*) 'ENDCOL: ',ENDCOL
      !WRITE(*,*) 'STRTROW: ',STRTROW
      !WRITE(*,*) 'ENDROW: ',ENDROW

      !WRITE(*,*) 'JDATE: ',JDATE
      !WRITE(*,*) 'JTIME: ',JTIME

      DO V = 1, NVARS_FRC
         IF ( .NOT. XTRACT3(ADJ_FORCE_FILE, VNAME3D(V),
     &        1, NLAYS_FRC, STRTROW, ENDROW, STRTCOL, ENDCOL,
     &        JDATE, JTIME, BUFFER) ) THEN
            XMSG = 'Could not read  ' // VNAME3D(1)
     &           // ' from ' // TRIM(ADJ_FORCE_FILE)
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
            
         ELSE
            
            LENROW = NROWS 
            LENCOL = NCOLS   
         !WRITE(*,*) 'LENROW: ',LENROW, ', LENCOL: ',LENCOL
         !WRITE(*,*) 'BUFFER(1,1,1): ',BUFFER(1,1,1),' on proc: ',MYPE
         !WRITE(*,*) 'BUFFER(LC,LR,1): ',BUFFER(LENCOL,LENROW,NLAYS_FRC),' on proc: ',MYPE
            DO R = 1, LENROW
               DO C = 1, LENCOL
                  DO L = 1, NLAYS_FRC
                     !RARRAY = R + STRTROW - 1
                     !CARRAY = C + STRTCOL - 1
                     !ARRAY(CARRAY, RARRAY, L, V) = BUFFER(C, R, L)
                     ARRAY(C, R, L, VMAP_FRC(V)) = BUFFER(C, R, L)
                  END DO
               END DO
            END DO
      
      !write(*,*) 'sum(array): ',sum(array)
      !write(*,*) 'sum(array(:,:,:,4)): ',sum(array(:,:,:,4))
      
         END IF
      END DO

      END SUBROUTINE RD_FORCE_FILE

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      END MODULE ADJOINT_FILES
