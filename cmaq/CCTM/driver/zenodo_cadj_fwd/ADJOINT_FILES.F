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
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

C-----------------------------------
C Data
C-----------------------------------

C Parameters
      CHARACTER(16), PARAMETER :: CTM_XFIRST_IN   = 'CTM_XFIRST_IN'
      CHARACTER(16), PARAMETER :: CTM_XFIRST_OUT  = 'CTM_XFIRST_OUT'

      CHARACTER(16), PARAMETER :: CTM_CONC_FWD    = 'CTM_CONC_FWD'
      CHARACTER(16), PARAMETER :: INIT_LGRID_1    = 'INIT_LGRID_1'
      CHARACTER(16), PARAMETER :: ADJ_LGRID       = 'ADJ_LGRID'

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      CONTAINS

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Create new adjoint file
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE OP_ADJ_FILE(JDATE, JTIME, TSTEP, ADJ_FILE)

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
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
      CHARACTER(16) :: PNAME = 'OP_ADJ_FILE'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: K, S           ! loop counters

C-----------------------------------------------------------------------

      IF ( MYPE == 0 ) THEN

         IF ( ADJ_FILE == ADJ_LGRID ) THEN

C Get file characteristics of CTM_CONC_FWD
            IF ( .NOT. OPEN3(CTM_CONC_FWD, FSREAD3, PNAME) ) THEN
               XMSG = 'Could not open ' // TRIM(CTM_CONC_FWD) //
     &              ' file as OLD:READ-ONLY'
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
            END IF

            IF ( .NOT. DESC3(CTM_CONC_FWD) ) THEN
               XMSG = 'Could not get ' // TRIM(CTM_CONC_FWD) // ' file description'
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
            END IF

            IF ( .NOT. CLOSE3(CTM_CONC_FWD) ) THEN
               XMSG = 'Could not close ' // TRIM(CTM_CONC_FWD)
               CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
            END IF

         ELSE

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

         END IF

C Set file characteristics that differ from the base file, CTM_CONC_FWD or
C CTM_CONC_1
         SELECT CASE ( ADJ_FILE )

         CASE ( ADJ_LGRID )

            IF ( CONST_SYNC_STEP ) THEN
               TSTEP3D = TSTEP(2) ! sync step
            ELSE
               TSTEP3D = TSTEP(1) ! output step
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

            DO S = 1, N_NR_SPC
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = NR_SPC(S)
               UNITS3D(NVARS3D) = 'CF/ppmV'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to ' //
     &              TRIM(VNAME3D(NVARS3D))
            END DO

            DO S = 1, N_TR_SPC
               NVARS3D = NVARS3D + 1
               VTYPE3D(NVARS3D) = M3REAL
               VNAME3D(NVARS3D) = TR_SPC(S)
               UNITS3D(NVARS3D) = 'CF/ppmV'
               VDESC3D(NVARS3D) = 'Sensitivity of cost function to ' //
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
#include      SUBST_IODECL      ! I/O definitions and declarations

C Arguments:
      REAL, INTENT(IN) :: ARRAY(*)
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

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_IOFDESC     ! file header data structure

C Arguments:
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
      DO V = 1, NVARS3D
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      END MODULE ADJOINT_FILES
