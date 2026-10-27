
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE ADJ_BC_SENS

C-----------------------------------------------------------------------
C Function:
C   Sensitivity of the cost function to the lateral boundary condition,
C   written to ADJ_LGRID_BC with the same structure as BNDY_GASC_1: one
C   variable per BNDY_GASC_1 variable, NTHIK = 1, the standard S, E, N, W
C   perimeter, and BNDY_GASC_1's own time axis over the records the run
C   reads. Units are CF/ppmV.
C
C   Horizontal advection is the only process that reads the boundary
C   condition. In the forward model a domain-edge PE puts the coupled
C   boundary value BCFAC * BBUF * RHOJ_BND into the ghost cell on an
C   inflow face, where BBUF is BNDY_GASC_1 interpolated in time to the
C   mid-point of the advection sub-step. XADV_BWD and YADV_BWD return the
C   adjoint of that ghost cell; BC_SENS_ACCUM scales it by
C   BCFAC * RHOJ_BND (from RDBCON) and spreads it over the two bracketing
C   BNDY_GASC_1 records with INTERPB's weights.
C
C   A BNDY_GASC_1 record is read by the sub-steps on both sides of it, so
C   the sensitivity is held for the whole run and written once at the end
C   by BC_SENS_WRITE, rather than one output step at a time.
C
C   The four corner perimeter slots are never read by advection, so their
C   sensitivity is identically zero.
C
C   This module lives with the driver, not with the horizontal advection
C   that feeds it, so that the driver builds whichever hadv variant is
C   selected. Only hadv/zenodo_dadj_bwd calls BC_SENS_ACCUM; with any other
C   variant ADJ_LGRID_BC is written as zeros.
C
C   Nothing is computed or written unless ADJ_LGRID_BC is set.
C
C Revision History:
C   Oct 2026: created
C-----------------------------------------------------------------------

      USE ADJOINT_FILES, ONLY : ADJ_LGRID_BC

      IMPLICIT NONE
      SAVE

      LOGICAL :: BC_SENS_ON = .FALSE.

      INTEGER :: BC_TSTEP       ! BNDY_GASC_1 time step (HHMMSS), 0 if
                                ! the file is time-independent
      INTEGER :: BC_SDATE       ! date and time of the first output record,
      INTEGER :: BC_STIME       ! which is a BNDY_GASC_1 record
      INTEGER :: BC_NREC        ! number of output records
      INTEGER :: BC_NVARS       ! number of BNDY_GASC_1 variables
      INTEGER :: GL_NBNDY_BC    ! global perimeter length (NTHIK = 1)

      INTEGER, ALLOCATABLE :: BC_VAR(:) ! advected species -> BNDY_GASC_1
                                        ! variable, 0 if not read from it
      INTEGER, ALLOCATABLE :: BND_MAP(:) ! local perimeter slot -> global
                                         ! perimeter slot, 0 if the slot is
                                         ! not on the domain boundary

      REAL, ALLOCATABLE :: BC_SENS(:, :, :, :) ! (GL_NBNDY_BC, NLAYS,
                                               !  BC_NVARS, BC_NREC)

      CONTAINS

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Decide whether ADJ_LGRID_BC is wanted and set up the accumulation for
C the run SDATE:STIME to EDATE:ETIME. Called on every PE.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BC_SENS_INIT(SDATE, STIME, EDATE, ETIME)

      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Parameters:
      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV
     &     + N_AE_ADV
     &     + N_NR_ADV
     &     + N_TR_ADV
     &     + 1

C Arguments:
      INTEGER, INTENT(IN) :: SDATE, STIME ! run start
      INTEGER, INTENT(IN) :: EDATE, ETIME ! run end

C External Functions:
      INTEGER, EXTERNAL :: SECSDIFF, TIME2SEC, SEC2TIME, INDEX1

C Local Variables:
      INTEGER :: LOGDEV

      CHARACTER(16) :: PNAME = 'BC_SENS_INIT'
      CHARACTER(96) :: XMSG = ' '
      CHARACTER(80) :: VARDESC
      CHARACTER(256) :: FPATH

      CHARACTER(16) :: CONCMIN = '1.00E-30' ! as RDBCON, for the log
      CHARACTER(16) :: BCNAME(N_SPC_ADV)
      REAL :: BCFAC(N_SPC_ADV)

      INTEGER :: STATUS, ALLOCSTAT
      INTEGER :: DT, SECS, SECE, REC1, RECN
      INTEGER :: C, R, S

C-----------------------------------------------------------------------

      LOGDEV = INIT3()

      VARDESC = 'Cost function sensitivity to BNDY_GASC_1 output file'
      CALL ENVSTR(ADJ_LGRID_BC, VARDESC, ' ', FPATH, STATUS)
      IF ( STATUS /= 0 ) THEN
         WRITE(LOGDEV, '(/5X, A)') 'ADJ_LGRID_BC not set: the ' //
     &        'boundary condition sensitivity will not be written'
         BC_SENS_ON = .FALSE.
         RETURN
      END IF

C The same advected species -> boundary variable map RDBCON uses. Only gas
C species are read from BNDY_GASC_1, so only they contribute.
      CALL ADVBC_MAP(CONCMIN, BCNAME, BCFAC)

      IF ( .NOT. OPEN3(BNDY_GASC_1, FSREAD3, PNAME) ) THEN
         XMSG = 'Could not open ' // BNDY_GASC_1
         CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT1)
      END IF

      IF ( .NOT. DESC3(BNDY_GASC_1) ) THEN
         XMSG = 'Could not get ' // TRIM(BNDY_GASC_1) // ' file description'
         CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT1)
      END IF

      IF ( FTYPE3D /= BNDARY3 ) THEN
         XMSG = TRIM(BNDY_GASC_1) // ' is not a boundary file'
         CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT2)
      END IF

C XADV_BWD and YADV_BWD index the perimeter assuming one cell thickness
      IF ( NTHIK3D /= 1 ) THEN
         WRITE(XMSG, '(A, A, I3)') TRIM(BNDY_GASC_1),
     &        ' must have NTHIK = 1 for ADJ_LGRID_BC, it has ', NTHIK3D
         CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT2)
      END IF

      IF ( NCOLS3D /= GL_NCOLS .OR. NROWS3D /= GL_NROWS .OR.
     &     NLAYS3D /= NLAYS ) THEN
         XMSG = TRIM(BNDY_GASC_1) // ' does not match the model grid'
         CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT2)
      END IF

      BC_NVARS = NVARS3D
      GL_NBNDY_BC = 2 * (GL_NCOLS + GL_NROWS) + 4

      ALLOCATE(BC_VAR(N_SPC_ADV), BND_MAP(NBNDY), STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating BC_VAR, BND_MAP'
         CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT1)
      END IF

      BC_VAR = 0
      DO S = 1, N_GC_ADV
         IF ( BCNAME(S) /= ' ' )
     &        BC_VAR(S) = INDEX1(BCNAME(S), NVARS3D, VNAME3D)
      END DO

C Records the run reads. Every advection sub-step interpolates to its own
C mid-point, which is strictly inside the run, so the first is the record
C at or before the run start and the last the record at or after the run
C end.
      IF ( TSTEP3D == 0 ) THEN
         BC_TSTEP = 0
         BC_SDATE = SDATE3D
         BC_STIME = STIME3D
         BC_NREC = 1
      ELSE
         BC_TSTEP = TSTEP3D
         DT = TIME2SEC(TSTEP3D)
         SECS = SECSDIFF(SDATE3D, STIME3D, SDATE, STIME)
         SECE = SECSDIFF(SDATE3D, STIME3D, EDATE, ETIME)
         IF ( SECS < 0 .OR. SECE > (MXREC3D - 1) * DT ) THEN
            XMSG = TRIM(BNDY_GASC_1) // ' does not cover the run'
            CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT2)
         END IF
         REC1 = SECS / DT + 1
         RECN = (SECE + DT - 1) / DT + 1
         BC_NREC = RECN - REC1 + 1
         BC_SDATE = SDATE3D
         BC_STIME = STIME3D
         CALL NEXTIME(BC_SDATE, BC_STIME, SEC2TIME((REC1 - 1) * DT))
      END IF

C Map this PE's perimeter onto the global one. The local slots are laid
C out as pario's BOUNDARY fills them from the file: S, E, N, W, each side
C with its corner, so on a side that lies on the domain boundary the local
C slot of column C (or row R) is the global slot of that column (row).
      BND_MAP = 0
      IF ( ROWSX_PE(1, MYPE+1) == 1 ) THEN ! south
         DO C = 1, NCOLS
            BND_MAP(C) = COLSX_PE(1, MYPE+1) - 1 + C
         END DO
      END IF
      IF ( COLSX_PE(2, MYPE+1) == GL_NCOLS ) THEN ! east
         DO R = 1, NROWS
            BND_MAP(NCOLS + 1 + R) =
     &           GL_NCOLS + 1 + ROWSX_PE(1, MYPE+1) - 1 + R
         END DO
      END IF
      IF ( ROWSX_PE(2, MYPE+1) == GL_NROWS ) THEN ! north
         DO C = 1, NCOLS
            BND_MAP(NCOLS + NROWS + 3 + C) =
     &           GL_NCOLS + GL_NROWS + 3 + COLSX_PE(1, MYPE+1) - 1 + C
         END DO
      END IF
      IF ( COLSX_PE(1, MYPE+1) == 1 ) THEN ! west
         DO R = 1, NROWS
            BND_MAP(2 * NCOLS + NROWS + 4 + R) =
     &           2 * GL_NCOLS + GL_NROWS + 4 + ROWSX_PE(1, MYPE+1) - 1 + R
         END DO
      END IF

      ALLOCATE(BC_SENS(GL_NBNDY_BC, NLAYS, BC_NVARS, BC_NREC),
     &     STAT = ALLOCSTAT)
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating BC_SENS'
         CALL M3EXIT(PNAME, SDATE, STIME, XMSG, XSTAT1)
      END IF
      BC_SENS = 0.0

      BC_SENS_ON = .TRUE.

      WRITE(LOGDEV, '(/5X, A, I4, A, I8, ":", I6.6, A, I6.6)')
     &     'ADJ_LGRID_BC will hold', BC_NREC,
     &     ' records from', BC_SDATE, BC_STIME, ' every ', BC_TSTEP

      END SUBROUTINE BC_SENS_INIT

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Add one advection sub-step's boundary adjoint into BC_SENS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BC_SENS_ACCUM(BNDY_ADJ, BCPL, ADATE, ATIME, ASTEP_L,
     &     LAYER)

      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

C Arguments:
      REAL, INTENT(IN) :: BNDY_ADJ(:, :) ! adjoint of the coupled ghost
                                         ! values, (NBNDY, N_SPC_ADV)
      REAL, INTENT(IN) :: BCPL(:, :) ! d ghost / d BNDY_GASC_1 value, from
                                     ! RDBCON, (NBNDY, N_SPC_ADV)
      INTEGER, INTENT(IN) :: ADATE, ATIME ! start of the sub-step, as
                                          ! passed to RDBCON
      INTEGER, INTENT(IN) :: ASTEP_L ! advection sub-step (HHMMSS)
      INTEGER, INTENT(IN) :: LAYER

C External Functions:
      INTEGER, EXTERNAL :: SECSDIFF, TIME2SEC, SEC2TIME

C Local Variables:
      CHARACTER(16) :: PNAME = 'BC_SENS_ACCUM'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: MDATE, MTIME   ! mid-point of the sub-step
      INTEGER :: DT, SECS
      INTEGER :: IREC           ! earlier bracketing record, in BC_SENS
      INTEGER :: B, G, S, V
      REAL :: P, Q              ! weights of records IREC and IREC+1
      REAL :: X

C-----------------------------------------------------------------------

C RDBCON interpolates BNDY_GASC_1 to the mid-point of the sub-step
      MDATE = ADATE
      MTIME = ATIME
      CALL NEXTIME(MDATE, MTIME, SEC2TIME(TIME2SEC(ASTEP_L) / 2))

C The same linear weights as INTERPB (pario's INTERPOL)
      IF ( BC_TSTEP == 0 ) THEN
         IREC = 1
         P = 1.0
         Q = 0.0
      ELSE
         DT = TIME2SEC(BC_TSTEP)
         SECS = SECSDIFF(BC_SDATE, BC_STIME, MDATE, MTIME)
         IF ( SECS < 0 ) THEN
            IREC = 0
         ELSE
            IREC = SECS / DT + 1
         END IF
         Q = REAL(MOD(SECS, DT)) / REAL(DT)
         P = 1.0 - Q
         IF ( IREC < 1 .OR. IREC > BC_NREC .OR.
     &        ( Q > 0.0 .AND. IREC + 1 > BC_NREC ) ) THEN
            XMSG = 'Sub-step reads a boundary record outside ADJ_LGRID_BC'
            CALL M3EXIT(PNAME, MDATE, MTIME, XMSG, XSTAT3)
         END IF
      END IF

      DO S = 1, N_GC_ADV
         V = BC_VAR(S)
         IF ( V == 0 ) CYCLE
         DO B = 1, NBNDY
            G = BND_MAP(B)
            IF ( G == 0 ) CYCLE
            X = BNDY_ADJ(B, S) * BCPL(B, S)
            BC_SENS(G, LAYER, V, IREC) = BC_SENS(G, LAYER, V, IREC)
     &           + P * X
            IF ( Q > 0.0 )
     &           BC_SENS(G, LAYER, V, IREC+1) =
     &           BC_SENS(G, LAYER, V, IREC+1) + Q * X
         END DO
      END DO

      END SUBROUTINE BC_SENS_ACCUM

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Gather BC_SENS onto PE 0 and write it to ADJ_LGRID_BC. Called on every
C PE once the backward run is complete. The file itself is created by
C OP_ADJ_FILE.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BC_SENS_WRITE

      USE GRID_CONF             ! horizontal & vertical domain specifications

      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
#ifdef parallel
      INCLUDE SUBST_MPICH       ! MPI definitions and parameters
#endif

C External Functions:
C plain I/O API WRITE3: pario's PWRITE3 does not handle boundary files
      LOGICAL, EXTERNAL :: WRITE3_BNDY

C Local Variables:
      INTEGER :: LOGDEV

      CHARACTER(16) :: PNAME = 'BC_SENS_WRITE'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: WDATE, WTIME
      INTEGER :: N

#ifdef parallel
      REAL, ALLOCATABLE :: GBUF(:, :, :, :)
      INTEGER :: IERR, ALLOCSTAT
#endif

C-----------------------------------------------------------------------

      IF ( .NOT. BC_SENS_ON ) RETURN

      LOGDEV = INIT3()

#ifdef parallel
C Each PE has filled only the global slots on its own part of the domain
C boundary and left the rest zero, so summing over PEs gathers the pieces
      IF ( MYPE == 0 ) THEN
         ALLOCATE(GBUF(GL_NBNDY_BC, NLAYS, BC_NVARS, BC_NREC),
     &        STAT = ALLOCSTAT)
      ELSE
         ALLOCATE(GBUF(1, 1, 1, 1), STAT = ALLOCSTAT)
      END IF
      IF ( ALLOCSTAT /= 0 ) THEN
         XMSG = 'Failure allocating GBUF'
         CALL M3EXIT(PNAME, BC_SDATE, BC_STIME, XMSG, XSTAT1)
      END IF

      CALL MPI_REDUCE(BC_SENS, GBUF, SIZE(BC_SENS), MPI_REAL, MPI_SUM,
     &     0, MPI_COMM_WORLD, IERR)
      IF ( IERR /= MPI_SUCCESS ) THEN
         XMSG = 'MPI_REDUCE of BC_SENS failed'
         CALL M3EXIT(PNAME, BC_SDATE, BC_STIME, XMSG, XSTAT1)
      END IF

      IF ( MYPE == 0 ) BC_SENS = GBUF
      DEALLOCATE(GBUF)
#endif

      IF ( MYPE == 0 ) THEN
         WDATE = BC_SDATE
         WTIME = BC_STIME
         DO N = 1, BC_NREC
            IF ( .NOT. WRITE3_BNDY(ADJ_LGRID_BC, ALLVAR3, WDATE, WTIME,
     &           BC_SENS(:, :, :, N)) ) THEN
               XMSG = 'Could not write to ' // TRIM(ADJ_LGRID_BC)
               CALL M3EXIT(PNAME, WDATE, WTIME, XMSG, XSTAT1)
            END IF
            WRITE(LOGDEV, '(5X, 3(A, :, 1X), I8, ":", I6.6)')
     &           'Timestep written to', ADJ_LGRID_BC,
     &           'for date and time', WDATE, WTIME
            IF ( BC_TSTEP /= 0 ) CALL NEXTIME(WDATE, WTIME, BC_TSTEP)
         END DO

         IF ( .NOT. CLOSE3(ADJ_LGRID_BC) ) THEN
            XMSG = 'Could not close ' // TRIM(ADJ_LGRID_BC)
            CALL M3WARN(PNAME, WDATE, WTIME, XMSG)
         END IF
      END IF

      END SUBROUTINE BC_SENS_WRITE

      END MODULE ADJ_BC_SENS
