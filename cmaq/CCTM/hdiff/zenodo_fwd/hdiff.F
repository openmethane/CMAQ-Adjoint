
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hdiff/multiscale_adj_fwd/hdiff.F,v 1.1 2010/11/17 21:06:38 sergey Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HDIFF(CGRID, JDATE, JTIME, TSTEP)

C-----------------------------------------------------------------------
C Function:
C   Horizontal diffusion with constant eddy diffusivity - gen. coord.
C   The process time step is set equal to TSTEP(2). Boundary concentrations
C   are set using a no flux condition.

C Revision history:
C   Jul 2010 by Peter Percell at UH-IMAQS: revised for adjoint/4dvar
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE SUBST_MODULES         ! stenex

#ifdef parallel
      USE SWAP_SANDIA
#endif

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_GC_DIFF     ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF     ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF     ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF     ! tracer diffusion species and map table
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

C Parameters:
      INTEGER, PARAMETER :: N_SPC_DIFF = N_GC_DIFF ! advected species dimension
     &     + N_AE_DIFF
     &     + N_NR_DIFF
     &     + N_TR_DIFF

C Arguments:
      REAL, POINTER :: CGRID(:, :, :, :)
      INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

C External Functions not declared in IODECL3.EXT:
      INTEGER, EXTERNAL :: TIME2SEC

C File Variables:
      REAL :: RHOJ(0:NCOLS+1, 0:NROWS+1, NLAYS) ! density X Jacobian
      REAL :: K11BAR3D(NCOLS+1, NROWS+1, NLAYS) ! ave. Cx11 eddy diff
      REAL :: K22BAR3D(NCOLS+1, NROWS+1, NLAYS) ! ave. Cx22 eddy diff

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      CHARACTER(16) :: PNAME = 'HDIFF'
      CHARACTER(96) :: XMSG = ' '

      REAL :: DX1               ! dx1 (meters)
      REAL :: DX2               ! dx2 (meters)
      REAL, SAVE :: RDX1S       ! reciprocal dx1*dx1
      REAL, SAVE :: RDX2S       ! reciprocal dx2*dx2

      INTEGER, SAVE :: DIFF_MAP(N_SPC_DIFF) ! global diff map to CGRID

      LOGICAL, SAVE :: S_BND_PE, N_BND_PE, W_BND_PE, E_BND_PE

      REAL :: DTSEC             ! model time step in seconds
      REAL :: DT                ! diffusion time step
      REAL :: DTDX1S            ! dt/dx1**2
      REAL :: DTDX2S            ! dt/dx2**2

      INTEGER :: NSTEPS         ! number of diffusion time steps

      INTEGER :: C, R, L, S, V, N ! loop counters
      INTEGER :: D2C

      REAL :: RK11(NCOLS+1, NROWS+1, NLAYS) ! initially used as RHOJ
                                ! at x1 cell face, then finally as 11 eddy diff. factor
      REAL :: RK22(NCOLS+1, NROWS+1, NLAYS) ! initially used as RHOJ
                                ! at x2 cell face, then finally as 22 eddy diff. factor

!slz      REAL, TARGET :: CONC(0:NCOLS+1, 0:NROWS+1, NLAYS, N_SPC_DIFF) ! conc working array
      REAL, ALLOCATABLE, TARGET :: CONC(:,:,:,:)   !slz
      INTEGER :: ALLOCSTAT

#ifdef parallel
      REAL, POINTER :: SEND(:, :, :, :)
      REAL, POINTER :: RECV(:, :, :, :)
#endif

C Required interfaces
      INTERFACE
         SUBROUTINE RHO_J(JDATE, JTIME, TSTEP, RHOJ)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: JDATE, JTIME, TSTEP(2)
            REAL, INTENT(OUT)   :: RHOJ(:, :, :)
         END SUBROUTINE RHO_J
         SUBROUTINE HCDIFF3D(JDATE, JTIME, K11BAR, K22BAR, DT)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: JDATE, JTIME
            REAL, INTENT(OUT)   :: K11BAR(:, :, :), K22BAR(:, :, :)
            REAL, INTENT(OUT)   :: DT
         END SUBROUTINE HCDIFF3D
      END INTERFACE

C-----------------------------------------------------------------------
      ALLOCATE(CONC(0:NCOLS+1, 0:NROWS+1, NLAYS, N_SPC_DIFF), STAT=ALLOCSTAT)  !SLZ
      IF ( ALLOCSTAT /=0 ) THEN  !SLZ
         XMSG = 'FAILURE ALLOCATING NADVSS'  !SLZ
         CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)   !SLZ
      END IF  !SLZ

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

C Get dx1 and dx2
         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX1 = DG2M * XCELL_GD
     &           * COS(PI180 * (YORIG_GD + YCELL_GD * REAL(GL_NROWS / 2))) ! in m.
            DX2 = DG2M * YCELL_GD ! in m.
         ELSE
            DX1 = XCELL_GD      ! in m.
            DX2 = YCELL_GD      ! in m.
         END IF

         RDX1S = 1.0 / (DX1 * DX1)
         RDX2S = 1.0 / (DX2 * DX2)

C Get CGRID offsets
         CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

C Create global map to CGRID
         S = 0
         DO V = 1, N_GC_DIFF
            S = S + 1
            DIFF_MAP(S) = GC_STRT - 1 + GC_DIFF_MAP(V)
         END DO
         DO V = 1, N_AE_DIFF
            S = S + 1
            DIFF_MAP(S) = AE_STRT - 1 + AE_DIFF_MAP(V)
         END DO
         DO V = 1, N_NR_DIFF
            S = S + 1
            DIFF_MAP(S) = NR_STRT - 1 + NR_DIFF_MAP(V)
         END DO
         DO V = 1, N_TR_DIFF
            S = S + 1
            DIFF_MAP(S) = TR_STRT - 1 + TR_DIFF_MAP(V)
         END DO

         CALL SUBST_HI_LO_BND_PE('C', W_BND_PE, E_BND_PE)
         CALL SUBST_HI_LO_BND_PE('R', S_BND_PE, N_BND_PE)
      END IF                    ! FIRSTIME

C Get the computational grid (rho X Jacobian) for this step
      CALL RHO_J(JDATE, JTIME, TSTEP, RHOJ)

#ifdef parallel
C swap3d interface ...
C     swap3d( 1st value of array to be sent,
C             1st value of array to be received,
C             number of values to send in 1st dimension of array,
C             number of values to send in 2nd dimension of array,
C             number of values to send in 3rd dimension of array,
C             size of 1st dimension in arrays,
C             size of 2nd dimension in arrays,
C             direction to receive from: NSEW)

      CALL SWAP3D(RHOJ(1, NROWS, 1), RHOJ(1, 0, 1), NCOLS, 1, NLAYS,
     &     NCOLS+2, NROWS+2, SOUTH)
      CALL SWAP3D(RHOJ(NCOLS, 1, 1), RHOJ(0, 1, 1), 1, NROWS, NLAYS,
     &     NCOLS+2, NROWS+2, WEST)
#endif

C Get FLUX point values for RHOJ
      DO L = 1, NLAYS
         DO R = 1, NROWS + 1
            DO C = 1, NCOLS + 1
               RK11(C, R, L) = 0.5 * (RHOJ(C, R, L) + RHOJ(C-1, R, L))
               RK22(C, R, L) = 0.5 * (RHOJ(C, R, L) + RHOJ(C, R-1, L))
            END DO
         END DO
      END DO

C Get the contravariant eddy diffusivities
      CALL HCDIFF3D(JDATE, JTIME, K11BAR3D, K22BAR3D, DT)

C Get number of steps based on eddy time
      DTSEC = REAL(TIME2SEC(TSTEP(2)))
      NSTEPS = INT(DTSEC / DT) + 1
      DT = DTSEC / REAL(NSTEPS)

      WRITE(LOGDEV, '(/ 5X, A, 1PE15.7, I8)')
     &     'H-eddy DT & integration steps: ', DT, NSTEPS

      DTDX1S = DT * RDX1S
      DTDX2S = DT * RDX2S

C Determine FLUX point values of contravariant diffusivity coupled with RHOJ
      DO L = 1, NLAYS
         DO R = 1, NROWS + 1
            DO C = 1, NCOLS + 1
               RK11(C, R, L) = RK11(C, R, L) * K11BAR3D(C, R, L)
               RK22(C, R, L) = RK22(C, R, L) * K22BAR3D(C, R, L)
            END DO
         END DO
      END DO

#ifdef parallel
      CALL SWAP3D(RK11(1, 1, 1), RK11(NCOLS+1, 1, 1), 1, NROWS, NLAYS,
     &     NCOLS+1, NROWS+1, EAST)
      CALL SWAP3D(RK22(1, 1, 1), RK22(1, NROWS+1, 1), NCOLS, 1, NLAYS,
     &     NCOLS+1, NROWS+1, NORTH)
#endif

C Enforce no flux boundary conditions by setting RK11 and RK22 to zero at flux
C points on cell edges along the external boundary of the global domain; make
C sure boundary values are initialized with a valid real number (value doesn't
C matter)
      IF ( W_BND_PE ) THEN      ! West boundary
         RK11(1, :, :) = 0.0
         CONC(0, :, :, :) = 0.0
      END IF
      IF ( E_BND_PE ) THEN      ! East boundary
         RK11(NCOLS+1, :, :) = 0.0
         CONC(NCOLS+1, :, :, :) = 0.0
      END IF
      IF ( S_BND_PE ) THEN      ! South boundary
         RK22(:, 1, :) = 0.0
         CONC(:, 0, :, :) = 0.0
      END IF
      IF ( N_BND_PE ) THEN      ! North boundary
         RK22(:, NROWS+1, :) = 0.0
         CONC(:, NROWS+1, :, :) = 0.0
      END IF

C NOTE: There has been a long-standing bug in this code that is fixed here.
C       The loop setting CONC from CGRID must be inside the loop on NSTEPS
C       because that loop updates CGRID, not CONC.

      DO N = 1, NSTEPS

C Load working array (CGRID is coupled, CONC is mixing ratio)
         DO S = 1, N_SPC_DIFF
            D2C = DIFF_MAP(S)
            CONC(1:NCOLS, 1:NROWS, 1:NLAYS, S) =
     &           CGRID(1:NCOLS, 1:NROWS, 1:NLAYS, D2C) /
     &           RHOJ(1:NCOLS, 1:NROWS, 1:NLAYS)
         END DO

C Get values of CONC in columns 0 & NCOLS+1 and rows 0 and NROWS+1 from
C other processors, as needed

#ifdef parallel
C swap4d interface ...
C     swap4d( 1st value of array to be sent,
C             1st value of array to be received,
C             number of values to send in 1st dimension of array,
C             number of values to send in 2nd dimension of array,
C             number of values to send in 3rd dimension of array,
C             number of values to send in 4th dimension of array,
C             size of 1st dimension in arrays,
C             size of 2nd dimension in arrays,
C             size of 3rd dimension in arrays,
C             direction to receive from: NSEW)

         SEND => CONC(1:NCOLS, 1:1, :, :)
         RECV => CONC(1:NCOLS, NROWS+1:NROWS+1, :, :)
         CALL SWAP4D(SEND, RECV, NCOLS, 1, NLAYS, N_SPC_DIFF,
     &        NCOLS, 1, NLAYS, NORTH)

         SEND => CONC(1:1, 1:NROWS, :, :)
         RECV => CONC(NCOLS+1:NCOLS+1, 1:NROWS, :, :)
         CALL SWAP4D(SEND, RECV, 1, NROWS, NLAYS, N_SPC_DIFF,
     &        1, NROWS, NLAYS, EAST)

         SEND => CONC(1:NCOLS, NROWS:NROWS, :, :)
         RECV => CONC(1:NCOLS, 0:0, :, :)
         CALL SWAP4D(SEND, RECV, NCOLS, 1, NLAYS, N_SPC_DIFF,
     &        NCOLS, 1, NLAYS, SOUTH)

         SEND => CONC(NCOLS:NCOLS, 1:NROWS, :, :)
         RECV => CONC(0:0, 1:NROWS, :, :)
         CALL SWAP4D(SEND, RECV, 1, NROWS, NLAYS, N_SPC_DIFF,
     &        1, NROWS, NLAYS, WEST)
#endif

C Do the gridded computation for horizontal diffusion
         DO S = 1, N_SPC_DIFF
            D2C = DIFF_MAP(S)

            DO L = 1, NLAYS
               DO R = 1, NROWS
                  DO C = 1, NCOLS
                     CGRID(C, R, L, D2C) = RHOJ(C, R, L) * CONC(C, R, L, S)
     &                    + DTDX1S
     &                    * ( RK11(C+1, R, L)
     &                    * ( CONC(C+1, R, L, S) - CONC(C, R, L, S) )
     &                    - RK11(C, R, L)
     &                    * ( CONC(C, R, L, S) - CONC(C-1, R, L, S) ) )
     &                    + DTDX2S
     &                    * ( RK22(C, R+1, L)
     &                    * ( CONC(C, R+1, L, S) - CONC(C, R, L, S) )
     &                    - RK22(C, R, L)
     &                    * ( CONC(C, R, L, S) - CONC(C, R-1, L, S) ) )
                  END DO
               END DO
            END DO
         END DO

      END DO                    ! NSTEPS

      DEALLOCATE(CONC)  !slz

      END SUBROUTINE HDIFF
