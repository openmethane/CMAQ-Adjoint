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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/EMIS_DEFN.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE EMIS_DEFN

C-----------------------------------------------------------------------
C Function: emissions interface to the chemistry-transport model

C Revision History:
C     28 Jul 2006 J.Young: initial implementation
C     18 Aug 2007 J.Young: move beis part to separate module; add plume rise
C     26 Jan 2010 J.Young: fix bug overwriting point source layer 1 NH3; inline rdemis

C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications

      IMPLICIT NONE

      REAL, ALLOCATABLE, SAVE :: VDEMIS( :,:,:,: ) ! total emissions array

      INTEGER, SAVE :: N_SPC_EMIS ! global emissions species
      INTEGER, PUBLIC, SAVE :: EMLAYS     ! no. of emis integration layers
                                  ! EMLAYS must be .LT. NLAYS
      INTEGER,          SAVE :: EMLYRS     ! no. of layers on base emis input file -
                                           ! a 2d file if inline, otherwise a 3d file
      INTEGER, SAVE :: NEMIS_AE   ! no. of aero emis species

      LOGICAL, SAVE :: EM_TRAC = .FALSE. ! do tracer emissions?
      REAL, SAVE :: CNVTE      ! combined conversion factor
      REAL,    PRIVATE, SAVE :: DT         ! TSTEP (output) in sec

      REAL,    PRIVATE, ALLOCATABLE, SAVE :: VDEMIS_AE( :,:,:,: ) ! aerosol emis
      REAL,    PRIVATE, ALLOCATABLE, SAVE :: BUFF( :,:,: )        ! read buffer

      INTEGER, SAVE :: STRT_GC, FINI_GC, STRT_AE, FINI_AE,
     &                          STRT_NR, FINI_NR, STRT_TR, FINI_TR
      INTEGER, PRIVATE, SAVE :: SDATE, STIME  ! scenario start date/time (beis)
      INTEGER, PRIVATE       :: LDATE, LTIME  ! step start date/time (beis)
      INTEGER, PRIVATE, SAVE :: NDATE, NTIME  ! step next date/time (beis)

      CONTAINS

C-----------------------------------------------------------------------
         FUNCTION EMIS_INIT ( JDATE, JTIME, TSTEP ) RESULT ( SUCCESS ) 

!slz         USE BEIS_DEFN           ! biogenic emissions
!slz         USE PT3D_DEFN           ! plume rise emissions
         USE PMEM_DEFN           ! defines emissions file(s) PM mapping

         IMPLICIT NONE
      
C Includes:

         INCLUDE SUBST_AE_SPC    ! aerosol species table
         INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names & map table
         INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names & map table
         INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names & map table
         INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names & map table
         INCLUDE SUBST_CONST     ! constants
         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include         SUBST_IODECL    ! I/O definitions and declarations

C Arguments:

         INTEGER JDATE, JTIME, TSTEP
         LOGICAL :: SUCCESS

C External Functions:

         INTEGER, EXTERNAL :: ENVINT

C Parameters:

         REAL, PARAMETER :: CMLMR = 1.0E+06       ! ppmV/Molar Mixing Ratio

C Local Variables:

         REAL    CONVEM          ! conversion for emissions rates to Kg/s
         REAL    DX1, DX2        ! CX x1- and x2-cell widths

         CHARACTER( 16 ), SAVE :: CTM_EMLAYS = 'CTM_EMLAYS'    ! env var for no. of
                                                               ! emiss. layers
         CHARACTER( 16 ) :: PNAME = 'EMIS_INIT'
         CHARACTER( 80 ) :: VARDESC   ! env variable description
         CHARACTER( 120 ) :: XMSG = ' '
         INTEGER V, K, LOGDEV, STATUS

C-----------------------------------------------------------------------

         SUCCESS = .TRUE.

         LOGDEV = INIT3()

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX1 = DG2M * XCELL_GD ! in m.
            DX2 = DG2M * YCELL_GD
     &          * COS( PI180*( YORIG_GD + YCELL_GD*FLOAT( GL_NROWS/2 )))! in m.
         ELSE
            DX1 = XCELL_GD        ! in m.
            DX2 = YCELL_GD        ! in m.
         END IF

C In-line biogenic emissions?

!slz         IF ( .NOT. BEIS_INIT( JDATE, JTIME, TSTEP ) ) THEN
!slz            XMSG = 'Failure initializing biogenics emissions module'
!slz            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!slz         END IF

C Open Emissions files

         CALL OPEMIS ( JDATE, JTIME, N_GC_EMIS, EM_TRAC, CONVEM, EMLAYS )

C Get number of emissions layers

         IF ( EMLAYS .LE. 1 ) THEN   ! better be inline plume rise run
            EMLAYS = NLAYS   ! env var default
            VARDESC = 'Number of emission layers'
            EMLAYS = ENVINT( CTM_EMLAYS, VARDESC, EMLAYS, STATUS )
            IF ( STATUS .NE. 0 ) THEN
               WRITE( LOGDEV, '(5X, A)' ) VARDESC
               IF ( STATUS .EQ. 1 ) THEN
                  XMSG = 'Environment variable improperly formatted'
                  CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
                  SUCCESS = .FALSE.; RETURN
               ELSE
                  XMSG = 'Environment variable not set or empty ... Using default:'
                  WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
               END IF
            END IF
         END IF

         EMLAYS = MIN ( EMLAYS, NLAYS )

         WRITE( LOGDEV,1009 ) EMLAYS, NLAYS
1009     FORMAT( / 10X, 'Number of Emissions Layers:         ', I3
     &           / 10X, 'out of total Number of Model Layers:', I3 )

         IF ( N_AE_SPC .GT. 0 ) THEN
            WRITE( LOGDEV,'( /5X, A )' ) 'Aerosol Emissions Processing in '
     &                                // 'Vertical diffusion ...'
            NEMIS_AE = N_AE_EMIS  ! from AE_EMIS.EXT
C set up PM species mapping from emissions file(s) for use in GET_AERO_EMIS
            CALL GET_PMEM_MAP ( JDATE, JTIME )   ! set in PMEM_DEFN
         ELSE
            NEMIS_AE = 0
         END IF

         N_SPC_EMIS = N_GC_EMIS
     &              + NEMIS_AE
     &              + N_NR_EMIS
     &              + N_TR_EMIS

!slzC 3D point source emissions (N_SPC_EMIS+1 serves as the null map for PTEM_MAP)
!slz         IF ( .NOT. PT3D_INIT( N_SPC_EMIS, EMLAYS, JDATE, JTIME, TSTEP ) ) THEN
!slz            XMSG = 'Failure initializing plume rise emissions module'
!slz            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!slz         END IF
         
!         WRITE(*,*) 'N_SPC_EMIS: ',N_SPC_EMIS,', EMLAYS: ',EMLAYS 
!         WRITE(*,*) 'MY_NCOLS: ',MY_NCOLS,', MY_NROWS: ',MY_NROWS
!         WRITE(*,*) 'VDEMIS: ',VDEMIS

         ALLOCATE ( VDEMIS ( N_SPC_EMIS+1,EMLAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = 'VDEMIS memory allocation failed'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
            SUCCESS = .FALSE.; RETURN
         END IF

!slz         IF ( PT3DEMIS ) THEN
!slz            EMLYRS = 1
!slz         ELSE
            EMLYRS = EMLAYS
!slz         END IF

         IF ( N_SPC_EMIS .LE. 0 ) THEN
            XMSG = 'There are no emissions: VDEMIS is set to zero' ! below
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
         END IF

         IF ( N_AE_SPC .GT. 0 ) THEN   ! all EMLAYS
            ALLOCATE ( VDEMIS_AE( NEMIS_AE,EMLAYS,MY_NCOLS,MY_NROWS ),
     &                 STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = 'VDEMIS_AE memory allocation failed'
               CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
               SUCCESS = .FALSE.; RETURN
            END IF
            VDEMIS_AE = 0.0   ! array assignment
         END IF

         ALLOCATE ( BUFF( NCOLS,NROWS,EMLYRS ), STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = 'BUFF memory allocation failed'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
            SUCCESS = .FALSE.; RETURN
         END IF
         BUFF = 0.0   ! array assignment

C combined gas emssions conversion factor

         CNVTE = CMLMR * CONVEM * MWAIR / ( DX1 * DX2 )

         STRT_GC = 1
         FINI_GC = N_GC_EMIS
         STRT_AE = N_GC_EMIS + 1
         FINI_AE = N_GC_EMIS + NEMIS_AE
         STRT_NR = N_GC_EMIS + NEMIS_AE + 1
         FINI_NR = N_GC_EMIS + NEMIS_AE + N_NR_EMIS
         STRT_TR = N_GC_EMIS + NEMIS_AE + N_NR_EMIS + 1
         FINI_TR = N_GC_EMIS + NEMIS_AE + N_NR_EMIS + N_TR_EMIS

         SUCCESS = .TRUE.; RETURN

         END FUNCTION EMIS_INIT
      END MODULE EMIS_DEFN
