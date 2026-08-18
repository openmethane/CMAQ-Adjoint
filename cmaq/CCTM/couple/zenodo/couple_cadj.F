
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/couple/gencoor/couple.F,v 1.3 2011/07/26 15:32:24 ppercell Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)couple.F	1.1 /project/mod3/CMAQ/src/convert/couple/SCCS/s.couple.F 03 Jun 1997 11:41:33

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE COUPLE_CADJ ( LGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   Convert units in LGRID for continuous adjoint of advection
 
C Preconditions:
 
C Subroutines and functions called:
 
C Revision History:
C    S. Zhao  May 2019 - adapted from couple.F
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE   

C Include files:

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

!     REAL      :: LGRID( :,:,:,: )  ! concentrations
      REAL, POINTER :: LGRID( :,:,:,: )  ! concentrations
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
C Parameters:

C File Variables:
      REAL        RHO  ( NCOLS,NROWS,NLAYS )  !air density

C External Functions (not already declared by IODECL3.EXT):

C Local Variables:

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: NQAE              ! number of micro-grams / m**3 species
      INTEGER, SAVE :: QAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""
      INTEGER, SAVE :: NNAE              ! number of  # / m**3 species
      INTEGER, SAVE :: NAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""
      INTEGER, SAVE :: NSAE              ! number of  m**2 / m**3 species
      INTEGER, SAVE :: SAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""

      INTEGER     GXOFF, GYOFF          ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

      INTEGER     NSPCS            ! number of species for subroutine calls
      INTEGER     OFF              ! loop offset to CGRID species
      INTEGER     C, R, L, SPC     ! loop counters

      INTEGER     LDATE, LTIME

      REAL, ALLOCATABLE, SAVE :: DS ( : ) ! dx3 (dimensionless in sigma coord.)
      INTEGER       ALLOCSTAT
      CHARACTER( 16 ) :: PNAME = 'DECOUPLE_CADJ'
      CHARACTER( 16 ) :: VNAME
      CHARACTER( 96 ) :: XMSG = ' '

      REAL, PARAMETER :: GPKG = 1.0E+03              ! g/kg
      REAL, PARAMETER :: MGPG = 1.0E+06              ! micro-g/g
      REAL            :: CONV

C-----------------------------------------------------------------------

      VNAME = 'DENS'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME,
     &                    RHO ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      NSPCS = N_GC_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = GC_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,OFF+SPC ) = LGRID( C,R,L,OFF+SPC )
     &                                      * DS (L) * RHO ( C,R,L)
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
	CONV = 1. / (GPKG * MGPG)
      NSPCS = NQAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,QAE( SPC ) ) = LGRID( C,R,L,QAE( SPC ) )
     &                                         * DS( L ) * CONV
                     END DO
                  END DO
               END DO
            END DO
         END IF

      NSPCS = NNAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,NAE( SPC ) ) = LGRID( C,R,L,NAE( SPC ) )
     &                                         * DS( L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
      NSPCS = NSAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,SAE( SPC ) ) = LGRID( C,R,L,SAE( SPC ) )
     &                                         * DS( L )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      NSPCS = N_NR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = NR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,OFF+SPC ) = LGRID( C,R,L,OFF+SPC )
     &                                      * DS( L )  * RHO( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
      NSPCS = N_TR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = TR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,OFF+SPC ) = LGRID( C,R,L,OFF+SPC )
     &                                      * DS( L ) * RHO( C,R,L )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      RETURN
 
C............................................................................
C entry DECOUPLE
 
      ENTRY DECOUPLE_CADJ ( LGRID, JDATE, JTIME, TSTEP )

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         ALLOCATE ( DS( NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DS'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Get dx3 from COORD include file
         DO L = 1, NLAYS
            DS ( L ) = ABS ( X3FACE_GD( L ) - X3FACE_GD( L-1 ) )
         END DO

C Get CGRID offsets
 
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C create aerosol species pointers to distinguish micro-grams / m**3,
C # / m**3 (number density), and m**2 / m**3 (surface area) species
 
         NQAE = 0       ! no. of micro-grams / m**3 species
         NNAE = 0       ! no. of  # / m**3 species
         NSAE = 0       ! no. of  m**2 / m**3 species
         OFF = AE_STRT - 1
         DO SPC = 1, N_AE_SPC
            IF ( AE_SPC( SPC )( 1:3 ) .EQ. 'NUM' ) THEN
               NNAE = NNAE + 1
               NAE( NNAE ) = OFF + SPC
               ELSE IF ( AE_SPC( SPC )( 1:3 ) .EQ. 'SRF' ) THEN
               NSAE = NSAE + 1
               SAE( NSAE ) = OFF + SPC
               ELSE
               NQAE = NQAE + 1
               QAE( NQAE ) = OFF + SPC
               END IF
            END DO

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

C For backward runs of CMAQ/ADJOINT, the function INTERPX needs to be told
C to use a backward interpolation interval for a first time call at the last
C date & time in the MET_CRO_3D file
         CALL LSTEPF( MET_CRO_3D, LDATE, LTIME )
         IF ( JDATE .EQ. LDATE .AND. JTIME .EQ. LTIME ) THEN
            CALL NEXTIME(LDATE, LTIME, -TSTEP(2))

            VNAME = 'DENS'
            IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                          LDATE, LTIME,
     &                          RHO ) ) THEN
               XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
               CALL M3EXIT( PNAME, LDATE, LTIME, XMSG, XSTAT1 )
               END IF
            END IF

         END IF       ! if firstime

      VNAME = 'DENS'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME,
     &                    RHO ) ) THEN
         XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      NSPCS = N_GC_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = GC_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,OFF+SPC ) = LGRID( C,R,L,OFF+SPC )
     &                                      / ( RHO( C,R,L) * DS( L ))
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
      NSPCS = NQAE
	CONV = 1. / (GPKG * MGPG)
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,QAE( SPC ) ) = LGRID( C,R,L,QAE( SPC ) )
     &                                         / ( CONV * DS( L ))
                     END DO
                  END DO
               END DO
            END DO
         END IF

      NSPCS = NNAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,NAE( SPC ) ) = LGRID( C,R,L,NAE( SPC ) )
     &                                         / DS( L )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
      NSPCS = NSAE
      IF ( NSPCS .GT. 0 ) THEN
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,SAE( SPC ) ) = LGRID( C,R,L,SAE( SPC ) )
     &                                         / DS( L )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      NSPCS = N_NR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = NR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,OFF+SPC ) = LGRID( C,R,L,OFF+SPC )
     &                                      / ( RHO( C,R,L ) * DS( L ) )
                     END DO
                  END DO
               END DO
            END DO
         END IF
 
      NSPCS = N_TR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = TR_STRT - 1
         DO SPC = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     LGRID( C,R,L,OFF+SPC ) = LGRID( C,R,L,OFF+SPC )
     &                                      / ( RHO( C,R,L) * DS( L ) )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      RETURN
      END
