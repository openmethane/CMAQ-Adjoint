      MODULE EMIS_DEFN_B

!slz adapted from emis_defn by shunliu
      USE EMIS_DEFN
      IMPLICIT NONE

      CONTAINS

C-----------------------------------------------------------------------
         SUBROUTINE GET_EMIS_B ( JDATE, JTIME, CONVPA,
     &                         RDX3F, RRHOJ, RJACM, VDEMISB, cgrid_em_bwd)

         USE AERO_EMIS_b           ! inherits GRID_CONF
         USE PMEM_DEFN           ! defines emissions file(s) PM mapping

         IMPLICIT NONE

C Includes:

         INCLUDE SUBST_AE_SPC    ! aerosol species table
         INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names & map table
         INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names & map table
         INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names & map table
         INCLUDE SUBST_FILES_ID  ! file name parameters
         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
         INCLUDE SUBST_IOFDESC   ! file header data structure
         INCLUDE SUBST_IODECL    ! I/O definitions and declarations

C Arguments:

         INTEGER JDATE, JTIME
         INTEGER TSTEP( 2 )                 ! time step vector (HHMMSS)
         REAL, POINTER :: CGRID_EM_BWD( :,:,:,: ) ! sensitivities of J w.r.t. emissions
         REAL :: vdemisb(:,:,:,:)
         REAL    CONVPA             ! conversion factor to pressure in Pascals
         REAL    RDX3F( NLAYS )             ! reciprocal layer thickness
         REAL    RRHOJ( NCOLS,NROWS,NLAYS ) ! reciprocal density X Jacobian
         REAL    RJACM( NCOLS,NROWS,NLAYS ) ! reciprocal mid-layer Jacobian

C External Functions:

         INTEGER, EXTERNAL :: TRIMLEN

C Local Variables:

         REAL             CNVTC         ! combined conversion factor
         REAL             CNVTI         ! intermediate combined conv. factor
         REAL,    SAVE :: CNVTP         ! intermediate combined conv. factor
         REAL             DELT          ! interpolation factor
         INTEGER          C, R, L, N, S, V ! loop induction variables
         INTEGER          isens !slz counter for sensitivities to emissions

         CHARACTER( 16 ) :: PNAME = 'GET_EMIS_B'
         CHARACTER( 120 ) :: XMSG = ' '
         LOGICAL, SAVE :: FIRSTIME = .TRUE.
         REAL, ALLOCATABLE, SAVE :: VDEMIS_AEB( :,:,:,: )
         INTEGER STATUS

C-----------------------------------------------------------------------

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            CNVTP = CONVPA * CNVTE

            IF ( N_AE_SPC .GT. 0 ) THEN   ! all EMLAYS
               ALLOCATE ( VDEMIS_AEB( NEMIS_AE,EMLAYS,MY_NCOLS,MY_NROWS ),
     &                     STAT = STATUS )
               IF ( STATUS .NE. 0 ) THEN
                  XMSG = 'VDEMIS_AEB memory allocation failed'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
            END IF
         END IF

         isens=0

         IF ( N_AE_SPC .GT. 0 ) THEN
            vdemis_aeb = 0.0
            DO r=1,my_nrows
              DO c=1,my_ncols
                DO l=1,emlays
                  DO v=strt_ae,fini_ae
                    s = v + 1 - strt_ae
                    vdemis_aeb(s, l, c, r) = vdemis_aeb(s, l, c, r) + vdemisb(v, l, c, r)
                  END DO
                END DO
              END DO
            END DO
            CALL GET_AERO_EMIS_B ( JDATE, JTIME, EMLAYS, RJACM,
     &                             VDEMISb, VDEMIS_AEb, cgrid_em_bwd, isens )
!slz vdemisb is updated: h2so4 vapor emission added to fine-PM sulfate emissions
         END IF

         DO V = STRT_GC, FINI_GC
            S = V + 1 - STRT_GC
            CNVTI = CNVTP * GC_EMIS_FAC( S )
            isens = isens + 1
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  DO L = 1, EMLYRS
                     CNVTC = CNVTI * RDX3F( L ) * RRHOJ( C,R,L )
                     VDEMISB( V,L,C,R ) = VDEMISB( V,L,C,R ) * CNVTC

                     cgrid_em_bwd(c,r,l,isens) = vdemisb(v,l,c,r)
                     vdemisb(v,l,c,r) = 0.0
                  END DO
               END DO
            END DO
         END DO   ! V


C non-reactive gases (conversion to ppmv/s) VDEMIS in this order from RDEMIS

         DO V = STRT_NR, FINI_NR
            S = V + 1 - STRT_NR
            CNVTI = CNVTP * NR_EMIS_FAC( S )
            isens = isens + 1
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  DO L = 1, EMLYRS
                     CNVTC = CNVTI * RDX3F( L ) * RRHOJ( C,R,L )
                     VDEMISB( V,L,C,R ) = VDEMISB( V,L,C,R ) * CNVTC

                     cgrid_em_bwd(c,r,l,isens) = vdemisb(v,l,c,r)
                     vdemisb(v,l,c,r) = 0.0
                  END DO
               END DO
            END DO
         END DO

C tracer gases (conversion to ppmv/s)

         IF ( EM_TRAC ) THEN
            DO V = STRT_TR, FINI_TR
               S = V + 1 - STRT_TR
               CNVTI = CNVTP * TR_EMIS_FAC( S )
               isens = isens + 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     DO L = 1, EMLYRS
                        CNVTC = CNVTI * RDX3F( L ) * RRHOJ( C,R,L )
                        VDEMISB( V,L,C,R ) = CNVTC * VDEMISB( V,L,C,R )

                        cgrid_em_bwd(c,r,l,isens) = vdemisb(v,l,c,r)
                        vdemisb(v,l,c,r) = 0.0
                     END DO
                  END DO
               END DO
            END DO
         END IF

        END SUBROUTINE GET_EMIS_B

      END MODULE EMIS_DEFN_B
