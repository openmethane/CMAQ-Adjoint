      MODULE ADJ_FDM_TEST

C-----------------------------------------------------------------------
C Function:
C   Evaluate the ADJOINT against the FDM on a point-to-point base
C   cell ghost region around CGRID

C Preconditions:
C   GRID_CONF
C   CGRID_SPCS

C Subroutines and functions called:

C Revision history:
C   May 07 - slz - initial implementation

C-----------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL :: FDM_TEST
      CHARACTER*16 :: PTB_SPC_NAME
      INTEGER :: PTB_SPC
      INTEGER :: PTB_TIME
      REAL    :: PTB
      LOGICAL :: PTB_ABS
      INTEGER :: PTB_COL1, PTB_COL2
      INTEGER :: PTB_ROW1, PTB_ROW2
      INTEGER :: PTB_LAY1, PTB_LAY2

      CONTAINS
         SUBROUTINE FDM_INIT (JDATE, JTIME)
         USE GRID_CONF           ! horizontal & vertical domain specifications
         USE CGRID_SPCS          ! CGRID species number and offsets

         INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
         INCLUDE SUBST_IODECL      ! I/O definitions and declarations
         INCLUDE SUBST_GC_SPC      ! gas chemistry species table
         INCLUDE SUBST_AE_SPC      ! aerosol species table
         INCLUDE SUBST_NR_SPC      ! non-reactive species table
         INCLUDE SUBST_TR_SPC      ! tracer species table

         INTEGER :: JDATE, JTIME
         INTEGER :: ALLOCSTAT
         INTEGER, SAVE :: LOGDEV
         CHARACTER*16  :: PNAME='FDM_INIT'
         CHARACTER( 120 ) :: XMSG = ' '

         INTEGER :: NVARS
         CHARACTER(16) :: VAL_LIST(3)
         INTEGER :: STATUS         ! ENVINT or ENVREAL status
         CHARACTER(80) :: VARDESC  ! environment variable description

         INTEGER :: COL_OFF, ROW_OFF ! offsets w.r.t. global domain

         LOGICAL, EXTERNAL :: ENVYN
         INTEGER, EXTERNAL :: ENVINT, INDEX1
         REAL, EXTERNAL    :: ENVREAL

         INTERFACE
            SUBROUTINE GET_ENVLIST ( ENV_VAR, NVARS, VAL_LIST )
               IMPLICIT NONE
               CHARACTER( * ),  INTENT ( IN )  :: ENV_VAR
               INTEGER,         INTENT ( OUT ) :: NVARS
               CHARACTER( 16 ), INTENT ( OUT ) :: VAL_LIST( : )
            END SUBROUTINE GET_ENVLIST
         END INTERFACE

         LOGDEV = INIT3()
         VARDESC='Retrieve FDM test flag'
         FDM_TEST = .false. !default value for no FDM test
         FDM_TEST = envyn('FDM_TEST', VARDESC , FDM_TEST, STATUS)
         IF ( STATUS > 0 ) THEN
            XMSG = 'Invalid value for FDM_TEST'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

         IF(.NOT.FDM_TEST) THEN
            WRITE (LOGDEV,*) "The current run is not for FDM test"
            RETURN; END IF

!adopted from define_receptor.F
! Get CGRID offsets
         CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

! Get name of perturbation species
         VARDESC = 'retrieve perturbation species name'
         CALL ENVSTR('PTB_SPC_NAME', VARDESC, PNAME, PTB_SPC_NAME, STATUS)
         IF ( STATUS /= 0 ) THEN
            XMSG = 'Invalid value for PTB_SPC_NAME'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

! Get perturbation type
         VARDESC='Retrieve PTB_ABS flag'
         PTB_ABS = .false. !by default, relative perturbation
         PTB_ABS = envyn('PTB_ABS', VARDESC , PTB_ABS, STATUS)
         IF ( STATUS == 1 ) THEN
            XMSG = 'Invalid value for PTB_ABS'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

! Get pertubation size
         PTB = 0          ! default
         VARDESC = 'retrieve perturbation size'
         PTB = ENVREAL('PTB_SIZE', VARDESC, PTB, STATUS)
         IF ( STATUS == 1 ) THEN
            XMSG = 'Invalid value for PTB_SIZE'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

! Get pertubation time
         PTB_TIME = 0          ! default
         VARDESC = 'retrieve perturbation time'
         PTB_TIME = ENVINT('PTB_TIME', VARDESC, PTB_TIME, STATUS)
         IF ( STATUS == 1 ) THEN
            XMSG = 'Invalid value for PTB_TIME'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

! Look for index of receptor's species in lists of known species
         PTB_SPC = INDEX1(PTB_SPC_NAME, N_GC_SPC, GC_SPC)
         IF ( PTB_SPC == 0 ) THEN
            PTB_SPC = INDEX1(PTB_SPC_NAME, N_AE_SPC, AE_SPC)
            IF ( PTB_SPC == 0 ) THEN
               PTB_SPC = INDEX1(PTB_SPC_NAME, N_NR_SPC, NR_SPC)
               IF ( PTB_SPC == 0 ) THEN
                  PTB_SPC = INDEX1(PTB_SPC_NAME, N_TR_SPC, TR_SPC)
                  IF ( PTB_SPC == 0 ) THEN
                     !slz IF (trim(ptb_spc_name).eq.'fine-pm') Then
                     !slz ptb_spc=-1 !specail cases use negative integers
               !slz to be done
                     XMSG = 'Perturbation species name not recognized: ' //
     &                      TRIM(PTB_SPC_NAME)
                     CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
                  ELSE
                     PTB_SPC = TR_STRT - 1 + PTB_SPC
                  END IF
               ELSE
                  PTB_SPC = NR_STRT - 1 + PTB_SPC
               END IF
            ELSE
               PTB_SPC = AE_STRT - 1 + PTB_SPC
            END IF
         ELSE
            PTB_SPC = GC_STRT - 1 + PTB_SPC
         END IF
         WRITE (LOGDEV,*) "Perturbation species "//trim(PTB_SPC_NAME), PTB_SPC

! Get perturbation location as ranges of columns, rows and levels
         CALL GET_ENVLIST('PTB_CRL1', NVARS, VAL_LIST)
         IF ( NVARS.LT.2 ) THEN
            XMSG = 'Error -- PTB_CRL1 = "COL ROW" or "COL ROW LAY"'
            CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
         ELSE IF( NVARS.EQ.2)  THEN
            READ(VAL_LIST(1), '(I4)') PTB_COL1
            READ(VAL_LIST(2), '(I4)') PTB_ROW1
            PTB_LAY1=1
         ELSE
            READ(VAL_LIST(1), '(I4)') PTB_COL1
            READ(VAL_LIST(2), '(I4)') PTB_ROW1
            READ(VAL_LIST(3), '(I4)') PTB_LAY1
         END IF

         CALL GET_ENVLIST('PTB_CRL2', NVARS, VAL_LIST)
         IF ( NVARS.LT.2 ) THEN
            XMSG = 'Error -- PTB_CRL2 = "COL ROW" or "COL ROW LAY"'
            CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
         ELSE IF( NVARS.EQ.2)  THEN
            READ(VAL_LIST(1), '(I4)') PTB_COL2
            READ(VAL_LIST(2), '(I4)') PTB_ROW2
            PTB_LAY1=1
         ELSE
            READ(VAL_LIST(1), '(I4)') PTB_COL2
            READ(VAL_LIST(2), '(I4)') PTB_ROW2
            READ(VAL_LIST(3), '(I4)') PTB_LAY2
         END IF

! Check
         IF ( PTB_COL1 .LT. 1 .OR. PTB_COL2 .GT. GL_NCOLS .OR.
     &        PTB_ROW1 .LT. 1 .OR. PTB_ROW2 .GT. GL_NROWS .OR.
     &        PTB_LAY1 .LT. 1 .OR. PTB_LAY2 .GT. NLAYS ) THEN
            XMSG = 'Part of the perturbation area is out of the CTM domain'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT2)
         END IF

! Find local indices
         COL_OFF = COLSX_PE(1, MYPE+1) - 1
         PTB_COL1 = MAX(PTB_COL1 - COL_OFF, 1)
         PTB_COL2 = MIN(PTB_COL2 - COL_OFF, NCOLS)

         ROW_OFF = ROWSX_PE(1, MYPE+1) - 1
         PTB_ROW1 = MAX(PTB_ROW1 - ROW_OFF, 1)
         PTB_ROW2 = MIN(PTB_ROW2 - ROW_OFF, NROWS)

         IF ( PTB_COL2 >= PTB_COL1 .AND. PTB_ROW2 >= PTB_ROW1 ) THEN
            WRITE (LOGDEV, *) "Perturbation to be applied to CPU #", MYPE
            WRITE (LOGDEV, *) "    to an area with bottom left corner of (", PTB_COL1,",",PTB_ROW1,")"
            WRITE (LOGDEV, *) "    and top right corner of (", PTB_COL2,",",PTB_ROW2,")"
            WRITE (LOGDEV, *) "    for layers from (", PTB_lay1,",",PTB_lay2,")"
         END IF

         END SUBROUTINE FDM_INIT

         SUBROUTINE FDM_PTB(CGRID, JDATE, JTIME)
         IMPLICIT NONE
         INCLUDE SUBST_IODECL      ! I/O definitions and declarations

         REAL, POINTER :: CGRID(:,:,:,:)
         INTEGER :: JDATE
         INTEGER :: JTIME
         INTEGER, SAVE :: LOGDEV

         IF ( PTB_COL2 >= PTB_COL1 .AND. PTB_ROW2 >= PTB_ROW1 ) THEN
           IF( PTB_ABS ) THEN
             CGRID(PTB_COL1:PTB_COL2, PTB_ROW1:PTB_ROW2, PTB_LAY1:PTB_LAY2, PTB_SPC)
     &       = PTB + CGRID(PTB_COL1:PTB_COL2, PTB_ROW1:PTB_ROW2, PTB_LAY1:PTB_LAY2, PTB_SPC)
           ELSE
             CGRID(PTB_COL1:PTB_COL2, PTB_ROW1:PTB_ROW2, PTB_LAY1:PTB_LAY2, PTB_SPC)
     &       = (1. + PTB) * CGRID(PTB_COL1:PTB_COL2, PTB_ROW1:PTB_ROW2, PTB_LAY1:PTB_LAY2, PTB_SPC)
           END IF
         END IF
         LOGDEV = INIT3()
         WRITE (LOGDEV,*) "Perturbation applied at", JDATE, JTIME
         END SUBROUTINE FDM_PTB

      END MODULE ADJ_FDM_TEST
