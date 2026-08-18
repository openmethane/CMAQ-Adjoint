
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

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE DEFINE_RECEPTOR(LGRID)

C-----------------------------------------------------------------------
C Function:
C   Define the receptor, i.e., initialize LGRID

C Revision History:
C   Apr 2007 Kumaresh Singh and Adrian Sandu at Viginia Tech: created
C   Aug 2007 Peter Percell at UH-IMAQS: read enviroment variables
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions (XSTAT*)

C Arguments:
      REAL, POINTER :: LGRID(:, :, :, :)

C External Functions:
      INTEGER, EXTERNAL :: INDEX1

C Local Variables:
      CHARACTER(16) :: PNAME = 'DEFINE_RECEPTOR'
      CHARACTER(96) :: XMSG = ' '

      INTEGER :: STATUS

      CHARACTER(16) :: SPC_NAME

      INTEGER :: NVARS
      CHARACTER(16) :: VAL_LIST(2)

      INTEGER :: SPC
      INTEGER :: BCOL, ECOL, BROW, EROW, BLEV, ELEV, NCELLS_GL
      INTEGER :: COL_OFF, ROW_OFF ! offsets w.r.t. global domain

      INTERFACE
         SUBROUTINE GET_ENVLIST(ENV_VAR, NVARS, VAL_LIST)
            IMPLICIT NONE
            CHARACTER(*),  INTENT(IN)  :: ENV_VAR
            INTEGER,       INTENT(OUT) :: NVARS
            CHARACTER(16), INTENT(OUT) :: VAL_LIST(:)
         END SUBROUTINE GET_ENVLIST
      END INTERFACE

C-----------------------------------------------------------------------

C Get CGRID offsets
      CALL CGRID_MAP(NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT)

C Get name of receptor's species
      CALL ENVSTR('RECPTR_SPC', 'Receptor species name', ' ', SPC_NAME, STATUS)
      IF ( STATUS /= 0 ) THEN
         XMSG = 'Invalid input for receptor species name'
         CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
      END IF

C Look for index of receptor's species in lists of known species
      SPC = INDEX1(SPC_NAME, N_GC_SPC, GC_SPC)
      IF ( SPC == 0 ) THEN
         SPC = INDEX1(SPC_NAME, N_AE_SPC, AE_SPC)
         IF ( SPC == 0 ) THEN
            SPC = INDEX1(SPC_NAME, N_NR_SPC, NR_SPC)
            IF ( SPC == 0 ) THEN
               SPC = INDEX1(SPC_NAME, N_TR_SPC, TR_SPC)
               IF ( SPC == 0 ) THEN
                  IF (trim(SPC_NAME).eq.'PMIJ') THEN
                     SPC = -1
                  ELSE
                     XMSG = 'Receptor species name not recognized: ' //
     &                    TRIM(SPC_NAME)
                     CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
                  END IF
               ELSE
                  SPC = TR_STRT - 1 + SPC
               END IF
            ELSE
               SPC = NR_STRT - 1 + SPC
            END IF
         ELSE
            SPC = AE_STRT - 1 + SPC
         END IF
      ELSE
         SPC = GC_STRT - 1 + SPC
      END IF

C Get receptor's location as ranges of columns, rows and levels
      CALL GET_ENVLIST('RECPTR_BCOL_ECOL', NVARS, VAL_LIST)
      IF ( NVARS /= 2 ) THEN
         XMSG = 'Error in column range for receptor'
         CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
      END IF
      READ(VAL_LIST(1), '(I4)') BCOL
      READ(VAL_LIST(2), '(I4)') ECOL

      CALL GET_ENVLIST('RECPTR_BROW_EROW', NVARS, VAL_LIST)
      IF ( NVARS /= 2 ) THEN
         XMSG = 'Error in row range for receptor'
         CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
      END IF
      READ(VAL_LIST(1), '(I4)') BROW
      READ(VAL_LIST(2), '(I4)') EROW

      CALL GET_ENVLIST('RECPTR_BLEV_ELEV', NVARS, VAL_LIST)
      IF ( NVARS /= 2 ) THEN
         XMSG = 'Error in level range for receptor'
         CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
      END IF
      READ(VAL_LIST(1), '(I4)') BLEV
      READ(VAL_LIST(2), '(I4)') ELEV

      IF ( BCOL .LT. 1 .OR. ECOL .GT. GL_NCOLS .OR.
     &     BROW .LT. 1 .OR. EROW .GT. GL_NROWS .OR.
     &     BLEV .LT. 1 .OR. ELEV .GT. NLAYS ) THEN
         XMSG = 'Receptor domain inconsistent with CTM domain'
         CALL M3EXIT(PNAME, 0, 0, XMSG, XSTAT2)
      END IF

C Number of receptor cells in global domain
      NCELLS_GL = (ECOL-BCOL+1) * (EROW-BROW+1) * (ELEV-BLEV+1)

      COL_OFF = COLSX_PE(1, MYPE+1) - 1
      BCOL = MAX(BCOL - COL_OFF, 1)
      ECOL = MIN(ECOL - COL_OFF, NCOLS)

      ROW_OFF = ROWSX_PE(1, MYPE+1) - 1
      BROW = MAX(BROW - ROW_OFF, 1)
      EROW = MIN(EROW - ROW_OFF, NROWS)

C Initialize adjoint variable, LGRID

      IF ( ECOL >= BCOL .AND. EROW >= BROW ) THEN
         IF (SPC.GT.0) THEN
            LGRID(BCOL:ECOL, BROW:EROW, BLEV:ELEV, SPC) = 1.0D0 / DBLE(NCELLS_GL)
         ELSE IF (SPC.EQ.-1) THEN
!slz PM IJ; summation of dry i- and j-mode aerosols
            LGRID(BCOL:ECOL, BROW:EROW, BLEV:ELEV, 74:101) = 1.0D0 / DBLE(NCELLS_GL)
            LGRID(BCOL:ECOL, BROW:EROW, BLEV:ELEV, 112:115) = 1.0D0 / DBLE(NCELLS_GL)
            LGRID(BCOL:ECOL, BROW:EROW, BLEV:ELEV, 122:124) = 1.0D0 / DBLE(NCELLS_GL)
         END IF
      END IF

!slz summation of the 72 gas species
!slzchem      LGRID(:, :, :, 1) = 1.
!slzchem      LGRID(:, :, :, 2) = 1.
!slzchem      LGRID(:, :, :, 3) = 1.
!slzchem      LGRID(:, :, :, 4) = 1.
!slzchem      LGRID(:, :, :, 5) = 1.
!slzchem      LGRID(:, :, :, 6) = 1.
!slzchem      LGRID(:, :, :, 7) = 1.
!slzchem      LGRID(:, :, :, 8) = 1.
!slzchem      LGRID(:, :, :, 9) = 1.
!slzchem      LGRID(:, :, :, 10) = 1.
!slzchem      LGRID(:, :, :, 11) = 1.
!slzchem      LGRID(:, :, :, 12) = 1.
!slzchem      LGRID(:, :, :, 13) = 1.
!slzchem      LGRID(:, :, :, 14) = 1.
!slzchem      LGRID(:, :, :, 15) = 1.
!slzchem      LGRID(:, :, :, 16) = 1.
!slzchem      LGRID(:, :, :, 17) = 1.
!slzchem      LGRID(:, :, :, 18) = 1.
!slzchem      LGRID(:, :, :, 19) = 1.
!slzchem      LGRID(:, :, :, 20) = 1.
!slzchem      LGRID(:, :, :, 21) = 1.
!slzchem      LGRID(:, :, :, 22) = 1.
!slzchem      LGRID(:, :, :, 23) = 1.
!slzchem      LGRID(:, :, :, 24) = 1.
!slzchem      LGRID(:, :, :, 25) = 1.
!slzchem      LGRID(:, :, :, 26) = 1.
!slzchem      LGRID(:, :, :, 27) = 1.
!slzchem      LGRID(:, :, :, 28) = 1.
!slzchem      LGRID(:, :, :, 29) = 1.
!slzchem      LGRID(:, :, :, 30) = 1.
!slzchem      LGRID(:, :, :, 31) = 1.
!slzchem      LGRID(:, :, :, 32) = 1.
!slzchem      LGRID(:, :, :, 33) = 1.
!slzchem      LGRID(:, :, :, 34) = 1.
!slzchem      LGRID(:, :, :, 35) = 1.
!slzchem      LGRID(:, :, :, 36) = 1.
!slzchem      LGRID(:, :, :, 37) = 1.
!slzchem      LGRID(:, :, :, 38) = 1.
!slzchem      LGRID(:, :, :, 39) = 1.
!slzchem      LGRID(:, :, :, 40) = 1.
!slzchem      LGRID(:, :, :, 41) = 1.
!slzchem      LGRID(:, :, :, 42) = 1.
!slzchem      LGRID(:, :, :, 43) = 1.
!slzchem      LGRID(:, :, :, 44) = 1.
!slzchem      LGRID(:, :, :, 45) = 1.
!slzchem      LGRID(:, :, :, 46) = 1.
!slzchem      LGRID(:, :, :, 47) = 1.
!slzchem      LGRID(:, :, :, 48) = 1.
!slzchem      LGRID(:, :, :, 49) = 1.
!slzchem      LGRID(:, :, :, 50) = 1.
!slzchem      LGRID(:, :, :, 51) = 1.
!slzchem      LGRID(:, :, :, 52) = 1.
!slzchem      LGRID(:, :, :, 53) = 1.
!slzchem      LGRID(:, :, :, 54) = 1.
!slzchem      LGRID(:, :, :, 55) = 1.
!slzchem      LGRID(:, :, :, 56) = 1.
!slzchem      LGRID(:, :, :, 57) = 1.
!slzchem      LGRID(:, :, :, 58) = 1.
!slzchem      LGRID(:, :, :, 59) = 1.
!slzchem      LGRID(:, :, :, 60) = 1.
!slzchem      LGRID(:, :, :, 61) = 1.
!slzchem      LGRID(:, :, :, 62) = 1.
!slzchem      LGRID(:, :, :, 63) = 1.
!slzchem      LGRID(:, :, :, 64) = 1.
!slzchem      LGRID(:, :, :, 65) = 1.
!slzchem      LGRID(:, :, :, 66) = 1.
!slzchem      LGRID(:, :, :, 67) = 1.
!slzchem      LGRID(:, :, :, 68) = 1.
!slzchem      LGRID(:, :, :, 69) = 1.
!slzchem      LGRID(:, :, :, 70) = 1.
!slzchem      LGRID(:, :, :, 71) = 1.
!slzchem      LGRID(:, :, :, 72) = 1.

!slz o3
!slzo3      LGRID(:, :, :, 4) = 1.
      RETURN
      END
