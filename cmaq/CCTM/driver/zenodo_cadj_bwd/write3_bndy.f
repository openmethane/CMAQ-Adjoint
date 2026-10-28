
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      FUNCTION WRITE3_BNDY ( FNAME, VNAME, JDATE, JTIME, DATA )
     &                       RESULT ( SUCCESS )

C-----------------------------------------------------------------------
C Function: Use I/O-API's WRITE3, not pario's PWRITE3, to write the
C           ADJ_LGRID_BC boundary file from PE 0. PWRITE3 only handles
C           gridded files: for any other type it logs a warning and
C           returns failure, so the output would silently go missing.

C Note: This function must live in a file that has a ".f" extension, so
C       that the build's -DWRITE3=PWRITE3 is not applied to it (see
C       vdiff/zenodo_bwd/write3_distr.f)

C Revision History:
C   Oct 2026: created for ADJ_LGRID_BC
C-----------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER( * ), INTENT( IN ) :: FNAME
      CHARACTER( * ), INTENT( IN ) :: VNAME
      INTEGER,        INTENT( IN ) :: JDATE, JTIME
      REAL,           INTENT( IN ) :: DATA( * )
      LOGICAL SUCCESS

      LOGICAL, EXTERNAL :: WRITE3

C-----------------------------------------------------------------------

      SUCCESS = WRITE3( FNAME, VNAME, JDATE, JTIME, DATA )

      RETURN
      END FUNCTION WRITE3_BNDY
