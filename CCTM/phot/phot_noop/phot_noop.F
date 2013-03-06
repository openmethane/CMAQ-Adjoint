
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/phot/phot_noop/phot_noop.F,v 1.2 2011/08/11 09:58:06 jaroslav Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)phot_noop.F       1.1 /project/mod3/CMAQ/src/phot/phot_noop/SCCS/s.phot_noop.F 03 Jun 1997 11:47:13

      SUBROUTINE PHOT ( MDATE, MTIME, JDATE, JTIME, NDARK, RJ )

C**********************************************************************

C FUNCTION: 
C    PHOT, adapted from RADM, calculates the photolysis rate constants
C    to be used by the chemical solver.
C    It uses linear interpolation in time of day, height, and latitude
C    from file tabular values and optionally adjusts photolysis rates 
C    above, below and in cumulus clouds.

C PRECONDITIONS: None

C KEY SUBROUTINES/FUNCTIONS CALLED: None

C REVISION HISTORY: 
C    prototype(adaptation from RADM), Rohit Mathur, April 1993.  
C    major mods, Jeff Young, May 1994 - annotated and/or "c" in col 1
C   30 Mar 01 J.Young: Use HGRD_DEFN
C   12 Feb 05 J.Young: vertical dyn alloc - Use VGRD_DEFN

C**********************************************************************

      USE HGRD_DEFN           ! horizontal domain specifications
      USE VGRD_DEFN             ! vertical layer specifications

      IMPLICIT NONE

C...include files

!     INCLUDE SUBST_HGRD_ID      ! horizontal dimensioning params
!     INCLUDE SUBST_VGRD_ID      ! vertical dimensioning params
      INCLUDE SUBST_RXCMMN       ! chemical mechamism reactions COMMON
      INCLUDE SUBST_IODECL       ! I/O definitions and declarations

C...arguments

      INTEGER MDATE              ! "centered" Julian date (YYYYDDD)
      INTEGER MTIME              ! "centered" time (HHMMSS)
      INTEGER JDATE              ! current Julian date (YYYYDDD)
      INTEGER JTIME              ! current time (HHMMSS)
      INTEGER NDARK              ! Number of level 1 cells in darkness
      
      REAL    RJ( NCOLS, NROWS, NLAYS, NPHOTAB )
                                     ! gridded J-values  (/min units)
 
C...local variables

      LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! .TRUE. if 1st call of subroutine

      INTEGER, SAVE :: LOGDEV

      INTEGER ROW
      INTEGER COL
      INTEGER LEV
      INTEGER JP                 ! loop indices

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
         END IF

      NDARK = GL_NCOLS * GL_NROWS
      DO JP = 1, NPHOTAB
         DO LEV = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  RJ( COL,ROW,LEV,JP ) = 0.0
                  END DO
               END DO
           END DO
         END DO

      WRITE( LOGDEV, 1003) JDATE, JTIME
1003  FORMAT( 8X, 'No-op photolysis: at', I9.7, ':', I6.6,
     &        1X, 'GMT - zero photolytic rate constants')
      RETURN
      END
