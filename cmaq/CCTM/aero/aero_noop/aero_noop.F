
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero_noop/aero_noop.F,v 1.1.1.1 2010/06/14 16:03:00 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)aero_driver_noop.F        1.1 /project/mod3/CMAQ/src/aero/aero_noop/SCCS/s.aero_driver_noop.F 13 Jun 1997 12:25:54

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AERO ( CGRID, JDATE, JTIME, TSTEP )

C   31 Mar 01 J.Young: generic CGRID declaration

      IMPLICIT NONE
 
C Includes:
 
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

c Arguments

C CGRID is incoming conc field (both gas and aerosol variables)
       
!     REAL      :: CGRID( :,:,:,: ) !  concentrations
      REAL, POINTER :: CGRID( :,:,:,: ) !  concentrations
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

C Local variables:
 
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

C-----------------------------------------------------------------------
 
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
         WRITE( LOGDEV,* ) '    '
         WRITE( LOGDEV,* ) '    Using no-op aerosol process'
         WRITE( LOGDEV,* ) '    '
         END IF
 
      RETURN
      END
