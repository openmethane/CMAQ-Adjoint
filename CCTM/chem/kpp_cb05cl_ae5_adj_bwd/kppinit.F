
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/kpp_cb05cl_ae5_adj_bwd/kppinit.F,v 1.1 2010/11/17 21:28:23 sergey Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%



      SUBROUTINE HRINIT
C***********************************************************************
C
C  FUNCTION: To initialize species tolerances, arrays, and indices 
C
C  PRECONDITIONS: For CB4 mechanisms only
C
C  RETURN VALUES: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, January, 2002
C
C***********************************************************************

      USE HRDATA

      IMPLICIT NONE

C.....INCLUDES:
      INCLUDE SUBST_GC_SPC    ! Gas chem species names and MWs
      INCLUDE SUBST_GC_EMIS   ! Gas chem emissions name and mapping tables
      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block

C.....ARGUMENTS: NONE
 
C.....PARAMETERS: NONE


C.....EXTERNAL FUNCTIONS:
      INTEGER  FINDEX         ! Finds location of a number in a list

C.....SAVED VARIABLES: NONE

 
C.....LOCAL VARIABLES: 
      CHARACTER*16  PNAME     ! Program name
      CHARACTER*132 MSG       ! Log message

      INTEGER IND             ! Species index
      INTEGER N               ! Loop index

C***********************************************************************
      DATA PNAME / 'HRINIT' / 


      N_SPEC = N_GC_SPC
      N_RXNS = NRXNS

      ALLOCATE( ICTRL( 20 ) )
      ALLOCATE( ISTAT( 20 ) )
      ALLOCATE( RCTRL( 20 ) )
      ALLOCATE( RSTAT( 20 ) )
      ALLOCATE( RKI( NRXNS ) )
      ALLOCATE( YCD( N_SPEC ) )
      ALLOCATE( LAMBDA( N_SPEC ) )
      ALLOCATE( RXRAT( NRXNS ) )
      ALLOCATE( RTOL( N_SPEC) )
      ALLOCATE( ATOL( N_SPEC) )
      ALLOCATE( YC(   N_SPEC) )

      ICTRL(:)  = 0
      RCTRL(:)  = 0.0d0
      ICTRL(1)=1
      
      RETURN

      END
