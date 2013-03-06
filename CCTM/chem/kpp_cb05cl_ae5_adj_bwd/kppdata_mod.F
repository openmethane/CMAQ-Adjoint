
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/chem/kpp_cb05cl_ae5_adj_bwd/kppdata_mod.F,v 1.1 2010/11/17 21:28:23 sergey Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%



      MODULE  HRDATA

C*************************************************************************
C
C  FUNCTION:  Mechanism & solver data for EBI solver
C             
C  PRECONDITIONS: For CB4 family of mechanisms only
C 
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, April, 2003
C                   C*************************************************************************

c..Mechanism specific variables
      INTEGER   N_SPEC       ! No. of species in mechanism
      INTEGER   N_RXNS       ! No. of reactions in mechanism

c..Control flags
      LOGICAL   L_AE_VRSN    ! Flag for aerosol version of mech
      LOGICAL   L_AQ_VRSN    ! Flag for aqueous chemistry version of mech

c..Miscellaneous variables
      INTEGER   LOGDEV       ! Unit number of output log

c...Allocatable arrays
      INTEGER, ALLOCATABLE, SAVE :: ICTRL( : )
      INTEGER, ALLOCATABLE :: ISTAT( : )

      REAL*8, ALLOCATABLE ::  RKI(   : )     ! Rate constants 
      REAL*8, ALLOCATABLE ::  YCD(   : )     ! concentrations
      REAL*8, ALLOCATABLE ::  LAMBDA( : )     ! adjoints
      REAL*8, ALLOCATABLE, SAVE ::  RCTRL( : )     ! Rate constants 
      REAL*8, ALLOCATABLE ::  RSTAT( : )     ! Rate constants 
      REAL, ALLOCATABLE ::  RXRAT( : )     ! Reaction rates 
      REAL*8, ALLOCATABLE, SAVE ::  RTOL(  : )     ! Species relative tolerances 
      REAL*8, ALLOCATABLE, SAVE ::  ATOL(  : )     ! Species absolute tolerances 
      REAL, ALLOCATABLE ::  YC(    : )     ! Species concentrations

      REAL*8 :: CHEMSTD
      INTEGER :: KPP_ERR

c..Gas species indices
      INTEGER   NO2
      INTEGER   NO
      INTEGER   O
      INTEGER   O3
      INTEGER   NO3
      INTEGER   O1D
      INTEGER   OH
      INTEGER   HO2
      INTEGER   N2O5
      INTEGER   HNO3
      INTEGER   HONO
      INTEGER   PNA
      INTEGER   H2O2
      INTEGER   XO2
      INTEGER   XO2N
      INTEGER   NTR
      INTEGER   ROOH
      INTEGER   FORM
      INTEGER   ALD2
      INTEGER   ALDX
      INTEGER   PAR
      INTEGER   CO
      INTEGER   MEO2
      INTEGER   MEPX
      INTEGER   MEOH
      INTEGER   HCO3
      INTEGER   FACD
      INTEGER   C2O3
      INTEGER   PAN
      INTEGER   PACD
      INTEGER   AACD
      INTEGER   CXO3
      INTEGER   PANX
      INTEGER   ROR
      INTEGER   OLE
      INTEGER   ETH
      INTEGER   IOLE
      INTEGER   TOL
      INTEGER   CRES
      INTEGER   TO2
      INTEGER   TOLRO2
      INTEGER   OPEN
      INTEGER   CRO
      INTEGER   MGLY
      INTEGER   XYL
      INTEGER   XYLRO2
      INTEGER   ISOP
      INTEGER   ISPD
      INTEGER   ISOPRXN
      INTEGER   TERP
      INTEGER   TRPRXN
      INTEGER   SO2
      INTEGER   SULF
      INTEGER   SULRXN
      INTEGER   ETOH
      INTEGER   ETHA
      INTEGER   CL2
      INTEGER   CL
      INTEGER   HOCL
      INTEGER   CLO
      INTEGER   FMCL
      INTEGER   HCL
      INTEGER   TOLNRXN
      INTEGER   TOLHRXN
      INTEGER   XYLNRXN
      INTEGER   XYLHRXN
      INTEGER   BENZENE
      INTEGER   BENZRO2
      INTEGER   BNZNRXN
      INTEGER   BNZHRXN
      INTEGER   SESQ
      INTEGER   SESQRXN

      END MODULE HRDATA
