C=======================================================================
C
C *** ISORROPIA CODE
C *** INCLUDE FILE 'ISRPIA_B.INC'
C *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
C     AND VARIABLES. 
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER nslds
      INTEGER nzsr
      INTEGER npair
      INTEGER ncomp
      INTEGER nions
      INTEGER ngasaq
      INTEGER nerrmx
C
C *** END OF INCLUDE FILE **********************************************
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** INCLUDE FILE 'ISRPIA.INC'
C *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
C     AND VARIABLES. 
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      PARAMETER (ncomp=5, nions=7, ngasaq=3, nslds=9, npair=13, nzsr=100
     +   , nerrmx=25)
C
C *** INPUT VARIABLES **************************************************
C
      INTEGER METSTBL, IPROB, NADJ
      DOUBLE PRECISION W, WAER, TEMP, RH
      COMMON /INPT/ W(NCOMP), WAER(NCOMP), TEMP, RH, IPROB, METSTBL,
     &              NADJ
C
C *** WATER ACTIVITIES OF PURE SALT SOLUTIONS **************************
C
      DOUBLE PRECISION AWAS, AWSS, AWAC, AWSC, AWAN, AWSN, AWSB, AWAB
      DOUBLE PRECISION AWSA, AWLC
      COMMON /ZSR / AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR),
     &              AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR),
     &              AWSA(NZSR), AWLC(NZSR)
C
C *** DELIQUESCENCE RELATIVE HUMIDITIES ********************************
C
      INTEGER WFTYP
      DOUBLE PRECISION DRH2SO4, DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3
      DOUBLE PRECISION DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
C
      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3, 
     &              DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
C
      DOUBLE PRECISION DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2
      DOUBLE PRECISION DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2
      DOUBLE PRECISION DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3
      DOUBLE PRECISION DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8
      DOUBLE PRECISION DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13
C
      COMMON /MDRH/ DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2,
     &              DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2,
     &              DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3,
     &              DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8,
     &              DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13,
     &              WFTYP
C
C *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************
C
      DOUBLE PRECISION MOLAL
      DOUBLE PRECISION MOLALR
      DOUBLE PRECISION M0
      DOUBLE PRECISION IONIC
      DOUBLE PRECISION GAMA, ZZ, Z, GAMOU, GAMIN, GASAQ, EPSACT
      DOUBLE PRECISION COH, CHNO3, CHCL, WATER
      INTEGER IACALC
      LOGICAL CALAOU, CALAIN, FRST, DRYF
      COMMON /IONS/ MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR), ZZ(NPAIR),
     &              Z(NIONS),     GAMOU(NPAIR),  GAMIN(NPAIR),M0(NPAIR),
     &              GASAQ(NGASAQ),
     &              EPSACT,       COH,           CHNO3,       CHCL,         
     &              WATER,        IONIC,         IACALC,      
     &              FRST,         CALAIN,        CALAOU,      DRYF
C
C *** VARIABLES FOR SOLID AEROSOL PHASE ********************************
C
      DOUBLE PRECISION CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4
      DOUBLE PRECISION CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC
      COMMON /SALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4, 
     &              CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC
C
C *** VARIABLES FOR GAS PHASE ******************************************
C
      DOUBLE PRECISION GNH3, GHNO3, GHCL 
      COMMON /GAS / GNH3, GHNO3, GHCL 
C
C *** EQUILIBRIUM CONSTANTS ********************************************
C
      DOUBLE PRECISION XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10
      DOUBLE PRECISION XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41
      DOUBLE PRECISION XK42
C
      COMMON /EQUK/ XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10,
     &              XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41,
     &              XK42
C
C *** MOLECULAR WEIGHTS ************************************************
C
      DOUBLE PRECISION IMW, R, WMW, SMW
      COMMON /OTHR/ R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)
C
C *** SOLUTION/INFO VARIABLES ******************************************
C
      CHARACTER(15) scase
      DOUBLE PRECISION SULRATW, SULRAT, SODRAT
      COMMON /CASE/ SULRATW, SULRAT, SODRAT, SCASE
C
      DOUBLE PRECISION EPS
      INTEGER MAXIT, NSWEEP, NDIV, ICLACT
      LOGICAL GINIT
      COMMON /SOLN/ EPS, MAXIT, NSWEEP, NDIV, ICLACT, GINIT
C
C *** ERROR SYSTEM *****************************************************
C
      CHARACTER(40) errmsg
      INTEGER   ERRSTK, NOFER   
      LOGICAL   STKOFL, NONPHYS, tayser      
      COMMON /EROR/ STKOFL, NOFER, ERRSTK(NERRMX), ERRMSG(NERRMX),
     &              tayser, NONPHYS
C
C *** GENERIC VARIABLES ************************************************
C
      CHARACTER(15) version
      DOUBLE PRECISION GREAT, TINY, TINY2, ZERO, ONE
      COMMON /CGEN/ GREAT, TINY, TINY2, ZERO, ONE, VERSION
C
C *** For differentiated routine ***************************************
C
      DOUBLE PRECISION psi1
      DOUBLE PRECISION psi2
      DOUBLE PRECISION psi3
      DOUBLE PRECISION psi4
      DOUBLE PRECISION psi5
      DOUBLE PRECISION psi6
      DOUBLE PRECISION psi7
      DOUBLE PRECISION psi8
      DOUBLE PRECISION a1
      DOUBLE PRECISION a2
      DOUBLE PRECISION chi1
      DOUBLE PRECISION a3
      DOUBLE PRECISION chi2
      DOUBLE PRECISION a4
      DOUBLE PRECISION chi3
      DOUBLE PRECISION a5
      DOUBLE PRECISION chi4
      DOUBLE PRECISION a6
      DOUBLE PRECISION chi5
      DOUBLE PRECISION a7
      DOUBLE PRECISION chi6
      DOUBLE PRECISION a8
      DOUBLE PRECISION chi7
      DOUBLE PRECISION chi8
      COMMON /solut/ chi1, chi2, chi3, chi4, chi5, chi6, chi7, chi8, 
     +psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8, a1, a2, a3, a4, a5
     +, a6, a7, a8
C
      DOUBLE PRECISION wb(ncomp)
      COMMON /inpt_b/ wb
C
C *** Case A adjoint ***************************************************
C
      DOUBLE PRECISION molalab(nions)
      DOUBLE PRECISION gamaab(npair)
      DOUBLE PRECISION waterab
      COMMON /ions_ab/ molalab, gamaab, waterab
      DOUBLE PRECISION wab(ncomp)
      COMMON /inpt_ab/ wab
      DOUBLE PRECISION gnh3ab
      COMMON /gas_ab/ gnh3ab
CC
CC *** Case B adjoint ***************************************************
CC
      DOUBLE PRECISION molalbb(nions)
      DOUBLE PRECISION gamabb(npair)
      DOUBLE PRECISION waterbb
      COMMON /ions_bb/ molalbb, gamabb, waterbb
      DOUBLE PRECISION gnh3bb
      COMMON /gas_bb/  gnh3bb
CC
CC *** Case C adjoint ***************************************************
CC
      DOUBLE PRECISION molalcb(nions)
      DOUBLE PRECISION gamacb(npair)
      DOUBLE PRECISION watercb
      COMMON /ions_cb/ molalcb, gamacb, watercb
      DOUBLE PRECISION gnh3cb
      COMMON /gas_cb/ gnh3cb
CC
CC *** Case D adjoint ***************************************************
CC
      DOUBLE PRECISION molaldnrd(nions)
      DOUBLE PRECISION gamadnrd(npair)
      DOUBLE PRECISION waterdnrd
      COMMON /ions_dnrd/ molaldnrd, gamadnrd, waterdnrd
C
      DOUBLE PRECISION chi3db
      DOUBLE PRECISION chi4db
      DOUBLE PRECISION psi1db
      DOUBLE PRECISION psi2db
      COMMON /solut_db/ chi3db, chi4db, psi1db, psi2db
      DOUBLE PRECISION molaldnrddb(nions)
      DOUBLE PRECISION gamadnrddb(npair)
      DOUBLE PRECISION waterdnrddb
      COMMON /ions_dnrd_db/ molaldnrddb, gamadnrddb, 
     +waterdnrddb
      DOUBLE PRECISION gnh3db, ghno3db
      COMMON /gas_db/ gnh3db, ghno3db
      DOUBLE PRECISION molaldb(nions)
      DOUBLE PRECISION molalrdb(npair)
      DOUBLE PRECISION gamadb(npair)  
      DOUBLE PRECISION waterdb
      COMMON /ions_db/ molaldb, molalrdb, gamadb, waterdb
CC
CC *** Case E adjoint ***************************************************
CC
      DOUBLE PRECISION molaleb(nions)
      DOUBLE PRECISION gamaeb(npair)
      DOUBLE PRECISION watereb
      COMMON /ions_eb/ molaleb, gamaeb, watereb
      DOUBLE PRECISION web(ncomp)
      COMMON /inpt_eb/ web
      DOUBLE PRECISION ghno3eb
      COMMON /gas_eb/ ghno3eb
CC
CC *** Case F adjoint ***************************************************
CC
      DOUBLE PRECISION molalfb(nions)
      DOUBLE PRECISION gamafb(npair)
      DOUBLE PRECISION waterfb
      COMMON /ions_fb/ molalfb, gamafb, waterfb
      DOUBLE PRECISION wfb(ncomp)
      COMMON /inpt_fb/ wfb
      DOUBLE PRECISION ghno3fb
      COMMON /gas_fb/ ghno3fb
CC
CC *** Case G adjoint ***************************************************
CC
      DOUBLE PRECISION molalgnrd(nions)
      DOUBLE PRECISION gamagnrd(npair)
      DOUBLE PRECISION watergnrd
      COMMON /ions_gnrd/ molalgnrd, gamagnrd, watergnrd   
      DOUBLE PRECISION molalgnrdgb(nions)
      DOUBLE PRECISION gamagnrdgb(npair)
      DOUBLE PRECISION watergnrdgb
      COMMON /ions_gnrd_gb/ molalgnrdgb, gamagnrdgb, watergnrdgb
      DOUBLE PRECISION chi4gb, chi5gb, chi6gb, psi2gb
      COMMON /solut_gb/ chi4gb, chi5gb, chi6gb, psi2gb
      DOUBLE PRECISION gnh3gb, ghno3gb, ghclgb
      COMMON /gas_gb/ gnh3gb, ghno3gb, ghclgb
      DOUBLE PRECISION molalgb(nions)
      DOUBLE PRECISION molalrgb(npair)
      DOUBLE PRECISION gamagb(npair)
      DOUBLE PRECISION watergb
      COMMON /ions_gb/ molalgb, molalrgb, gamagb, watergb
      DOUBLE PRECISION wgb(5)
      COMMON /inpt_gb/ wgb
CC
CC *** Case H adjoint ***************************************************
CC
      DOUBLE PRECISION molalhnrd(nions)
      DOUBLE PRECISION gamahnrd(npair)
      DOUBLE PRECISION waterhnrd
      COMMON /ions_hnrd/ molalhnrd, gamahnrd, waterhnrd
      DOUBLE PRECISION molalhnrdhb(nions)
      DOUBLE PRECISION gamahnrdhb(npair)
      DOUBLE PRECISION waterhnrdhb
      COMMON /ions_hnrd_hb/ molalhnrdhb, gamahnrdhb, waterhnrdhb
      DOUBLE PRECISION chi1hb, chi4hb, chi5hb, chi6hb, chi7hb, chi8hb
      COMMON /solut_hb/ chi1hb, chi4hb, chi5hb, chi6hb, chi7hb, chi8hb
      DOUBLE PRECISION gnh3hb, ghno3hb, ghclhb
      COMMON /gas_hb/ gnh3hb, ghno3hb, ghclhb
      DOUBLE PRECISION molalhb(nions)
      DOUBLE PRECISION molalrhb(npair)
      DOUBLE PRECISION gamahb(npair)
      DOUBLE PRECISION waterhb
      COMMON /ions_hb/ molalhb, molalrhb, gamahb, waterhb
CC
CC *** Case I adjoint ***************************************************
CC
      DOUBLE PRECISION gnh3ib, ghno3ib, ghclib
      COMMON /gas_ib/ gnh3ib, ghno3ib, ghclib
      DOUBLE PRECISION molalib(nions)
      DOUBLE PRECISION gamaib(npair)
      DOUBLE PRECISION waterib
      COMMON /ions_ib/ molalib, gamaib, waterib
      DOUBLE PRECISION wib(5)
      COMMON /inpt_ib/ wib
CC
CC *** Case J adjoint ***************************************************
CC
      DOUBLE PRECISION gnh3jb, ghno3jb, ghcljb
      COMMON /gas_jb/ gnh3jb, ghno3jb, ghcljb
      DOUBLE PRECISION molaljb(nions)
      DOUBLE PRECISION gamajb(npair)
      DOUBLE PRECISION waterjb
      COMMON /ions_jb/ molaljb, gamajb, waterjb
      DOUBLE PRECISION wjb(5)
      COMMON /inpt_jb/ wjb
C
C *** END OF INCLUDE FILE **********************************************
C