
!***********************************************************************
!   Portions of Models-3/CMAQ software were developed or based on      *
!   information from various groups: Federal Government employees,     *
!   contractors working on a United States Government contract, and    *
!   non-Federal sources (including research institutions).  These      *
!   research institutions have given the Government permission to      *
!   use, prepare derivative works, and distribute copies of their      *
!   work in Models-3/CMAQ to the public and to permit others to do     *
!   so.  EPA therefore grants similar permissions for use of the       *
!   Models-3/CMAQ software, but users are requested to provide copies  *
!   of derivative works to the Government without restrictions as to   *
!   use by others.  Users are responsible for acquiring their own      *
!   copies of commercial software associated with Models-3/CMAQ and    *
!   for complying with vendor requirements.  Software copyrights by    *
!   the MCNC Environmental Modeling Center are used with their         *
!   permissions subject to the above restrictions.                     *
!***********************************************************************

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE INTRACOAG_GH_ADJ  ( LAMDA, KFM, KFM_ADJ, KNC, DG, 
     &                             DG_ADJ, XLNSIG, XLNSIG_ADJ, QUADS11, 
     &                             QUADS11_ADJ, QUADN11, QUADN11_ADJ )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of subroutine that does Gauss-Hermite numerical quadrature
!   to calculate intramodal coagulation rates for number and 2nd moment
!
! INPUTS:
!   LAMDA, KFM, KNC, DG, XLNSIG, QUADS11_ADJ, QUADN11_ADJ
!
! OUTPUTS:
!   KFM_ADJ, DG_ADJ, XLNSIG_ADJ
!
! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL*8      LAMDA   ! mean free path [ m ]
  
      INTEGER     I, J

      REAL*8      KFM, KNC
      REAL        DG, XLNSIG
      REAL*8      QUADS11, QUADN11
      REAL*8, PARAMETER   ::  PI = 3.14159265358979324d0
      REAL*8, PARAMETER   ::  TWO3RDS = 2.0D0 / 3.0D0
      REAL*8, PARAMETER   ::  SQRT2 = 1.414213562373095D0
      REAL*8      SUM1SFM, SUM2SFM, SUM1NFM, SUM2NFM
      REAL*8      SUM1SNC, SUM2SNC, SUM1NNC, SUM2NNC
      REAL*8      XI, WXI, XF, DP1P, DP1M, DP1PSQ, DP1MSQ
      REAL*8      V1P, V1M, V2P, V2M, A2P, A2M
      REAL*8      YI, WYI, YF, DP2P, DP2M, DP2PSQ, DP2MSQ
      REAL*8      DSPP, DSMP, DSPM, DSMM
      REAL*8      BPPFM, BMPFM, BPMFM, BMMFM
      REAL*8      BPPNC, BMPNC, BPMNC, BMMNC
      REAL*8      xx1, xx2
      REAL*8      XBSFM, XBSNC, XBNFM, XBNNC
      REAL*8      BETAFM, BETANC

      REAL*8, PARAMETER   ::  A = 1.246D0  ! approx Cunningham corr. factor
      
      REAL*8, PARAMETER   ::  TWOA = 2.0D0 * A  

! *** Has a fixed number of Gauss-Herimite abscissas ( n )
      INTEGER, PARAMETER   ::  N = 5   ! one-half the number of abscissas
      REAL*8      GHXI(N)  ! Gauss-Hermite abscissas
      REAL*8      GHWI(N)  ! Gauss-Hermite weights

! adjoint variables
      REAL*8 KFM_ADJ    ! intent( out )
      REAL   DG_ADJ, XLNSIG_ADJ  ! intent( out )
      REAL*8 QUADS11_ADJ, QUADN11_ADJ    ! intent( in )
      REAL*8 SUM1SFM_ADJ, SUM2SFM_ADJ, SUM1NFM_ADJ, SUM2NFM_ADJ
      REAL*8 SUM1SNC_ADJ, SUM2SNC_ADJ, SUM1NNC_ADJ, SUM2NNC_ADJ
      REAL*8 XF_ADJ, DP1P_ADJ, DP1M_ADJ, DP1PSQ_ADJ, DP1MSQ_ADJ
      REAL*8 V1P_ADJ, V1M_ADJ, A2P_ADJ, A2M_ADJ, V2P_ADJ, V2M_ADJ
      REAL*8 YF_ADJ, DP2P_ADJ, DP2M_ADJ, DP2PSQ_ADJ, DP2MSQ_ADJ
      REAL*8 DSPP_ADJ, DSMP_ADJ, DSPM_ADJ, DSMM_ADJ
      REAL*8 BPPFM_ADJ, BMPFM_ADJ, BPMFM_ADJ, BMMFM_ADJ
      REAL*8 BPPNC_ADJ, BMPNC_ADJ, BPMNC_ADJ, BMMNC_ADJ
      REAL*8 XBSFM_ADJ, XBSNC_ADJ, XBNFM_ADJ, XBNNC_ADJ

! ** Values from Table 25.10 (page 924) of Abramowitz and Stegun,
!    Handbook of Mathematical Functions, National Bureau of Standards,
!  December 1965.
!
!    breaks in number to facilitate comparison with printed table
!
! *** tests show that 10 point is adquate.

      DATA GHXI/0.34290 13272 23705D0,
     &          1.03661 08297 89514D0,
     &          1.75668 36492 99882D0,
     &          2.53273 16742 32790D0,
     &          3.43615 91188 37738D0/,

     
     1     GHWI/6.10862 63373 53D-1,
     &          2.40138 61108 23D-1,
     &          3.38743 94455 48D-2,
     &          1.34364 57467 81D-3,
     &          7.64043 28552 33D-6/

!----------------------------------------------------

! *** The following expressions are from Binkowski & Shanker
!     Jour. Geophys. Research. Vol. 100,no. d12, pp 26,191-26,209
!     December 20, 1995
      
! ***  for Free Molecular Eq. A5
        betafm(xx1, xx2) = KFM *
     &       sqrt(1.D0 / xx1**3  + 1.D0 / xx2**3 ) * (xx1 + xx2)**2
      
! ***  for Near Continuum  Eq. A6
        betanc(xx1, xx2) =  KNC * (xx1 + xx2) *
     &                       ( 1.0D0 / xx1 + 1.0D0 / xx2  +
     &                     TWOA * LAMDA * ( 1.0D0 / xx1 ** 2
     &                                    + 1.0D0 / xx2 **2 ) )

      ! Initialize adj variables:
      KFM_ADJ = 0.0
      DG_ADJ = 0.0
      XLNSIG_ADJ = 0.0
      SUM1SFM_ADJ = 0.0
      SUM2SFM_ADJ = 0.0
      SUM1NFM_ADJ = 0.0
      SUM2NFM_ADJ = 0.0
      SUM1SNC_ADJ = 0.0
      SUM2SNC_ADJ = 0.0
      SUM1NNC_ADJ = 0.0
      SUM2NNC_ADJ = 0.0
      XF_ADJ = 0.0
      DP1P_ADJ = 0.0
      DP1M_ADJ = 0.0
      DP1PSQ_ADJ = 0.0
      DP1MSQ_ADJ = 0.0
      V1P_ADJ = 0.0
      V1M_ADJ = 0.0
      A2P_ADJ = 0.0
      A2M_ADJ = 0.0
      V2P_ADJ = 0.0
      V2M_ADJ = 0.0
      YF_ADJ = 0.0
      DP2P_ADJ = 0.0
      DP2M_ADJ = 0.0
      DP2PSQ_ADJ = 0.0
      DP2MSQ_ADJ = 0.0
      DSPP_ADJ = 0.0
      DSMP_ADJ = 0.0
      DSPM_ADJ = 0.0
      DSMM_ADJ = 0.0
      BPPFM_ADJ = 0.0
      BMPFM_ADJ = 0.0
      BPMFM_ADJ = 0.0
      BMMFM_ADJ = 0.0
      BPPNC_ADJ = 0.0
      BPMNC_ADJ = 0.0
      BMPNC_ADJ = 0.0
      BMMNC_ADJ = 0.0
      XBSFM_ADJ = 0.0
      XBSNC_ADJ = 0.0
      XBNFM_ADJ = 0.0
      XBNNC_ADJ = 0.0

! Forward Code:

      sum1sfm = 0.D0
      sum1snc = 0.D0

      sum1nfm = 0.D0
      sum1nnc = 0.D0
      do 1 i=1,n

        sum2sfm = 0.D0
        sum2snc = 0.D0
        sum2nfm = 0.D0
        sum2nnc = 0.D0

        xi = ghxi(i)
        wxi = ghwi(i)
        xf = exp( sqrt2 * xi *xlnsig)
        dp1p = dg*xf
        dp1m = dg/xf
        dp1psq = dp1p*dp1p
        dp1msq = dp1m*dp1m
        v1p = dp1p*dp1psq
        v1m = dp1m*dp1msq
      do 11 j=1,n
        yi = ghxi(j)
        wyi = ghwi(j)
        yf = exp( sqrt2 * yi * xlnsig)
        dp2p = dg*yf
        dp2m = dg/yf
        dp2psq = dp2p*dp2p
        dp2msq = dp2m*dp2m
        a2p = dp2psq
        a2m = dp2msq
        v2p =  dp2p*dp2psq
        v2m =dp2m*dp2msq
        dspp = 0.5D0*(v1p+v2p)**two3rds - a2p
        dsmp = 0.5D0*(v1m+v2p)**two3rds - a2p
        dspm = 0.5D0*(v1p+v2m)**two3rds - a2m
        dsmm = 0.5D0*(v1m+v2m)**two3rds - a2m

        bppfm = betafm(dp1p,dp2p)
        bmpfm = betafm(dp1m,dp2p)
        bpmfm = betafm(dp1p,dp2m)
        bmmfm = betafm(dp1m,dp2m)

        bppnc = betanc(dp1p,dp2p)
        bmpnc = betanc(dp1m,dp2p)
        bpmnc = betanc(dp1p,dp2m)
        bmmnc = betanc(dp1m,dp2m)
        
        sum2sfm = sum2sfm + wyi*(dspp * bppfm + dspm * bpmfm
     &               +   dsmp * bmpfm + dsmm * bmmfm )

        sum2nfm = sum2nfm + wyi*(bppfm + bmpfm + bpmfm + bmmfm)
        sum2snc = sum2snc + wyi*(dspp * bppnc + dspm * bpmnc
     &               +   dsmp * bmpnc + dsmm * bmmnc ) 
        sum2nnc = sum2nnc + wyi*(bppnc + bmpnc + bpmnc + bmmnc)
   
   11 continue
      sum1sfm = sum1sfm + wxi * sum2sfm
      sum1nfm = sum1nfm + wxi * sum2nfm
      
      sum1snc = sum1snc + wxi * sum2snc
      sum1nnc = sum1nnc + wxi * sum2nnc
    
    1 continue
      
      xbsfm   = -sum1sfm  / pi
      xbsnc   = -sum1snc  / pi
      
      quads11 =  xbsfm * xbsnc / ( xbsfm + xbsnc )

c *** quads11 is the intra-modal coagulation term for 2nd moment
      
      xbnfm   = 0.5D0 * sum1nfm  / pi
      xbnnc   = 0.5D0 * sum1nnc  / pi

      
      quadn11 =  xbnfm * xbnnc / ( xbnfm + xbnnc )

c *** quadn11 is the intra-modal coagulation term for number

      !------
      ! fwd code:
      ! quadn11 = xbnfm * xbnnc / ( xbnfm + xbnnc )
      ! adj code:
      XBNFM_ADJ = ( ( XBNNC / ( XBNFM + XBNNC ) )
     &          - ( XBNFM * XBNNC / ( ( XBNFM + XBNNC ) ** 2.0 ) )
     &          ) * QUADN11_ADJ
      XBNNC_ADJ = ( ( XBNFM / ( XBNFM + XBNNC ) )
     &          - ( XBNFM * XBNNC / ( ( XBNFM + XBNNC ) ** 2.0 ) )
     &          ) * QUADN11_ADJ
      QUADN11_ADJ = 0.0

      !------
      ! fwd code:
      ! xbnfm   = 0.5D0 * sum1nfm  / pi
      ! xbnnc   = 0.5D0 * sum1nnc  / pi
      ! adj code:
      SUM1NNC_ADJ = ( 0.5D0 / PI ) * XBNNC_ADJ
      XBNNC_ADJ = 0.0
      SUM1NFM_ADJ = ( 0.5D0 / PI ) * XBNFM_ADJ
      XBNFM_ADJ = 0.0

      !------
      ! fwd code:
      ! quads11 =  xbsfm * xbsnc / ( xbsfm + xbsnc )
      ! adj code:
      XBSFM_ADJ = ( XBSNC / ( XBSFM + XBSNC ) 
     &          - XBSFM * XBSNC / ( ( XBSFM + XBSNC ) ** 2.0 ) 
     &          ) * QUADS11_ADJ
      XBSNC_ADJ = ( XBSFM / ( XBSFM + XBSNC ) 
     &          - XBSFM * XBSNC / ( ( XBSFM + XBSNC ) ** 2.0 ) 
     &          ) * QUADS11_ADJ
      QUADS11_ADJ = 0.0

      !------
      ! fwd code:
      ! xbsfm   = -sum1sfm  / pi
      ! xbsnc   = -sum1snc  / pi
      ! adj code:
      SUM1SNC_ADJ = - 1 / PI * XBSNC_ADJ
      XBSNC_ADJ = 0.0
      SUM1SFM_ADJ = - 1 / PI * XBSFM_ADJ
      XBSFM_ADJ = 0.0
 
      DO I = N, 1, -1
        
        XI = GHXI( I )
        XF = EXP( SQRT2 * XI * XLNSIG )
        DP1P = DG * XF
        DP1M = DG / XF
        DP1PSQ = DP1P * DP1P
        DP1MSQ = DP1M * DP1M
        V1P = DP1P * DP1PSQ 
        V1M = DP1M * DP1MSQ 
        WXI = GHWI( I ) 

        !------
        ! fwd code:
        ! sum1sfm = sum1sfm + wxi * sum2sfm
        ! sum1nfm = sum1nfm + wxi * sum2nfm
        ! sum1snc = sum1snc + wxi * sum2snc
        ! sum1nnc = sum1nnc + wxi * sum2nnc
        ! adj code:
        SUM2NNC_ADJ = WXI * SUM1NNC_ADJ
        SUM2SNC_ADJ = WXI * SUM1SNC_ADJ
        SUM2NFM_ADJ = WXI * SUM1NFM_ADJ
        SUM2SFM_ADJ = WXI * SUM1SFM_ADJ

        DO J = N, 1, -1
           
           YI = GHXI( J ) 
           YF = EXP( SQRT2 * YI * XLNSIG )
           DP2P = DG * YF 
           DP2M = DG / YF
           DP2PSQ = DP2P * DP2P
           DP2MSQ = DP2M * DP2M
           A2P = DP2PSQ
           A2M = DP2MSQ
           V2P = DP2P * DP2PSQ
           V2M = DP2M * DP2MSQ
           WYI = GHWI( J ) 
           DSPP = 0.5D0*(V1P+V2P)**TWO3RDS - A2P
           DSMP = 0.5D0*(V1M+V2P)**TWO3RDS - A2P
           DSPM = 0.5D0*(V1P+V2M)**TWO3RDS - A2M
           DSMM = 0.5D0*(V1M+V2M)**TWO3RDS - A2M

           BPPFM = BETAFM(DP1P,DP2P)
           BMPFM = BETAFM(DP1M,DP2P)
           BPMFM = BETAFM(DP1P,DP2M)
           BMMFM = BETAFM(DP1M,DP2M)

           BPPNC = BETANC(DP1P,DP2P)
           BMPNC = BETANC(DP1M,DP2P)
           BPMNC = BETANC(DP1P,DP2M)
           BMMNC = BETANC(DP1M,DP2M)

           !------
           ! fwd code:
           ! sum2nnc = sum2nnc + wyi*(bppnc + bmpnc + bpmnc + bmmnc)
           ! adj code:
           BPPNC_ADJ = WYI * SUM2NNC_ADJ
           BMPNC_ADJ = WYI * SUM2NNC_ADJ
           BPMNC_ADJ = WYI * SUM2NNC_ADJ
           BMMNC_ADJ = WYI * SUM2NNC_ADJ

           !------
           ! fwd code:
           ! sum2snc = sum2snc + wyi*(dspp * bppnc + dspm * bpmnc
           ! &                        + dsmp * bmpnc + dsmm * bmmnc )
           ! adj code:
           DSPP_ADJ = WYI * BPPNC * SUM2SNC_ADJ
           BPPNC_ADJ = WYI * DSPP * SUM2SNC_ADJ
           DSPM_ADJ = WYI * BPMNC * SUM2SNC_ADJ
           BPMNC_ADJ = WYI * DSPM * SUM2SNC_ADJ
           DSMP_ADJ = WYI * BMPNC * SUM2SNC_ADJ
           BMPNC_ADJ = WYI * DSMP * SUM2SNC_ADJ
           DSMM_ADJ = WYI * BMMNC * SUM2SNC_ADJ
           BMMNC_ADJ = WYI * DSMM * SUM2SNC_ADJ
    
           !------
           ! fwd code:
           ! sum2nfm = sum2nfm + wyi*(bppfm + bmpfm + bpmfm + bmmfm)
           ! adj code:
           BPPFM_ADJ = WYI * SUM2NFM_ADJ
           BMPFM_ADJ = WYI * SUM2NFM_ADJ
           BPMFM_ADJ = WYI * SUM2NFM_ADJ
           BMMFM_ADJ = WYI * SUM2NFM_ADJ

           !------
           ! fwd code:
           ! sum2sfm = sum2sfm + wyi*(dspp * bppfm + dspm * bpmfm
           ! &                        + dsmp * bmpfm + dsmm * bmmfm )
           ! adj code:
           DSPP_ADJ = DSPP_ADJ + WYI * BPPFM * SUM2SFM_ADJ
           BPPFM_ADJ = BPPFM_ADJ + WYI * DSPP * SUM2SFM_ADJ
           DSPM_ADJ = DSPM_ADJ + WYI * BPMFM * SUM2SFM_ADJ
           BPMFM_ADJ = BPMFM_ADJ + WYI * DSPM * SUM2SFM_ADJ
           DSMP_ADJ = DSMP_ADJ + WYI * BMPFM * SUM2SFM_ADJ
           BMPFM_ADJ = BMPFM_ADJ + WYI * DSMP * SUM2SFM_ADJ
           DSMM_ADJ = DSMM_ADJ + WYI * BMMFM * SUM2SFM_ADJ
           BMMFM_ADJ = BMMFM_ADJ + WYI * DSMM * SUM2SFM_ADJ

           !------
           ! fwd code:
           ! bppnc = betanc(dp1p,dp2p)
           ! bmpnc = betanc(dp1m,dp2p)
           ! bpmnc = betanc(dp1p,dp2m)
           ! bmmnc = betanc(dp1m,dp2m)
           ! adj code:
           DP1M_ADJ = DP1M_ADJ + ( ( DP1M + DP2M ) * KNC * 
     &                ( (- 1.0 / ( DP1M ** 2.0 ) ) - ( 2.0 * LAMDA * 
     &                TWOA / ( DP1M ** 3.0 ) ) ) + KNC * ( ( 1.0 / DP1M
     &                ) + ( 1.0 / DP2M ) + LAMDA * TWOA * ( ( 
     &                1.0 / ( DP1M ** 2.0 ) ) + ( 1.0 / ( DP2M ** 2.0 ) 
     &                ) ) ) ) * BMMNC_ADJ
           DP2M_ADJ = DP2M_ADJ + ( KNC * ( ( 1.0 / DP1M ) + 
     &                ( 1.0 / DP2M ) + ( ( 1.0 / ( DP1M ** 2.0 ) ) +
     &                ( 1.0 / ( DP2M ** 2.0 ) ) ) * TWOA * LAMDA ) + 
     &                ( DP1M + DP2M ) * KNC * ( ( -1.0 / ( DP2M ** 2.0 
     &                ) ) - ( 2.0 * LAMDA * TWOA / ( DP2M ** 2.0 ) ) ) 
     &                ) * BMMNC_ADJ
           BMMNC_ADJ = 0.0
           DP1P_ADJ = DP1P_ADJ + ( ( DP1P + DP2M ) * KNC * 
     &                ( (- 1.0 / ( DP1P ** 2.0 ) ) - ( 2.0 * LAMDA *
     &                TWOA / ( DP1P ** 3.0 ) ) ) + KNC * ( ( 1.0 / DP1P
     &                ) + ( 1.0 / DP2M ) + LAMDA * TWOA * ( ( 
     &                1.0 / ( DP1P ** 2.0 ) ) + ( 1.0 / ( DP2M ** 2.0 ) 
     &                ) ) ) ) * BPMNC_ADJ
           DP2M_ADJ = DP2M_ADJ + ( KNC * ( ( 1.0 / DP1P ) + 
     &                ( 1.0 / DP2M ) + ( ( 1.0 / ( DP1P ** 2.0 ) ) +
     &                ( 1.0 / ( DP2M ** 2.0 ) ) ) * TWOA * LAMDA ) + 
     &                ( DP1P + DP2M ) * KNC * ( ( -1.0 / ( DP2M ** 2.0 
     &                ) ) - ( 2.0 * LAMDA * TWOA / ( DP2M ** 2.0 ) ) ) 
     &                ) * BPMNC_ADJ
           BPMNC_ADJ = 0.0
           DP1M_ADJ = DP1M_ADJ + ( ( DP1M + DP2P ) * KNC * 
     &                ( (- 1.0 / ( DP1M ** 2.0 ) ) - ( 2.0 * LAMDA * 
     &                TWOA / ( DP1M ** 3.0 ) ) ) + KNC * ( ( 1.0 / DP1M
     &                ) + ( 1.0 / DP2P ) + LAMDA * TWOA * ( ( 
     &                1.0 / ( DP1M ** 2.0 ) ) + ( 1.0 / ( DP2P ** 2.0 ) 
     &                ) ) ) ) * BMPNC_ADJ
           DP2P_ADJ = DP2P_ADJ + ( KNC * ( ( 1.0 / DP1M ) + 
     &                ( 1.0 / DP2P ) + ( ( 1.0 / ( DP1M ** 2.0 ) ) +
     &                ( 1.0 / ( DP2P ** 2.0 ) ) ) * TWOA * LAMDA ) + 
     &                ( DP1M + DP2P ) * KNC * ( ( -1.0 / ( DP2P ** 2.0 
     &                ) ) - ( 2.0 * LAMDA * TWOA / ( DP2P ** 2.0 ) ) ) 
     &                ) * BMPNC_ADJ
           BMPNC_ADJ = 0.0
           DP1P_ADJ = DP1P_ADJ + ( ( DP1P + DP2P ) * KNC * 
     &                ( (- 1.0 / ( DP1P ** 2.0 ) ) - ( 2.0 * LAMDA * 
     &                TWOA / ( DP1P ** 3.0 ) ) ) + KNC * ( ( 1.0 / DP1P
     &                ) + ( 1.0 / DP2P ) + LAMDA * TWOA * ( ( 
     &                1.0 / ( DP1P ** 2.0 ) ) + ( 1.0 / ( DP2P ** 2.0 ) 
     &                ) ) ) ) * BPPNC_ADJ
           DP2P_ADJ = DP2P_ADJ + ( KNC * ( ( 1.0 / DP1P ) + 
     &                ( 1.0 / DP2P ) + ( ( 1.0 / ( DP1P ** 2.0 ) ) +
     &                ( 1.0 / ( DP2P ** 2.0 ) ) ) * TWOA * LAMDA ) + 
     &                ( DP1P + DP2P ) * KNC * ( ( -1.0 / ( DP2P ** 2.0 
     &                ) ) - ( 2.0 * LAMDA * TWOA / ( DP2P ** 2.0 ) ) ) 
     &                ) * BPPNC_ADJ
           BPPNC_ADJ = 0.0

           !------
           ! fwd code:
           ! bppfm = betafm(dp1p,dp2p)
           ! bmpfm = betafm(dp1m,dp2p)
           ! bpmfm = betafm(dp1p,dp2m)
           ! bmmfm = betafm(dp1m,dp2m)
           ! adj code:
           DP1M_ADJ = DP1M_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1M ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) * ( DP1M + DP2M ) 
     &              * KFM - ( 3.0 * ( ( DP1M + DP2M ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP1M ** 4.0 ) * SQRT( ( 1.0 / ( DP1M 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) ) ) ) 
     &              * BMMFM_ADJ
           DP2M_ADJ = DP2M_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1M ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) * ( DP1M + DP2M ) 
     &              * KFM - ( 3.0 * ( ( DP1M + DP2M ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP2M ** 4.0 ) * SQRT( ( 1.0 / ( DP1M 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) ) ) ) 
     &              * BMMFM_ADJ
           KFM_ADJ = KFM_ADJ + ( SQRT( ( 1.0 / ( DP1M ** 3.0 ) ) 
     &             + ( 1.0 / ( DP2M ** 3.0 ) ) ) * ( ( DP1M + DP2M 
     &             ) ** 2.0 ) ) * BMMFM_ADJ
           BMMFM_ADJ = 0.0
           DP1P_ADJ = DP1P_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1P ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) * ( DP1P + DP2M ) 
     &              * KFM - ( 3.0 * ( ( DP1P + DP2M ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP1P ** 4.0 ) * SQRT( ( 1.0 / ( DP1P 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) ) ) ) 
     &              * BPMFM_ADJ
           DP2M_ADJ = DP2M_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1P ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) * ( DP1P + DP2M )
     &              * KFM - ( 3.0 * ( ( DP1P + DP2M ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP2M ** 4.0 ) * SQRT( ( 1.0 / ( DP1P 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2M ** 3.0 ) ) ) ) ) ) 
     &              * BPMFM_ADJ
           KFM_ADJ = KFM_ADJ + ( SQRT( ( 1.0 / ( DP1P ** 3.0 ) ) 
     &             + ( 1.0 / ( DP2M ** 3.0 ) ) ) * ( ( DP1P + DP2M 
     &             ) ** 2.0 ) ) * BPMFM_ADJ
           BPMFM_ADJ = 0.0
           DP1M_ADJ = DP1M_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1M ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) * ( DP1M + DP2P ) 
     &              * KFM - ( 3.0 * ( ( DP1M + DP2P ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP1M ** 4.0 ) * SQRT( ( 1.0 / ( DP1M 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) ) ) ) 
     &              * BMPFM_ADJ
           DP2P_ADJ = DP2P_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1M ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) * ( DP1M + DP2P ) 
     &              * KFM - ( 3.0 * ( ( DP1M + DP2P ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP2P ** 4.0 ) * SQRT( ( 1.0 / ( DP1M 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) ) ) ) 
     &              * BMPFM_ADJ
           KFM_ADJ = KFM_ADJ + ( SQRT( ( 1.0 / ( DP1M ** 3.0 ) ) 
     &             + ( 1.0 / ( DP2P ** 3.0 ) ) ) * ( ( DP1M + DP2P 
     &             ) ** 2.0 ) ) * BMPFM_ADJ
           BMPFM_ADJ = 0.0
           DP1P_ADJ = DP1P_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1P ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) * ( DP1P + DP2P ) 
     &              * KFM - ( 3.0 * ( ( DP1P + DP2P ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP1P ** 4.0 ) * SQRT( ( 1.0 / ( DP1P 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) ) ) ) 
     &              * BPPFM_ADJ
           DP2P_ADJ = DP2P_ADJ + ( 2.0 * SQRT( ( 1.0 / ( DP1P ** 3.0 ) 
     &              ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) * ( DP1P + DP2P ) 
     &              * KFM - ( 3.0 * ( ( DP1P + DP2P ) ** 2.0 ) * KFM 
     &              / ( 2.0 * ( DP2P ** 4.0 ) * SQRT( ( 1.0 / ( DP1P 
     &              ** 3.0 ) ) + ( 1.0 / ( DP2P ** 3.0 ) ) ) ) ) ) 
     &              * BPPFM_ADJ
           KFM_ADJ = KFM_ADJ + ( SQRT( ( 1.0 / ( DP1P ** 3.0 ) ) 
     &             + ( 1.0 / ( DP2P ** 3.0 ) ) ) * ( ( DP1P + DP2P 
     &             ) ** 2.0 ) ) * BPPFM_ADJ
           BPPFM_ADJ = 0.0

           !------
           ! fwd code:
           ! dsmm = 0.5D0*(v1m+v2m)**two3rds - a2m
           ! adj code:
           A2M_ADJ = A2M_ADJ - DSMM_ADJ
           V1M_ADJ = V1M_ADJ + TWO3RDS * 0.5D0 * ( ( V1M + V2M ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSMM_ADJ
           V2M_ADJ = V2M_ADJ + TWO3RDS * 0.5D0 * ( ( V1M + V2M ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSMM_ADJ 
           DSMM_ADJ = 0.0

           !------ 
           ! fwd code:
           ! dspm = 0.5D0*(v1p+v2m)**two3rds - a2m
           ! adj code:
           A2M_ADJ = A2M_ADJ - DSPM_ADJ
           V1P_ADJ = V1P_ADJ + TWO3RDS * 0.5D0 * ( ( V1P + V2M ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSPM_ADJ
           V2M_ADJ = V2M_ADJ + TWO3RDS * 0.5D0 * ( ( V1P + V2M ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSPM_ADJ
           DSPM_ADJ = 0.0

           !------
           ! fwd code:
           ! dsmp = 0.5D0*(v1m+v2p)**two3rds - a2p
           ! adj code:
           A2P_ADJ = A2P_ADJ - DSMP_ADJ
           V1M_ADJ = V1M_ADJ + TWO3RDS * 0.5D0 * ( ( V1M + V2P ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSMP_ADJ
           V2P_ADJ = V2P_ADJ + TWO3RDS * 0.5D0 * ( ( V1M + V2P ) ** (
     &               TWO3RDS - 1.0 ) ) * DSMP_ADJ
           DSMP_ADJ = 0.0
 
           !------
           ! fwd code:
           ! dspp = 0.5D0*(v1p+v2p)**two3rds - a2p
           ! adj code:
           A2P_ADJ = A2P_ADJ - DSPP_ADJ
           V1P_ADJ = V1P_ADJ + TWO3RDS * 0.5D0 * ( ( V1P + V2P ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSPP_ADJ
           V2P_ADJ = V2P_ADJ + TWO3RDS * 0.5D0 * ( ( V1P + V2P ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSPP_ADJ
           DSPP_ADJ = 0.0
 
           !------
           ! fwd code:
           ! v2p =  dp2p*dp2psq
           ! v2m =dp2m*dp2msq
           ! adj code:
           DP2M_ADJ = DP2M_ADJ + DP2MSQ * V2M_ADJ
           DP2MSQ_ADJ = DP2MSQ_ADJ + DP2M * V2M_ADJ
           V2M_ADJ = 0.0
           DP2P_ADJ = DP2P_ADJ + DP2PSQ * V2P_ADJ
           DP2PSQ_ADJ = DP2PSQ_ADJ + DP2P * V2P_ADJ
           V2P_ADJ = 0.0
  
           !------
           ! fwd code:
           ! a2p = dp2psq
           ! a2m = dp2msq
           ! adj code:
           DP2MSQ_ADJ = DP2MSQ_ADJ + A2M_ADJ
           A2M_ADJ = 0.0
           DP2PSQ_ADJ = DP2PSQ_ADJ + A2P_ADJ
           A2P_ADJ = 0.0

           !------
           ! fwd code:
           ! dp2psq = dp2p*dp2p
           ! dp2msq = dp2m*dp2m
           ! adj code:
           DP2M_ADJ = DP2M_ADJ + 2.0 * DP2M * DP2MSQ_ADJ
           DP2MSQ_ADJ = 0.0
           DP2P_ADJ = DP2P_ADJ + 2.0 * DP2P * DP2PSQ_ADJ
           DP2PSQ_ADJ = 0.0

           !------
           ! fwd code:
           ! dp2p = dg*yf
           ! dp2m = dg/yf
           ! adj code:
           DG_ADJ = DG_ADJ + ( 1 / YF ) * DP2M_ADJ
           YF_ADJ = YF_ADJ - ( DG / ( YF ** 2.0 ) ) * DP2M_ADJ
           DP2M_ADJ = 0.0
           DG_ADJ = DG_ADJ + YF * DP2P_ADJ
           YF_ADJ = YF_ADJ + DG * DP2P_ADJ
           DP2P_ADJ = 0.0

           !------
           ! fwd code:
           ! yf = exp( sqrt2 * yi * xlnsig)
           ! adj code:
           XLNSIG_ADJ = XLNSIG_ADJ + SQRT2 * YI * EXP( SQRT2 * YI * 
     &                  XLNSIG ) * YF_ADJ
           YF_ADJ = 0.0

        END DO

        !------
        ! fwd code:
        ! v1p = dp1p*dp1psq
        ! v1m = dp1m*dp1msq
        ! adj code:
        DP1M_ADJ = DP1M_ADJ + DP1MSQ * V1M_ADJ
        DP1MSQ_ADJ = DP1MSQ_ADJ + DP1M * V1M_ADJ
        V1M_ADJ = 0.0
        DP1P_ADJ = DP1P_ADJ + DP1PSQ * V1P_ADJ
        DP1PSQ_ADJ = DP1PSQ_ADJ + DP1P * V1P_ADJ
        V1P_ADJ = 0.0

        !------
        ! fwd code:
        ! dp1psq = dp1p*dp1p
        ! dp1msq = dp1m*dp1m
        ! adj code:
        DP1M_ADJ = DP1M_ADJ + 2.0 * DP1M * DP1MSQ_ADJ
        DP1MSQ_ADJ = 0.0
        DP1P_ADJ = DP1P_ADJ + 2.0 * DP1P * DP1PSQ_ADJ
        DP1PSQ_ADJ = 0.0

        !------
        ! fwd code:
        ! dp1p = dg*xf
        ! dp1m = dg/xf
        ! adj code:
        DG_ADJ = DG_ADJ + ( 1.0 / XF ) * DP1M_ADJ
        XF_ADJ = XF_ADJ - ( DG / ( XF ** 2.0 ) ) * DP1M_ADJ
        DP1M_ADJ = 0.0
        DG_ADJ = DG_ADJ + XF * DP1P_ADJ
        XF_ADJ = XF_ADJ + DG * DP1P_ADJ
        DP1P_ADJ = 0.0

        !------
        ! fwd code:
        ! xf = exp( sqrt2 * xi *xlnsig)
        ! adj code:
        XLNSIG_ADJ = XLNSIG_ADJ + SQRT2 * XI * EXP( SQRT2 * XI * 
     &               XLNSIG ) * XF_ADJ 
        XF_ADJ = 0.0

        !------
        ! fwd code:
        ! sum2sfm = 0.D0
        ! sum2snc = 0.D0
        ! sum2nfm = 0.D0
        ! sum2nnc = 0.D0
        ! adj code:
        SUM2NNC_ADJ = 0.0
        SUM2NFM_ADJ = 0.0
        SUM2SNC_ADJ = 0.0
        SUM2SFM_ADJ = 0.0

      END DO

      END SUBROUTINE INTRACOAG_GH_ADJ
