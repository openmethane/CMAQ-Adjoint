
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
      SUBROUTINE INTERCOAG_GH_ADJ  ( LAMDA, KFM, KFM_ADJ, KNC, DG1, 
     &                             DG1_ADJ, DG2, DG2_ADJ, XLNSIG1, 
     &                             XLNSIG1_ADJ, XLNSIG2, XLNSIG2_ADJ, 
     &                             QUADS12, QUADS12_ADJ, QUADS21, 
     &                             QUADS21_ADJ, QUADN12, QUADN12_ADJ, 
     &                             QUADV12, QUADV12_ADJ )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of subroutine that does Gauss-Hermite numerical quadrature
!   to calculate intermodal coagulation for number, 2nd, and 3rd moment
!
! INPUTS:
!   LAMDA, KFM, KNC, DG1, DG2, XLNSIG1, XLNSIG2, QUADS12_ADJ, 
!   QUADS21_ADJ, QUADN12_ADJ, QUADV12_ADJ
!
! OUTPUTS:
!   KFM_ADJ, DG1_ADJ, DG2_ADJ, XLNSIG1_ADJ, XLNSIG2_ADJ
!
! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL*8      LAMDA   ! mean free path [ m ]

      INTEGER     I, J

      REAL*8      KFM, KNC
      REAL        DG1, DG2, XLNSIG1, XLNSIG2
      REAL*8      QUADS12, QUADS21, QUADN12, QUADV12
      REAL*8, PARAMETER   ::  PI = 3.14159265358979324d0
      REAL*8, PARAMETER   ::  TWO3RDS = 2.0D0 / 3.0D0
      REAL*8, PARAMETER   ::  SQRT2 = 1.414213562373095D0
      REAL*8      SUM1S12FM, SUM1S21FM, SUM2S12FM, SUM2S21FM
      REAL*8      SUM1NFM, SUM2NFM
      REAL*8      SUM1VFM, SUM2VFM
      REAL*8      SUM1S12NC, SUM1S21NC, SUM2S12NC, SUM2S21NC
      REAL*8      SUM1NNC, SUM2NNC, SUM1VNC, SUM2VNC
      REAL*8      XI, WXI, XF, DP1P, DP1M, DP1PSQ, DP1MSQ
      REAL*8      A1P, A1M, V1P, V1M
      REAL*8      A2P, A2M, V2P, V2M
      REAL*8      YI, WYI, YF, DP2P ,DP2M, DP2PSQ, DP2MSQ
      REAL*8      DSPP, DSMP, DSPM, DSMM
      REAL*8      BPPFM, BMPFM, BPMFM, BMMFM
      REAL*8      BPPNC, BMPNC, BPMNC, BMMNC
      REAL*8      xx1, xx2
      REAL*8      XBSFM, XBSNC, XBNFM, XBNNC, XBVFM, XBVNC
      REAL*8      BETAFM, BETANC

      REAL*8, PARAMETER   ::  A = 1.246D0  ! approx Cunningham corr. factor
      
      REAL*8, PARAMETER   ::  TWOA = 2.0D0 * A

! *** Has a fixed number of Gauss-Herimite abscissas ( n )
      INTEGER, PARAMETER   ::  N = 5   ! one-half the number of abscissas
      REAL*8      GHXI(N)  ! Gauss-Hermite abscissas
      REAL*8      GHWI(N)  ! Gauss-Hermite weights

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

! adjoint variables

      REAL*8 KFM_ADJ   !INTENT OUT
      REAL   DG1_ADJ, XLNSIG1_ADJ, DG2_ADJ, XLNSIG2_ADJ   !INTENT OUT
      REAL*8 QUADS12_ADJ, QUADS21_ADJ, QUADN12_ADJ, QUADV12_ADJ   ! INTENT IN
      REAL*8 SUM1S12FM_ADJ, SUM1S21FM_ADJ, SUM2S12FM_ADJ, SUM2S21FM_ADJ
      REAL*8 SUM1NFM_ADJ, SUM2NFM_ADJ
      REAL*8 SUM1VFM_ADJ, SUM2VFM_ADJ
      REAL*8 SUM1S12NC_ADJ, SUM1S21NC_ADJ, SUM2S12NC_ADJ, SUM2S21NC_ADJ
      REAL*8 SUM1NNC_ADJ, SUM2NNC_ADJ, SUM1VNC_ADJ, SUM2VNC_ADJ
      REAL*8 XF_ADJ, DP1P_ADJ, DP1M_ADJ, DP1PSQ_ADJ, DP1MSQ_ADJ
      REAL*8 A1P_ADJ, A1M_ADJ, V1P_ADJ, V1M_ADJ
      REAL*8 A2P_ADJ, A2M_ADJ, V2P_ADJ, V2M_ADJ
      REAL*8 YF_ADJ, DP2P_ADJ ,DP2M_ADJ, DP2PSQ_ADJ, DP2MSQ_ADJ
      REAL*8 DSPP_ADJ, DSMP_ADJ, DSPM_ADJ, DSMM_ADJ
      REAL*8 BPPFM_ADJ, BMPFM_ADJ, BPMFM_ADJ, BMMFM_ADJ
      REAL*8 BPPNC_ADJ, BMPNC_ADJ, BPMNC_ADJ, BMMNC_ADJ
      REAL*8 XBSFM_ADJ, XBSNC_ADJ, XBNFM_ADJ, XBNNC_ADJ, XBVFM_ADJ
      REAL*8 XBVNC_ADJ

! *** The following expressions are from Binkowski & Shanker
!     Jour. Geophys. Research. Vol. 100,no. d12, pp 26,191-26,209
!     December 20, 1995

! ***  for Free Molecular Eq. A5
        BETAFM(xx1, xx2) = KFM *
     &       SQRT(1.D0 / xx1**3  + 1.D0 / xx2**3 ) * (xx1 + xx2)**2
     
! ***  for Near Continuum  Eq. A6
        BETANC(xx1, xx2) =  KNC * (xx1 + xx2) *
     &                       ( 1.0D0 / xx1 + 1.0D0 / xx2  +
     &                     TWOA * LAMDA * ( 1.0D0 / xx1 ** 2
     &                                    + 1.0D0 / xx2 **2 ) )

!---------------------------

! *** Initialize Adjoint Variables
      KFM_ADJ = 0.0
      DG1_ADJ = 0.0
      DG2_ADJ = 0.0
      XLNSIG1_ADJ = 0.0
      XLNSIG2_ADJ = 0.0
      SUM1S12FM_ADJ = 0.0
      SUM1S21FM_ADJ = 0.0
      SUM2S12FM_ADJ = 0.0
      SUM2S21FM_ADJ = 0.0
      SUM1NFM_ADJ = 0.0
      SUM2NFM_ADJ = 0.0
      SUM1VFM_ADJ = 0.0
      SUM2VFM_ADJ = 0.0
      SUM1S12NC_ADJ = 0.0
      SUM1S21NC_ADJ = 0.0
      SUM2S12NC_ADJ = 0.0
      SUM2S21NC_ADJ = 0.0
      SUM1NNC_ADJ = 0.0
      SUM2NNC_ADJ = 0.0
      SUM1VNC_ADJ = 0.0
      SUM2VNC_ADJ = 0.0
      XF_ADJ = 0.0
      DP1P_ADJ = 0.0
      DP1M_ADJ = 0.0
      DP1PSQ_ADJ = 0.0
      DP1MSQ_ADJ = 0.0
      A1P_ADJ = 0.0
      A1M_ADJ = 0.0
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
      BMPNC_ADJ = 0.0
      BPMNC_ADJ = 0.0
      BMMNC_ADJ = 0.0
      XBSFM_ADJ = 0.0
      XBSNC_ADJ = 0.0
      XBNFM_ADJ = 0.0
      XBNNC_ADJ = 0.0
      XBVFM_ADJ = 0.0
      XBVNC_ADJ = 0.0

! Forward Code:

      sum1s12fm = 0.D0
      sum1s12nc = 0.D0
      sum1s21fm = 0.D0
      sum1s21nc = 0.D0
      sum1vnc = 0.D0
      sum1vfm = 0.D0
      sum1nfm = 0.D0
      sum1nnc = 0.D0
      do 1 i=1,n

        sum2s12fm = 0.D0
        sum2s12nc = 0.D0
        sum2s21fm = 0.D0
        sum2s21nc = 0.D0
        sum2nfm = 0.D0
        sum2nnc = 0.D0
        sum2vnc = 0.D0
        sum2vfm = 0.D0
        xi = ghxi(i)
        wxi = ghwi(i)
        xf = exp( sqrt2 * xi *xlnsig1)
        dp1p = dg1*xf
        dp1m = dg1/xf
        dp1psq = dp1p*dp1p
        dp1msq = dp1m*dp1m
        a1p = dp1psq
        a1m = dp1msq
        v1p = dp1p*dp1psq
        v1m = dp1m*dp1msq

      do 11 j=1,n
        yi  = ghxi(j)
        wyi = ghwi(j)
        yf = exp( sqrt2 * yi * xlnsig2)
        dp2p = dg2*yf
        dp2m = dg2/yf
        dp2psq = dp2p*dp2p
        dp2msq = dp2m*dp2m
        a2p  = dp2psq
        a2m  = dp2msq
        v2p  =  dp2p*dp2psq
        v2m  = dp2m*dp2msq
        dspp = (v1p+v2p)**two3rds - a2p
        dsmp = (v1m+v2p)**two3rds - a2p
        dspm = (v1p+v2m)**two3rds - a2m
        dsmm = (v1m+v2m)**two3rds - a2m

        bppfm = betafm(dp1p,dp2p)
        bmpfm = betafm(dp1m,dp2p)
        bpmfm = betafm(dp1p,dp2m)
        bmmfm = betafm(dp1m,dp2m)
        
        bppnc = betanc(dp1p,dp2p)
        bmpnc = betanc(dp1m,dp2p) 
        bpmnc = betanc(dp1p,dp2m)         
        bmmnc = betanc(dp1m,dp2m)

        sum2s12fm = sum2s12fm + wyi*(a1p * bppfm + a1p * bpmfm
     &               +   a1m * bmpfm + a1m * bmmfm )

        sum2s21fm = sum2s21fm + wyi*(dspp * bppfm + dspm * bpmfm
     &               +   dsmp * bmpfm + dsmm * bmmfm )


        sum2s12nc = sum2s12nc + wyi*(a1p * bppnc + a1p * bpmnc
     &               +   a1m * bmpnc + a1m * bmmnc )

        sum2s21nc = sum2s21nc + wyi*(dspp * bppnc + dspm * bpmnc
     &               +   dsmp * bmpnc + dsmm * bmmnc )

        sum2nfm = sum2nfm + wyi*(bppfm + bmpfm + bpmfm + bmmfm)

        sum2nnc = sum2nnc + wyi*(bppnc + bmpnc + bpmnc + bmmnc)

        sum2vfm = sum2vfm + wyi*(v1p*(bppfm + bpmfm) +
     &                           v1m*(bmpfm + bmmfm) )

        sum2vnc = sum2vnc + wyi*(v1p*(bppnc + bpmnc) +
     &                           v1m*(bmpnc + bmmnc) )

   11 continue

      sum1s12fm = sum1s12fm + wxi * sum2s12fm
      sum1s21fm = sum1s21fm + wxi * sum2s21fm
      sum1nfm   = sum1nfm + wxi * sum2nfm
      sum1vfm   = sum1vfm + wxi * sum2vfm

      sum1s12nc = sum1s12nc + wxi * sum2s12nc
      sum1s21nc = sum1s21nc + wxi * sum2s21nc
      sum1nnc   = sum1nnc + wxi * sum2nnc
      sum1vnc   = sum1vnc + wxi * sum2vnc

    1 continue

C *** Second moment intermodal coagulation coefficients

c FSB NOTE: the transfer of second moment is not symmetric.
c     See equations A3 & A4 of Binkowski & Shankar (1995)

c ***  to accumulation mode from Aitken mode

      xbsfm   = sum1s21fm  / pi
      xbsnc   = sum1s21nc  / pi

      quads21 =  xbsfm * xbsnc / ( xbsfm + xbsnc )

c *** from Aitken mode to accumulation mode

      xbsfm   = sum1s12fm  / pi
      xbsnc   = sum1s12nc  / pi

      quads12 =  xbsfm * xbsnc / ( xbsfm + xbsnc )

      xbnfm   = sum1nfm  / pi 
      xbnnc   = sum1nnc  / pi

      quadn12 =  xbnfm * xbnnc / ( xbnfm + xbnnc )
     
c *** quadn12 is the intermodal coagulation coefficient for number

        
       xbvfm = sum1vfm / pi
       xbvnc = sum1vnc / pi
        
       quadv12 = xbvfm * xbvnc / ( xbvfm + xbvnc )

c *** quadv12 is the intermodal coagulation coefficient for 3rd moment

! Adjoint Code:

      !------
      ! fwd code:
      ! quadv12 = xbvfm * xbvnc / ( xbvfm + xbvnc )
      ! adj code:
      XBVFM_ADJ = XBVFM_ADJ + ( XBVNC / ( XBVFM + XBVNC ) 
     &          - XBVFM * XBVNC / ( ( XBVFM + XBVNC ) ** 2.0 ) 
     &          ) * QUADV12_ADJ
      XBVNC_ADJ = XBVNC_ADJ + ( XBVFM / ( XBVFM + XBVNC ) 
     &          - XBVFM * XBVNC / ( ( XBVFM + XBVNC ) ** 2.0 ) 
     &          ) * QUADV12_ADJ
      QUADV12_ADJ = 0.0

      !------
      ! fwd code:
      ! xbvfm = sum1vfm / pi
      ! xbvnc = sum1vnc / pi
      ! adj code:
      SUM1VNC_ADJ = SUM1VNC_ADJ + ( 1.0 / PI ) * XBVNC_ADJ
      XBVNC_ADJ = 0.0
      SUM1VFM_ADJ = SUM1VFM_ADJ + ( 1.0 / PI ) * XBVFM_ADJ
      XBVFM_ADJ = 0.0

      !------
      ! fwd code:
      ! quadn12 =  xbnfm * xbnnc / ( xbnfm + xbnnc )
      ! adj code:
      XBNFM_ADJ = XBNFM_ADJ + ( XBNNC / ( XBNFM + XBNNC ) 
     &          - XBNFM * XBNNC / ( ( XBNFM + XBNNC ) ** 2.0 ) 
     &          ) * QUADN12_ADJ
      XBNNC_ADJ = XBNNC_ADJ + ( XBNFM / ( XBNFM + XBNNC ) 
     &          - XBNFM * XBNNC / ( ( XBNFM + XBNNC ) ** 2.0 )
     &          ) * QUADN12_ADJ
      QUADN12_ADJ = 0.0

      !------
      ! fwd code:
      ! xbnfm   = sum1nfm  / pi
      ! xbnnc   = sum1nnc  / pi
      ! adj code:
      SUM1NNC_ADJ = SUM1NNC_ADJ + ( 1.0 / PI ) * XBNNC_ADJ
      XBNNC_ADJ = 0.0
      SUM1NFM_ADJ = SUM1NFM_ADJ + ( 1.0 / PI ) * XBNFM_ADJ
      XBNFM_ADJ = 0.0

      !------
      ! fwd code:
      ! quads12 =  xbsfm * xbsnc / ( xbsfm + xbsnc )
      ! adj code:
      XBSFM_ADJ = XBSFM_ADJ + ( XBSNC / ( XBSFM + XBSNC ) 
     &          - XBSFM * XBSNC / ( ( XBSFM + XBSNC ) ** 2.0 ) 
     &          ) * QUADS12_ADJ
      XBSNC_ADJ = XBSNC_ADJ + ( XBSFM / ( XBSFM + XBSNC ) 
     &          - XBSFM * XBSNC / ( ( XBSFM + XBSNC ) ** 2.0 ) 
     &          ) * QUADS12_ADJ
      QUADS12_ADJ = 0.0

      !------
      ! fwd code:
      ! xbsfm   = sum1s12fm  / pi
      ! xbsnc   = sum1s12nc  / pi
      ! adj code:
      SUM1S12NC_ADJ = SUM1S12NC_ADJ + ( 1.0 / PI ) * XBSNC_ADJ
      XBSNC_ADJ = 0.0
      SUM1S12FM_ADJ = SUM1S12FM_ADJ + ( 1.0 / PI ) * XBSFM_ADJ
      XBSFM_ADJ = 0.0

      !------
      ! fwd code:
      ! quads21 =  xbsfm * xbsnc / ( xbsfm + xbsnc )
      ! adj code:
      XBSFM_ADJ = XBSFM_ADJ + ( XBSNC / ( XBSFM + XBSNC ) 
     &          - XBSFM * XBSNC / ( ( XBSFM + XBSNC ) ** 2.0 ) 
     &          ) * QUADS21_ADJ
      XBSNC_ADJ = XBSNC_ADJ + ( XBSFM / ( XBSFM + XBSNC ) 
     &          - XBSFM * XBSNC / ( ( XBSFM + XBSNC ) ** 2.0 ) 
     &          ) * QUADS21_ADJ
      QUADS21_ADJ = 0.0

      !------
      ! fwd code:
      ! xbsfm   = sum1s21fm  / pi
      ! xbsnc   = sum1s21nc  / pi
      ! adj code:
      SUM1S21NC_ADJ = SUM1S21NC_ADJ + ( 1.0 / PI ) * XBSNC_ADJ
      XBSNC_ADJ = 0.0
      SUM1S21FM_ADJ = SUM1S21FM_ADJ + ( 1.0 / PI ) * XBSFM_ADJ 
      XBSFM_ADJ = 0.0

      DO I = N, 1, -1
   
         XI = GHXI(I)
         WXI = GHWI(I)
         XF = EXP( SQRT2 * XI * XLNSIG1 )
         DP1P = DG1*XF
         DP1M = DG1/XF
         DP1PSQ = DP1P*DP1P
         DP1MSQ = DP1M*DP1M
         A1P = DP1PSQ
         A1M = DP1MSQ
         V1P = DP1P*DP1PSQ
         V1M = DP1M*DP1MSQ

         !------
         ! fwd code:
         ! sum1s12nc = sum1s12nc + wxi * sum2s12nc
         ! sum1s21nc = sum1s21nc + wxi * sum2s21nc
         ! sum1nnc   = sum1nnc + wxi * sum2nnc
         ! sum1vnc   = sum1vnc + wxi * sum2vnc
         ! adj code:
         SUM2VNC_ADJ = SUM2VNC_ADJ + WXI * SUM1VNC_ADJ
         SUM2NNC_ADJ = SUM2NNC_ADJ + WXI * SUM1NNC_ADJ
         SUM2S21NC_ADJ = SUM2S21NC_ADJ + WXI * SUM1S21NC_ADJ
         SUM2S12NC_ADJ = SUM2S12NC_ADJ + WXI * SUM1S12NC_ADJ

         !------
         ! fwd code:
         ! sum1s12fm = sum1s12fm + wxi * sum2s12fm
         ! sum1s21fm = sum1s21fm + wxi * sum2s21fm
         ! sum1nfm   = sum1nfm + wxi * sum2nfm
         ! sum1vfm   = sum1vfm + wxi * sum2vfm
         ! adj code:
         SUM2VFM_ADJ = SUM2VFM_ADJ + WXI * SUM1VFM_ADJ
         SUM2NFM_ADJ = SUM2NFM_ADJ + WXI * SUM1NFM_ADJ
         SUM2S21FM_ADJ = SUM2S21FM_ADJ + WXI * SUM1S21FM_ADJ
         SUM2S12FM_ADJ = SUM2S12FM_ADJ + WXI * SUM1S12FM_ADJ

         DO J = N, 1, -1
  
            YI  = GHXI(J)
            WYI = GHWI(J)
            YF = EXP( SQRT2 * YI * XLNSIG2 )
            DP2P = DG2*YF
            DP2M = DG2/YF
            DP2PSQ = DP2P*DP2P
            DP2MSQ = DP2M*DP2M
            A2P  = DP2PSQ
            A2M  = DP2MSQ
            V2P  =  DP2P*DP2PSQ
            V2M  = DP2M*DP2MSQ
            DSPP = (V1P+V2P)**TWO3RDS - A2P
            DSMP = (V1M+V2P)**TWO3RDS - A2P
            DSPM = (V1P+V2M)**TWO3RDS - A2M
            DSMM = (V1M+V2M)**TWO3RDS - A2M
      
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
            ! sum2vnc = sum2vnc + wyi*(v1p*(bppnc + bpmnc) +
            ! &                        v1m*(bmpnc + bmmnc) )
            ! adj code:
            V1P_ADJ = V1P_ADJ + WYI * ( BPMNC + BPPNC ) * SUM2VNC_ADJ
            V1M_ADJ = V1M_ADJ + WYI * ( BMPNC + BMMNC ) * SUM2VNC_ADJ
            BPPNC_ADJ = BPPNC_ADJ + V1P * WYI * SUM2VNC_ADJ
            BPMNC_ADJ = BPMNC_ADJ + V1P * WYI * SUM2VNC_ADJ
            BMPNC_ADJ = BMPNC_ADJ + V1M * WYI * SUM2VNC_ADJ
            BMMNC_ADJ = BMMNC_ADJ + V1M * WYI * SUM2VNC_ADJ

            !------
            ! fwd code:
            ! sum2vfm = sum2vfm + wyi*(v1p*(bppfm + bpmfm) +
            ! &                        v1m*(bmpfm + bmmfm) )
            ! adj code:
            V1P_ADJ = V1P_ADJ + WYI * ( BPPFM + BPMFM ) * SUM2VFM_ADJ
            V1M_ADJ = V1M_ADJ + WYI * ( BMPFM + BMMFM ) * SUM2VFM_ADJ
            BPPFM_ADJ = BPPFM_ADJ + V1P * WYI * SUM2VFM_ADJ
            BPMFM_ADJ = BPMFM_ADJ + V1P * WYI * SUM2VFM_ADJ
            BMPFM_ADJ = BMPFM_ADJ + V1M * WYI * SUM2VFM_ADJ
            BMMFM_ADJ = BMMFM_ADJ + V1M * WYI * SUM2VFM_ADJ

            !------
            ! fwd code:
            ! sum2nfm = sum2nfm + wyi*(bppfm + bmpfm + bpmfm + bmmfm)
            ! sum2nnc = sum2nnc + wyi*(bppnc + bmpnc + bpmnc + bmmnc)
            ! adj code:
            BPPNC_ADJ = BPPNC_ADJ + WYI * SUM2NNC_ADJ
            BMPNC_ADJ = BMPNC_ADJ + WYI * SUM2NNC_ADJ
            BPMNC_ADJ = BPMNC_ADJ + WYI * SUM2NNC_ADJ
            BMMNC_ADJ = BMMNC_ADJ + WYI * SUM2NNC_ADJ
            BPPFM_ADJ = BPPFM_ADJ + WYI * SUM2NFM_ADJ
            BMPFM_ADJ = BMPFM_ADJ + WYI * SUM2NFM_ADJ
            BPMFM_ADJ = BPMFM_ADJ + WYI * SUM2NFM_ADJ
            BMMFM_ADJ = BMMFM_ADJ + WYI * SUM2NFM_ADJ

            !------
            ! fwd code:
            ! sum2s21nc = sum2s21nc + wyi*(dspp * bppnc + dspm * bpmnc
            ! &                        +   dsmp * bmpnc + dsmm * bmmnc )
            ! adj code:
            DSPP_ADJ = DSPP_ADJ + WYI * BPPNC * SUM2S21NC_ADJ
            BPPNC_ADJ = BPPNC_ADJ + WYI * DSPP * SUM2S21NC_ADJ
            DSPM_ADJ = DSPM_ADJ + WYI * BPMNC * SUM2S21NC_ADJ
            BPMNC_ADJ = BPMNC_ADJ + WYI * DSPM * SUM2S21NC_ADJ
            DSMP_ADJ = DSMP_ADJ + WYI * BMPNC * SUM2S21NC_ADJ
            BMPNC_ADJ = BMPNC_ADJ + WYI * DSMP * SUM2S21NC_ADJ
            DSMM_ADJ = DSMM_ADJ + WYI * BMMNC * SUM2S21NC_ADJ
            BMMNC_ADJ = BMMNC_ADJ + WYI * DSMM * SUM2S21NC_ADJ

            !------
            ! fwd code:
            ! sum2s12nc = sum2s12nc + wyi*(a1p * bppnc + a1p * bpmnc
            ! &               +   a1m * bmpnc + a1m * bmmnc ) 
            ! adj code:
            BPPNC_ADJ = BPPNC_ADJ + WYI * A1P * SUM2S12NC_ADJ
            BPMNC_ADJ = BPMNC_ADJ + WYI * A1P * SUM2S12NC_ADJ
            BMPNC_ADJ = BMPNC_ADJ + WYI * A1M * SUM2S12NC_ADJ
            BMMNC_ADJ = BMMNC_ADJ + WYI * A1M * SUM2S12NC_ADJ
            A1P_ADJ = A1P_ADJ + WYI * ( BPMNC + BPPNC ) * SUM2S12NC_ADJ
            A1M_ADJ = A1M_ADJ + WYI * ( BMPNC + BMMNC ) * SUM2S12NC_ADJ

            !------
            ! fwd code:
            ! sum2s21fm = sum2s21fm + wyi*(dspp * bppfm + dspm * bpmfm
            ! &               +   dsmp * bmpfm + dsmm * bmmfm )
            ! adj code:
            DSPP_ADJ = DSPP_ADJ + WYI * BPPFM * SUM2S21FM_ADJ
            BPPFM_ADJ = BPPFM_ADJ + WYI * DSPP * SUM2S21FM_ADJ
            DSPM_ADJ = DSPM_ADJ + WYI * BPMFM * SUM2S21FM_ADJ
            BPMFM_ADJ = BPMFM_ADJ + WYI * DSPM * SUM2S21FM_ADJ
            DSMP_ADJ = DSMP_ADJ + WYI * BMPFM * SUM2S21FM_ADJ
            BMPFM_ADJ = BMPFM_ADJ + WYI * DSMP * SUM2S21FM_ADJ
            DSMM_ADJ = DSMM_ADJ + WYI * BMMFM * SUM2S21FM_ADJ
            BMMFM_ADJ = BMMFM_ADJ + WYI * DSMM * SUM2S21FM_ADJ

            !------
            ! fwd code:
            ! sum2s12fm = sum2s12fm + wyi*(a1p * bppfm + a1p * bpmfm
            ! &               +   a1m * bmpfm + a1m * bmmfm )
            ! adj code:
            BPPFM_ADJ = BPPFM_ADJ + WYI * A1P * SUM2S12FM_ADJ 
            BPMFM_ADJ = BPMFM_ADJ + WYI * A1P * SUM2S12FM_ADJ
            BMPFM_ADJ = BMPFM_ADJ + WYI * A1M * SUM2S12FM_ADJ
            BMMFM_ADJ = BMMFM_ADJ + WYI * A1M * SUM2S12FM_ADJ
            A1P_ADJ = A1P_ADJ + WYI * ( BPPFM + BPMFM ) * SUM2S12FM_ADJ
            A1M_ADJ = A1M_ADJ + WYI * ( BMPFM + BMMFM ) * SUM2S12FM_ADJ

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
     &                ( 1.0 / DP2P ) + ( ( 1.0 / ( DP1M ** 2.0 ) )  +
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
           ! dsmm = (v1m+v2m)**two3rds - a2m
           ! adj code:
           A2M_ADJ = A2M_ADJ - DSMM_ADJ
           V1M_ADJ = V1M_ADJ + TWO3RDS * ( ( V1M + V2M ) ** ( 
     ^               TWO3RDS - 1.0 ) ) * DSMM_ADJ
           V2M_ADJ = V2M_ADJ + TWO3RDS * ( ( V1M + V2M ) ** ( 
     ^               TWO3RDS - 1.0 ) ) * DSMM_ADJ
           DSMM_ADJ = 0.0

           !------ 
           ! fwd code:
           ! dspm = 0.5D0*(v1p+v2m)**two3rds - a2m
           ! adj code:
           A2M_ADJ = A2M_ADJ - DSPM_ADJ
           V1P_ADJ = V1P_ADJ + TWO3RDS * ( ( V1P + V2M ) ** (
     &               TWO3RDS - 1.0 ) ) * DSPM_ADJ
           V2M_ADJ = V2M_ADJ + TWO3RDS * ( ( V1P + V2M ) ** (
     &               TWO3RDS - 1.0 ) ) * DSPM_ADJ
           DSPM_ADJ = 0.0

           !------
           ! fwd code:
           ! dsmp = 0.5D0*(v1m+v2p)**two3rds - a2p
           ! adj code:
           A2P_ADJ = A2P_ADJ - DSMP_ADJ
           V1M_ADJ = V1M_ADJ + TWO3RDS * ( ( V1M + V2P ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSMP_ADJ
           V2P_ADJ = V2P_ADJ + TWO3RDS * ( ( V1M + V2P ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSMP_ADJ
           DSMP_ADJ = 0.0

           !------
           ! fwd code:
           ! dspp = 0.5D0*(v1p+v2p)**two3rds - a2p
           ! adj code:
           A2P_ADJ = A2P_ADJ - DSPP_ADJ
           V1P_ADJ = V1P_ADJ + TWO3RDS * ( ( V1P + V2P ) ** ( 
     &               TWO3RDS - 1.0 ) ) * DSPP_ADJ
           V2P_ADJ = V2P_ADJ + TWO3RDS * ( ( V1P + V2P ) ** ( 
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
           ! dp2p = dg2*yf
           ! dp2m = dg2/yf
           ! adj code:
           DG2_ADJ = DG2_ADJ + ( 1 / YF ) * DP2M_ADJ
           YF_ADJ = YF_ADJ - ( DG2 / ( YF ** 2.0 ) ) * DP2M_ADJ
           DP2M_ADJ = 0.0
           DG2_ADJ = DG2_ADJ + YF * DP2P_ADJ
           YF_ADJ = YF_ADJ + DG2 * DP2P_ADJ
           DP2P_ADJ = 0.0

           !------
           ! fwd code:
           ! yf = exp( sqrt2 * yi * xlnsig2)
           ! adj code:
           XLNSIG2_ADJ = XLNSIG2_ADJ + SQRT2 * YI * EXP( SQRT2 * YI * 
     &                   XLNSIG2 ) * YF_ADJ
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
          ! a1p = dp1psq
          ! a1m = dp1msq
          ! adj code:
          DP1MSQ_ADJ = DP1MSQ_ADJ + A1M_ADJ
          A1M_ADJ = 0.0
          DP1PSQ_ADJ = DP1PSQ_ADJ + A1P_ADJ 
          A1P_ADJ = 0.0
  
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
          ! dp1p = dg1*xf
          ! dp1m = dg1/xf
          ! adj code:
          DG1_ADJ = DG1_ADJ + ( 1.0 / XF ) * DP1M_ADJ
          XF_ADJ = XF_ADJ - ( DG1 / ( XF ** 2.0 ) ) * DP1M_ADJ
          DP1M_ADJ = 0.0
          DG1_ADJ = DG1_ADJ + XF * DP1P_ADJ
          XF_ADJ = XF_ADJ + DG1 * DP1P_ADJ
          DP1P_ADJ = 0.0
 
          !------
          ! fwd code:
          ! xf = exp( sqrt2 * xi *xlnsig1)
          ! adj code:
          XLNSIG1_ADJ = XLNSIG1_ADJ + SQRT2 * XI * EXP( SQRT2 * XI * 
     &                XLNSIG1 ) * XF_ADJ
          XF_ADJ = 0.0
  
        END DO

      END SUBROUTINE INTERCOAG_GH_ADJ
