C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP1F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY 
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISRP1F (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION WI, RHI, TEMPI
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
      DOUBLE PRECISION DC
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = (W(3))/(W(2))
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (2.0.LE.SULRAT) THEN 
C         WP = W
         CALL ISRP1FA
C         DC   = W(3) - 2.001D0*W(2)  ! For numerical stability
C         W(3) = W(3) + MAX(-DC, ZERO)
CC
CC      IF(METSTBL.EQ.1) THEN
C         SCASE = 'A2'
C         CALL CALCA2                 ! Only liquid (metastable)
C      ELSE
C
C         IF (RH.LT.DRNH42S4) THEN    
C            SCASE = 'A1'
C            CALL CALCA1              ! NH42SO4              ; case A1
CC
C         ELSEIF (DRNH42S4.LE.RH) THEN
C            SCASE = 'A2'
C            CALL CALCA2              ! Only liquid          ; case A2
C         ENDIF
C      ENDIF
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'B4'
C         WP = W
         CALL CALCB4                 ! Only liquid (metastable)
C      ELSE
CC
C         IF (RH.LT.DRNH4HS4) THEN         
C            SCASE = 'B1'
C            CALL CALCB1              ! NH4HSO4,LC,NH42SO4   ; case B1
CC
C         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
C            SCASE = 'B2'
C            CALL CALCB2              ! LC,NH42S4            ; case B2
CC
C         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
C            SCASE = 'B3'
C            CALL CALCB3              ! NH42S4               ; case B3
CC
C         ELSEIF (DRNH42S4.LE.RH) THEN         
C            SCASE = 'B4'
C            CALL CALCB4              ! Only liquid          ; case B4
C         ENDIF
C      ENDIF
      CALL CALCACT3F              !
      CALL CALCNH3
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'C2'
         CALL CALCC2                 ! Only liquid (metastable)
C      ELSE
CC
C         IF (RH.LT.DRNH4HS4) THEN         
C            SCASE = 'C1'
C            CALL CALCC1              ! NH4HSO4              ; case C1
CC
C         ELSEIF (DRNH4HS4.LE.RH) THEN         
C            SCASE = 'C2'
C            CALL CALCC2              ! Only liquid          ; case C2
CC
C         ENDIF
C      ENDIF
      CALL CALCACT3F              !
      CALL CALCNH3
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE ISRP1F *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP1F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY 
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISRP1FA 
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION WI, RHI, TEMPI
      DOUBLE PRECISION WP(NCOMP), MOLALP(NIONS)
      DOUBLE PRECISION DC, GAS(3), AERLIQ(7)
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
C      CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
C      SULRAT = (W(3))/(W(2))
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
C      IF (2.0.LE.SULRAT) THEN 
C         W = WP
         DC   = W(3) - 2.001D0*W(2)  ! For numerical stability
         W(3) = W(3) + MAX(-DC, ZERO)
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'A2'
         CALL CALCA2                 ! Only liquid (metastable)
C         GAS(1) = GNH3                ! Gaseous aerosol species
C         GAS(2) = GHNO3
C         GAS(3) = GHCL
CC
C         DO I=1,NIONS              ! Liquid aerosol species
C            AERLIQ(I) = MOLAL(I)
C         ENDDO 
C      ELSE
C
C         IF (RH.LT.DRNH42S4) THEN    
C            SCASE = 'A1'
C            CALL CALCA1              ! NH42SO4              ; case A1
CC
C         ELSEIF (DRNH42S4.LE.RH) THEN
C            SCASE = 'A2'
C            CALL CALCA2              ! Only liquid          ; case A2
C         ENDIF
C      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE ISRP1F *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP2F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
C     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISRP2F (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION WI, RHI, TEMPI
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = (W(3))/(W(2))
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (2.0.LE.SULRAT) THEN                
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'D3'
         CALL CALCD3                   ! Only liquid (metastable)
C      ELSE
CC
C         IF (RH.LT.DRNH4NO3) THEN    
C            SCASE = 'D1'
C            CALL CALCD1              ! NH42SO4,NH4NO3       ; case D1
CC
C         ELSEIF (DRNH4NO3.LE.RH .AND. RH.LT.DRNH42S4) THEN         
C            SCASE = 'D2'
C            CALL CALCD2              ! NH42S4               ; case D2
CC
C         ELSEIF (DRNH42S4.LE.RH) THEN
C            SCASE = 'D3'
C            CALL CALCD3              ! Only liquid          ; case D3
C         ENDIF
C      ENDIF
C
C *** SULFATE RICH (NO ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
C     SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
C      IF(METSTBL.EQ.1) THEN
C         WP = W
         SCASE = 'E4'
         CALL CALCB4E                    ! Only liquid (metastable)
         SCASE = 'E4'
C      ELSE
CC
C         IF (RH.LT.DRNH4HS4) THEN         
C            SCASE = 'B1'
C            CALL CALCB1              ! NH4HSO4,LC,NH42SO4   ; case E1
C            SCASE = 'E1'
CC
C         ELSEIF (DRNH4HS4.LE.RH .AND. RH.LT.DRLC) THEN         
C            SCASE = 'B2'
C            CALL CALCB2              ! LC,NH42S4            ; case E2
C            SCASE = 'E2'
CC
C         ELSEIF (DRLC.LE.RH .AND. RH.LT.DRNH42S4) THEN         
C            SCASE = 'B3'
C            CALL CALCB3              ! NH42S4               ; case E3
C            SCASE = 'E3'
CC
C         ELSEIF (DRNH42S4.LE.RH) THEN         
C            SCASE = 'B4'
C            CALL CALCB4              ! Only liquid          ; case E4
C            SCASE = 'E4'
C         ENDIF
C      ENDIF
C
      CALL CALCACT3F              !
      CALL CALCNA                 ! HNO3(g) DISSOLUTION
C
C *** SULFATE RICH (FREE ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
C     SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
C      IF(METSTBL.EQ.1) THEN
C         WP = W
         SCASE = 'F2'
         CALL CALCC2F                  ! Only liquid (metastable)
         SCASE = 'F2'
C      ELSE
CC
C         IF (RH.LT.DRNH4HS4) THEN         
C            SCASE = 'C1'
C            CALL CALCC1              ! NH4HSO4              ; case F1
C            SCASE = 'F1'
CC
C         ELSEIF (DRNH4HS4.LE.RH) THEN         
C            SCASE = 'C2'
C            CALL CALCC2              ! Only liquid          ; case F2
C            SCASE = 'F2'
C         ENDIF
C      ENDIF
C
      CALL CALCACT3F              !
      CALL CALCNA                 ! HNO3(g) DISSOLUTION
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE ISRP2F *****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP3F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM 
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISRP3F (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DIMENSION WI(NCOMP)
C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
C
      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
      ENDIF
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL ISOINIT3 (WI, RHI, TEMPI)
C
C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
C
      REST = 2.D0*W(2) + W(4) + W(5) 
      IF (W(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
         W(1) = (ONE-1D-6)*REST         ! Adjust Na amount
         CALL PUSHERR (0050, 'ISRP3F')  ! Warning error: Na adjusted
      ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
      SULRAT = (W(1)+W(3))/W(2)
      SODRAT = W(1)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

C *** SULFATE POOR ; SODIUM POOR
C
      IF (2.0.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'G5'
         CALL CALCG5                 ! Only liquid (metastable)
C      ENDIF
C
C *** SULFATE POOR ; SODIUM RICH
C
      ELSE IF (SULRAT.GE.2.0 .AND. SODRAT.GE.2.0) THEN                
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'H6'
         CALL CALCH6                 ! Only liquid (metastable)
C      ENDIF
C
C *** SULFATE RICH (NO ACID) 
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'I6'
         CALL CALCI6                 ! Only liquid (metastable)
C      ENDIF
C                                    
      CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL CALCACT3F              !
      CALL CALCNH3                !                NH3 
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
C      IF(METSTBL.EQ.1) THEN
         SCASE = 'J3'
         CALL CALCJ3                 ! Only liquid (metastable)
C      ENDIF
C                                    
      CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL CALCACT3F              !
      CALL CALCNH3                !                NH3 
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE ISRP3F *****************************************
C
      END
C
C
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCA2
C *** CASE A2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A !!!!!BISECTION IS PERFORMED TOWARDS X, THE
C     AMOUNT OF HYDROGEN IONS (H+) FOUND IN THE LIQUID PHASE.
C     FOR EACH ESTIMATION OF H+, FUNCTION FUNCB2A CALCULATES THE
C     CONCENTRATION OF IONS FROM THE NH3(GAS) - NH4+(LIQ) EQUILIBRIUM.
C     ELECTRONEUTRALITY IS USED AS THE OBJECTIVE FUNCTION.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCA2
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION DELTA
C
C *** CREATE ACTIVE VARIABLES FOR DIFFERENTIATION *********************
C
C      W(2) = XiSO4
C      W(3) = XiNH3
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU    =.TRUE.              ! Outer loop activity calculation flag
C      UCONLO    = TINY               ! Low  limit: No excess NH3 dissolves
C      UCONHI    = W(3) - 2.0D0*W(2)  ! High limit: All NH3 remaining in gas dissolves
C
C *** CALCULATE WATER CONTENT *****************************************
C
      MOLAL(5) = W(2)
      MOLAL(6) = ZERO
C
C      CALL CALCMR
C
         MOLALR(4) = MOLAL(5)+MOLAL(6) ! (NH4)2SO4 - CORRECT FOR SO4 TO HSO4
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
      WATER = ZERO
      DO I=1,NPAIR
         WATER = WATER + MOLALR(I)/M0(I)
      ENDDO
      WATER = MAX(WATER, TINY)
C
C *** CREATE ITERATION FOR ACTIVITY COEFFICIENTS
C
      CALL FUNCA2P
C
      IF ((MOLAL(1)).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCA2H ****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCA2P
C *** CASE A2 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE A2 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA2P.
C
C=======================================================================
C
      SUBROUTINE FUNCA2P
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION LAMDA, DISC, SQDR, THRSHHI, THRSHLO
      DOUBLE PRECISION NCON, QCON, UCON, UCONOLD
      DOUBLE PRECISION :: W2, W3
      LOGICAL TST, TST2
      INTEGER I
C
C *** SETUP PARAMETERS ************************************************
C
      TST    = .TRUE.
      TST2   = .TRUE.
      FRST   = .TRUE.
      CALAIN = .TRUE.
      W2    = W(2)         ! INITIAL AMOUNT OF (NH4)2SO4 IN SOLUTION
      W3    = W(3)
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      I=1
      UCON = 0.D0
      DO WHILE ((I.LE.14).AND. TST .AND. TST2)
C      DO I=1,14 !NSWEEP
         UCONOLD = UCON
         A2    = XK2*R*TEMP/XKW*(GAMA(8)/GAMA(9))**2.
C
         AA   = -A2
         BB   = A2*W3 - 2.D0*A2*W2 + 1
         CC   = 2.d0*W2
         DISC = BB*BB - 4.D0*AA*CC
         SQDR = SQRT(DISC)
C
         RT1  = (-BB + SQDR)/2.D0/AA
         RT2  = (-BB - SQDR)/2.D0/AA
C         WRITE(*,*) 'ROOTS', RT1, RT2
C
         IF ((RT1).LT.ZERO .AND. (RT2).GE.ZERO) THEN
            UCON = RT1
         ELSEIF ((RT2).LT.ZERO .AND. (RT1).GE.ZERO) THEN
            UCON = RT2
         ELSE
             TST2 = .FALSE.
         ENDIF
C
         QCON = -UCON
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = QCON                        ! HI
C         MOLAL (3) = MAX(W(3)/(ONE/A2/OMEGI + ONE), 2.*MOLAL(5)) ! NH4I
C         MOLAL (3) = MAX(2.D0*W2 + UCON, TINY)   ! NH4I
         IF (TINY .GT. (2.D0*W2 + UCON)) THEN
            MOLAL(3) = TINY
         ELSE
            MOLAL(3) = 2.D0*W2 + UCON
         ENDIF
         MOLAL (5) = W2                          ! SO4I
         MOLAL (6) = ZERO                        ! HSO4I
C         GNH3      = MAX(W(3)-MOLAL(3), TINY)   ! NH3GI
         IF (TINY .GT. (W(3)-MOLAL(3))) THEN
            GNH3 = TINY
         ELSE
            GNH3 = W(3)-MOLAL(3)
         ENDIF
         COH       = XKW/QCON                    ! OHI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         THRSHLO = UCONOLD - UCONOLD*1.0D-15
         THRSHHI = UCONOLD + UCONOLD*1.0D-15 
         IF (((UCON).LE.(THRSHLO)).AND.
     &       ((UCON).GE.(THRSHHI))) THEN
          	TST = .FALSE.
            CALL CALCACT3F
         ELSE 
            TST = .TRUE.
            CALL CALCACT3P     
         ENDIF
C         CALL CALCACT     
C         ELSE
C            GOTO 20
C         ENDIF
         I = I + 1
CC
CC *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
CC
CC         WRITE(*,*) 'THRESHOLDS', THRSHHI, THRSHLO
C         IF ((UCON.LE.THRSHLO).AND.(UCON.GE.THRSHHI)) THEN
CC            WRITE(*,*) 'At convergence, I:', I
C            RETURN
C         ELSE
C            CALL CALCACT     
C         ENDIF
      ENDDO
C10    CONTINUE
C
C *** END OF FUNCTION FUNCA2 ********************************************
C
      END SUBROUTINE FUNCA2P
C


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB4
C *** CASE B4 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
C     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
C     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCB4!(WP, GAS, AERLIQ)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION X,Y, SO4I, HSO4I, BB, CC, DD
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
      INTEGER          I
C
C *** SOLVE EQUATIONS **************************************************
C
C     W = WP
C
      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.
C
C *** CALCULATE WATER CONTENT ******************************************
C
C
C      CALL CALCB1A         ! GET DRY SALT CONTENT, AND USE FOR WATER.
C
C *** SETUP PARAMETERS ************************************************
C
      X = 2.d0*W(2)-W(3)       ! Equivalent NH4HSO4
      Y = W(3)-W(2)         ! Equivalent (NH4)2SO4
C
C *** CALCULATE COMPOSITION *******************************************
C
      IF ((X).LE.(Y)) THEN      ! LC is the MIN(x,y)
         CLC     = 2.D0*W(2)-W(3) !X        ! NH4HSO4 >= (NH4)2S04
         CNH4HS4 = ZERO
         CNH42S4 = 2.D0*W(3) - 3.D0*W(2) !Y-X
      ELSE
         CLC     = W(3)-W(2)      !Y        ! NH4HSO4 <  (NH4)2S04
         CNH4HS4 = 3.D0*W(2) - 2.D0*W(3)   !X-Y
         CNH42S4 = ZERO
      ENDIF
C
      MOLALR(13) = CLC       
      MOLALR(9)  = CNH4HS4   
      MOLALR(4)  = CNH42S4   
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)
C
      MOLAL(3)   = W(3)   ! NH4I
C
      I = 1
C      NSWEEP = 50
      DO WHILE ((I.LE.NSWEEP).AND.(CALAIN))
C         IF (I.GT.1) CALL CALCACT3
         AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
         BET   = W(2)
         GAM   = MOLAL(3)
C
         BB    = BET + AK1 - GAM
         CC    =-AK1*BET
         DD    = BB*BB - 4.D0*CC
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (5) = MAX(MIN(0.5*(-BB + SQRT(DD)), W(2)),TINY) ! SO4I
         MOLAL (6) = MAX(MIN(W(2)-MOLAL(5), W(2)), TINY)         ! HSO4I
         MOLAL (1) = MAX(MIN(AK1*MOLAL(6)/MOLAL(5), W(2)), TINY) ! HI
C
C         CALL CALCMR                                           ! Water content
C
         SO4I  = MOLAL(5)-MOLAL(1)     ! CORRECT FOR HSO4 DISSOCIATION 
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF ((SO4I).LT.(HSO4I)) THEN                
            MOLALR(13) = SO4I                   ! [LC] = [SO4]       
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)  ! NH4HSO4
         ELSE                                   
            MOLALR(13) = HSO4I                  ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)  ! (NH4)2SO4
         ENDIF
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
C         IF (.NOT.CALAIN) GOTO 30
         I = I + 1
         CALL CALCACT3  !*** slc.11.2009 moved to beginning of loop
C20    CONTINUE
      ENDDO
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCB4')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C      CALL CALCNH3
C
C         GAS(1) = GNH3                ! Gaseous aerosol species
C         GAS(2) = GHNO3
C         GAS(3) = GHCL
CC
C         DO I=1,NIONS              ! Liquid aerosol species
C            AERLIQ(I) = MOLAL(I)
C         ENDDO 
C30    RETURN
C
C *** END OF SUBROUTINE CALCB4 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB4E
C *** CASE B4 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
C     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
C     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCB4E!(WP, GAS, AERLIQ)   
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: X,Y, SO4I, HSO4I, BB, CC, DD
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
      INTEGER    :: I
C
C *** SOLVE EQUATIONS **************************************************
C
C      W = WP
C
      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.
C
C *** CALCULATE WATER CONTENT ******************************************
C
C      CALL CALCB1A         ! GET DRY SALT CONTENT, AND USE FOR WATER.
C
C *** SETUP PARAMETERS ************************************************
C
      X = 2.d0*W(2)-W(3)       ! Equivalent NH4HSO4
      Y = W(3)-W(2)         ! Equivalent (NH4)2SO4
C
C *** CALCULATE COMPOSITION *******************************************
C
      IF ((X).LE.(Y)) THEN       ! LC is the MIN(x,y)
         CLC     = 2.d0*W(2)-W(3)	 !X        ! NH4HSO4 >= (NH4)2S04
         CNH4HS4 = ZERO
         CNH42S4 = 2.d0*W(3)-3.d0*W(2) !Y-X
      ELSE
         CLC     = W(3)-W(2)     !Y        ! NH4HSO4 <  (NH4)2S04
         CNH4HS4 = 3.d0*W(2)-2.d0*W(3) !X-Y
         CNH42S4 = ZERO
      ENDIF
C
      MOLALR(13) = CLC       
      MOLALR(9)  = CNH4HS4   
      MOLALR(4)  = CNH42S4   
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)
C
      MOLAL(3)   = W(3)   ! NH4I
C
      I = 1
C      NSWEEP = 50
      DO WHILE ((I.LE.NSWEEP).AND.(CALAIN))
C         IF (I.GT.1) CALL CALCACT3
         AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
         BET   = W(2)
         GAM   = MOLAL(3)
C
         BB    = BET + AK1 - GAM
         CC    =-AK1*BET
         DD    = BB*BB - 4.D0*CC
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (5) = MAX(MIN(0.5d0*(-BB + SQRT(DD)), W(2)),TINY) ! SO4I
C         WRITE(*,*) 'MOLAL(5) ', molal(5), 'w2 ', w(2)
         MOLAL (6) = MAX(MIN(W(2)-MOLAL(5), W(2)), TINY)         ! HSO4I
C         WRITE(*,*) 'MOLAL(6) ', MOLAL(6)
         MOLAL (1) = MAX(MIN(AK1*MOLAL(6)/MOLAL(5), W(2)), TINY) ! HI
C         WRITE(*,*) 'MOLAL(5, 6) ', MOLAL(5), MOLAL(6)
C
C         CALL CALCMR                                           ! Water content
C      slc.1.2011 - calling CALCMR for case E rather than B
C
         SO4I  = MOLAL(5)-MOLAL(1)     ! CORRECT FOR HSO4 DISSOCIATION  as from B4
C         SO4I  = MAX(MOLAL(5)-MOLAL(1),ZERO)      ! FROM HSO4 DISSOCIATION 
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF ((SO4I).LT.(HSO4I)) THEN                
            MOLALR(13) = SO4I                     ! [LC] = [SO4] 
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)    ! NH4HSO4
         ELSE                                   
            MOLALR(13) = HSO4I                    ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)    ! (NH4)2SO4
         ENDIF
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
C         IF (.NOT.CALAIN) GOTO 30
         I = I + 1
         CALL CALCACT3  !*** slc.11.2009 moved to beginning of loop
C20    CONTINUE
C         WRITE(*,*) 'Inside E4 loop, i: ',i
C         WRITE(*,*) ' MOLAL',MOLAL
C         WRITE(*,*) 'gama ', gama(4), gama(5), gama(7), gama(8)
C         WRITE(*,*) 'gama ', gama(9), gama(10), gama(13)
C         WRITE(*,*) 'water ',WATER
C         PAUSE
      ENDDO
C      WRITE(*,*) 'W',W,'RH',RH,'TEMP',TEMP
C      WRITE(*,*) ' MOLAL',MOLAL
C      WRITE(*,*) 'gama ', gama(4), gama(5), gama(7), gama(8)
C      WRITE(*,*) 'gama ', gama(9), gama(10), gama(13)
C      WRITE(*,*) 'water ',WATER
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCB4E')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C
C      CALL CALCNA
C
C         GAS(1) = GNH3                ! Gaseous aerosol species
C         GAS(2) = GHNO3
C         GAS(3) = GHCL
CC
C         DO I=1,NIONS              ! Liquid aerosol species
C            AERLIQ(I) = MOLAL(I)
C         ENDDO 
C
C30    RETURN
C
C *** END OF SUBROUTINE CALCB4 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCC2
C *** CASE C2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCC2!(WP, GAS, AERLIQ)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION LAMDA, KAPA, PSI, PARM
      DOUBLE PRECISION BB, CC
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
      INTEGER          I
C
C      W = WP
C
      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** SOLVE EQUATIONS **************************************************
C
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION
      I = 1
C      NSWEEP = 50
      DO WHILE ((I.LE.NSWEEP).AND.(CALAIN))
C         IF (I.GT.1) CALL CALCACT3 
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         BB    = PSI+PARM
         CC    =-PARM*(LAMDA+PSI)
         KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL(1) = PSI+KAPA                               ! HI
         MOLAL(3) = LAMDA                                  ! NH4I
         MOLAL(5) = KAPA                                   ! SO4I
         MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
         CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  ! Free H2SO4
C
C         CALL CALCMR                                       ! Water content
C
         MOLALR(9) = MOLAL(3)                     ! NH4HSO4  *** As in ISORROPIA 1.7
         MOLALR(7) = MAX(W(2)-W(3), ZERO)         ! H2SO4
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
C         IF (.NOT.CALAIN) GOTO 30
         I = I + 1
         CALL CALCACT3   !*** slc.11.2009 Moved to beginning of loop     
C20    CONTINUE
      ENDDO
C      WRITE(*,*) 'I',I,'MOLAL(1,5,6)',MOLAL(1),MOLAL(5),MOLAL(6)
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCC2')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C      CALL CALCNH3
C
C         GAS(1) = GNH3                ! Gaseous aerosol species
C         GAS(2) = GHNO3
C         GAS(3) = GHCL
CC
C         DO I=1,NIONS              ! Liquid aerosol species
C            AERLIQ(I) = MOLAL(I)
C         ENDDO 
C 
C30    RETURN
C    
C *** END OF SUBROUTINE CALCC2 *****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCC2
C *** CASE C2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCC2F!(WP, GAS, AERLIQ)  
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: LAMDA, KAPA, PSI, PARM
      DOUBLE PRECISION :: BB, CC
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
      INTEGER :: I
C
C      W = WP
C
      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** SOLVE EQUATIONS **************************************************
C
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION
      I = 1
C      NSWEEP = 50
      DO WHILE ((I.LE.NSWEEP).AND.(CALAIN))
C         IF (I.GT.1) CALL CALCACT3 
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         BB    = PSI+PARM
         CC    =-PARM*(LAMDA+PSI)
         KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL(1) = PSI+KAPA                               ! HI
         MOLAL(3) = LAMDA                                  ! NH4I
         MOLAL(5) = KAPA                                   ! SO4I
         MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
         CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  ! Free H2SO4
C
C         CALL CALCMR                                       ! Water content
C
C      slc.1.2011 - calling CALCMR for case F rather than C
C
         MOLALR(9) = MOLAL(3)                    ! NH4HSO4 - slc.1.2011 - from ISORROPIA 1.7
         MOLALR(7) = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3),ZERO)  ! H2SO4
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C         WRITE(*,*) 'Iteration: i', I
C         WRITE(*,*) 'MOLAL ',MOLAL(1), MOLAL(3), MoLAL(5), MOLAL(6)
C         WRITE(*,*) 'MOLALR ', (MOLALR(7)), (MOLALR(4)) 
C         WRITE(*,*) 'M0 ',(M0(7)), (M0(4))
C         WRITE(*,*) 'GAMA ', (GAMA)
C         WRITE(*,*) 'water', water
C         PAUSE
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
C         IF (.NOT.CALAIN) GOTO 30
         I = I + 1
         CALL CALCACT3   !*** slc.11.2009 Moved to beginning of loop     
C20    CONTINUE
      ENDDO
C      WRITE(*,*) 'I',I,'MOLAL(1,5,6)',MOLAL(1),MOLAL(5),MOLAL(6)
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCC2F')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C
C      CALL CALCNA
C
C         GAS(1) = GNH3                ! Gaseous aerosol species
C         GAS(2) = GHNO3
C         GAS(3) = GHCL
CC
C         DO I=1,NIONS              ! Liquid aerosol species
C            AERLIQ(I) = MOLAL(I)
C         ENDDO 
C 
C30    RETURN
C    
C *** END OF SUBROUTINE CALCC2 *****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD3
C *** CASE D3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS OLNY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
C      SUBROUTINE CALCD3!(WI2,WI3,WI4,RHIi,TEMPIi)
      SUBROUTINE CALCD3
      INCLUDE 'isrpia_b.inc'
C
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
C     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
C
      INTEGER :: NDIVOLD
      LOGICAL :: CHNDIVF, BISECT, EARLY, REX, LDIFFX
      DOUBLE PRECISION :: X1, X2, Y1, Y2, X3, Y3, YF, YLO
      DOUBLE PRECISION :: THRSH, DIFF, TSTSIGN, PSI4LO, PSI4HI, P4
      DOUBLE PRECISION :: DIFFX, DIFFXQ
      DOUBLE PRECISION :: FEPS
      LOGICAL :: TST1, TST2, TST
      CHARACTER(40) :: ERRINF
      INTEGER  :: ERRSTKI(25), J
      LOGICAL    :: DEXS,       IEXS,       EOF
      CHARACTER(40)  :: ERRMSGI(25)
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
C
C *** FIND DRY COMPOSITION **********************************************
C
      FEPS   = 1.d-5
      REX = .FALSE.
      CALL CALCD1AL
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3
C
      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CHI2
      PSI3 = ZERO   
      PSI4 = ZERO  
C
      MOLAL(5) = PSI2              ! sc.7.2010  - include dissolved sulfate in initial water calc
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL CALCMR                  ! Initial water
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      TST1   = .TRUE.
      TST2   = .TRUE.
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
60    X1 = PSI4LO
      CALL RSTGAMP
      CALL FUNCD3 (X1, Y1)
      IF (ABS(Y1).LE.(EPS)) THEN
         X3 = X1
         GOTO 50
      ENDIF
      YLO = Y1                 ! Save Y-value at HI position
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      X2 = X1
      Y2 = Y1
      I  = 1
      DO WHILE ((I.LE.NDIV) .AND. TST1)
         X1 = X2
         Y1 = Y2
         X2 = X1+DX 
         CALL FUNCD3 (X2, Y2)
C         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) THEN
         IF (((Y1) .LT. ZERO) .AND. ((Y2) .GT. ZERO)) THEN
             TST1 = .FALSE.! (Y1*Y2.LT.ZERO)
         ENDIF
         I = I + 1
      ENDDO
      IF (.NOT.TST1) GOTO 20
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YHI= Y1                      ! Save Y-value at Hi position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         X3 = X2
         Y3 = Y2
         GOTO 50
C
C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
C Physically I dont know when this might happen, but I have put this
C branch in for completeness. I assume there is no solution; all NO3 goes to the
C gas phase.
C
      ELSE IF ((YLO).LT.ZERO .AND. (YHI).LT.ZERO) THEN
         P4 = TINY ! PSI4LO ! CHI4
         CALL RSTGAMP
         CALL FUNCD3(P4, Y3)
         X3 = P4
         GOTO 50
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
C This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
C and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
C and proceed again with root tracking.
C
      ELSE IF ((YLO).GT.ZERO .AND. (YHI).GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
         IF ((PSI4LO).LT.(-1.D0*(PSI1+PSI2))) THEN
C            WRITE(*,*) 'Error'
            CALL PUSHERR (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
            RETURN
         ELSE
            MOLAL(5) = PSI2              ! so4 included in water calc
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
C            WRITE(*,*) 'Re-executing'
            REX = .TRUE.
            GOTO 60                        ! Redo root tracking
         ENDIF
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    I = 1
      TST2 = .TRUE.
      Y3 = Y2
      DO WHILE ((I.LE.MAXIT) .AND. TST2)
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP
         CALL FUNCD3 (X3,Y3)
         IF (SIGN(1.d0,(Y1))*SIGN(1.d0,(Y3)) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
C         IF (((ABS(X2-X1) .LE. EPS*X1) .OR. (ABS(X2-X1) .LE. TINY))
C     +      .AND. (ABS(Y3).LT.FEPS)) THEN
         IF ((ABS(X2-X1) .LE. EPS*ABS(X1)) .AND. 
     +       (ABS(Y3).LT.FEPS)) THEN
C            WRITE(*,*) 'ABS(X2-X1)',ABS(X2-X1)
C            WRITE(*,*) 'EPS*X1',EPS*X1,'Y3',Y3
            TST2 = .FALSE.
         ENDIF 
C      WRITE(*, '(A,E12.5,A,E12.5)') 'In loop: X3',(X3),'Y3',(Y3)
         I = I + 1
      ENDDO
      IF ((I.GT.MAXIT+1) .AND. TST2) THEN
         CALL PUSHERR (0002, 'CALCD3')    ! WARNING ERROR: NO CONVERGENCE
      ENDIF
C
C *** CONVERGED ; RETURN **********************************************
C
      X3 = 0.5*(X1+X2)
      CALL RSTGAMP
      CALL FUNCD3 (X3, Y3)
C      WRITE(*, '(A,E12.5,A,E12.5)') 'Final X3',X3, 'Y3',Y3
C
50    CONTINUE
C      IF (ABS(Y3).GT. FEPS) THEN
C         WRITE(*, '(A,E12.5,A,E12.5)')'Error: X3',(X3),'Y3',(Y3)
C         WRITE(ERRINF, '(A,E12.5,A)') 'CALCD3 (',(Y3),')'
C         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
C         RETURN
C      ENDIF
C
C      WP = W
      CALL FUNCD3P(X3,YF)!,WP,GAS,AERLIQ)
C
      RETURN
C
C *** END OF SUBROUTINE CALCD3 ******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD3
C *** CASE D3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS OLNY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCD3_B(wpb, gasb, aerliqb)
      INCLUDE 'isrpia_b.inc'
C
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
C     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
C
      INTEGER NDIVOLD
      LOGICAL CHNDIVF, BISECT, EARLY, REX, LDIFFX
      DOUBLE PRECISION X1, X2, Y1, Y2, X3, Y3, YF, YLO
      DOUBLE PRECISION THRSH, DIFF, TSTSIGN, PSI4LO, PSI4HI, P4
      DOUBLE PRECISION DIFFX, DIFFXQ
      DOUBLE PRECISION FEPS
      LOGICAL TST1, TST2, TST
      CHARACTER(40) ERRINF
      INTEGER  ERRSTKI(25), J
      LOGICAL    DEXS,       IEXS,       EOF
      CHARACTER(40)  ERRMSGI(25)
      DOUBLE PRECISION WP(NCOMP), GAS(3), AERLIQ(7)
      DOUBLE PRECISION WPB(NCOMP), GASB(3), AERLIQB(7)
      DOUBLE PRECISION WPDB(NCOMP), GASDB(3), AERLIQDB(7)
C
C *** FIND DRY COMPOSITION **********************************************
C
C      WRITE(*,*) 'W',W, 'gasb',gasb
C      WRITE(*,*) 'aerliqb ',aerliqb
      FEPS   = 1d-5
      REX = .FALSE.
      CALL CALCD1AL
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3
C
      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CHI2
      PSI3 = ZERO   
      PSI4 = ZERO  
C
      MOLAL(5) = PSI2              ! sc.7.2010  - include dissolved sulfate in initial water calc
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL CALCMR                  ! Initial water
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      TST1   = .TRUE.
      TST2   = .TRUE.
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
60    X1 = PSI4LO
      CALL RSTGAMP
      CALL FUNCD3 (X1, Y1)
      IF (ABS(Y1).LE.(EPS)) THEN
         X3 = X1
         GOTO 50
      ENDIF
      YLO = Y1                 ! Save Y-value at HI position
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      X2 = X1
      Y2 = Y1
      I  = 1
      DO WHILE ((I.LE.NDIV) .AND. TST1)
         X1 = X2
         Y1 = Y2
         X2 = X1+DX 
         CALL FUNCD3 (X2, Y2)
C         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) THEN
         IF (((Y1) .LT. ZERO) .AND. ((Y2) .GT. ZERO)) THEN
             TST1 = .FALSE.! (Y1*Y2.LT.ZERO)
         ENDIF
         I = I + 1
      ENDDO
      IF (.NOT.TST1) GOTO 20
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YHI= Y1                      ! Save Y-value at Hi position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         X3 = X2
         Y3 = Y2
         GOTO 50
C
C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
C Physically I dont know when this might happen, but I have put this
C branch in for completeness. I assume there is no solution; all NO3 goes to the
C gas phase.
C
      ELSE IF ((YLO).LT.ZERO .AND. (YHI).LT.ZERO) THEN
         P4 = TINY ! PSI4LO ! CHI4
         CALL RSTGAMP
         CALL FUNCD3(P4, Y3)
         X3 = P4
         GOTO 50
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
C This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
C and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
C and proceed again with root tracking.
C
      ELSE IF ((YLO).GT.ZERO .AND. (YHI).GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
         IF ((PSI4LO).LT.(-1.D0*(PSI1+PSI2)) .OR. PSI4LO.LE.TINY) THEN
C            WRITE(*,*) 'Error'
            CALL PUSHERR (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
            GOTO 50
         ELSE
            MOLAL(5) = PSI2              ! so4 included in water calc
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
C            WRITE(*,*) 'Re-executing: psi4lo',psi4lo,psi1,psi2
            REX = .TRUE.
            GOTO 60                        ! Redo root tracking
         ENDIF
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    I = 1
      TST2 = .TRUE.
      Y3 = Y2
      DO WHILE ((I.LE.MAXIT) .AND. TST2)
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP
         CALL FUNCD3 (X3,Y3)
         IF (SIGN(1.d0,(Y1))*SIGN(1.d0,(Y3)) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
C         IF (((ABS(X2-X1) .LE. EPS*X1) .OR. (ABS(X2-X1) .LE. TINY))
C     +      .AND. (ABS(Y3).LT.FEPS)) THEN
         IF ((ABS(X2-X1) .LE. EPS*ABS(X1)) .AND. 
     +       (ABS(Y3).LT.FEPS)) THEN
C            WRITE(*,*) 'ABS(X2-X1)',ABS(X2-X1)
C            WRITE(*,*) 'EPS*X1',EPS*X1,'Y3',Y3
            TST2 = .FALSE.
         ENDIF 
C      WRITE(*, '(A,E12.5,A,E12.5)') 'In loop: X3',(X3),'Y3',(Y3)
         I = I + 1
      ENDDO
      IF ((I.GT.MAXIT+1) .AND. TST2) THEN
         CALL PUSHERR (0002, 'CALCD3')    ! WARNING ERROR: NO CONVERGENCE
      ENDIF
C
C *** CONVERGED ; RETURN **********************************************
C
C      WRITE(*,*) 'X1',X1,'X2',X2,'X3',X3
      X3 = 0.5*(X1+X2)
      CALL RSTGAMP
      CALL FUNCD3 (X3, Y3)
C      WRITE(*, '(A,E12.5,A,E12.5)') 'Final X3',X3, 'Y3',Y3
C
50    CONTINUE
C      IF (ABS(Y3).GT. FEPS) THEN
C         WRITE(*, '(A,E12.5,A,E12.5)')'Error: X3',(X3),'Y3',(Y3)
C         PAUSE
C         WRITE(ERRINF, '(A,E12.5,A)') 'CALCD3 (',(Y3),')'
C         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
C         npflag = 1
C         ncase = 4
C         DO i = 1,5
C            wpb(i) = 0.D0
C         ENDDO
C         OPEN (199, FILE='adj_sens.csv',STATUS='UNKNOWN',
C     &             POSITION='APPEND')
C         WRITE(199,888) w,rh,temp,wpb,npflag,ncase
C         CLOSE (199, STATUS='KEEP')
C         RETURN
C      ENDIF
C
C      WP = W
C      ncase = 4
C      WRITE(*,*) 'gasb',gasb, ' aerliqb',aerliqb

      wpb = 0.d0
      CALL FUNCD3P_DB(x3, y1, wpdb, gasb, aerliqb)
C      WRITE(*,*) 'wpdb',wpdb
C
      CALL ISERRINF (ERRSTKI, ERRMSGI, NOFER, STKOFL) ! Obtain error stack
      IF (NOFER.EQ.0) THEN ! No errors
          NONPHYS = .FALSE.
C          npflag = 0
C          WRITE(*,*) 'Setting NONPHYS to FALSE'
C          PAUSE
      ELSE
          NONPHYS = .TRUE.
C          npflag = 1
CC          WRITE(*,*) 'Setting NONPHYS to TRUE'
CC          PAUSE
CC          WRITE(*,*) 'After pause'
C      ENDIF
C      OPEN (199, FILE='adj_sens.csv',STATUS='UNKNOWN',
C     &             POSITION='APPEND')
C      WRITE(199,888) w,rh,temp,wpdb,npflag,ncase
C      CLOSE (199, STATUS='KEEP')
C 888  FORMAT (12(1PE11.4,","),I2,",",I2)
CC
      RETURN
C
C *** END OF SUBROUTINE CALCD3 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCD3
C *** CASE D3 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
C
C=======================================================================
C
      SUBROUTINE FUNCD3(P4, FD3)
      INCLUDE 'isrpia_b.inc'
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
C     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      LOGICAL :: TST
      DOUBLE PRECISION :: GMAX, GTHRESH, P4, FD3, BB, DENM
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      TST    = .TRUE.
      PSI4   = P4
      I      = 1
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO WHILE ((I.LE.NSWEEP) .AND. TST)
C
         IF (I.GT.1) CALL CALCACT3
C
         A2   = XK7*(WATER/GAMA(4))**3.0
         A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7   = XKW *RH*WATER*WATER
C
         PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         PSI3 = MIN(MAX(PSI3, ZERO), CHI3)
C
         BB   = PSI4 - PSI3
CCCOLD         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
CCC         AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0
         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF ((DENM).LE.TINY) THEN       ! Avoid overflow when HI->0
C            WRITE(*,*) 'TINY .GT. DENM: ',DENM
            ABB  = ABS(BB)
            DENM = BB + ABB + 2.0*A7/ABB - (2.0*A7*A7)/ABB**3.0 ! Taylor expansion of SQRT
C            WRITE(*,*) 'FUNCD3, TS approx. DENM: ',DENM
         ENDIF
         AHI = 2.0*A7/DENM
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = AHI                             ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         ! NH4I
         MOLAL (5) = PSI2                            ! SO4I
         MOLAL (6) = ZERO                            ! HSO4I
         MOLAL (7) = PSI3 + PSI1                     ! NO3I
         CNH42S4   = CHI2 - PSI2                     ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                            ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                     ! Gas HNO3
         GNH3      = CHI4 - PSI4                     ! Gas NH3
         CALL CALCMR                                 ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            TST = .TRUE.
         ELSE
            TST = .FALSE.
         ENDIF
         I = I+1
      ENDDO
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
CCC      FUNCD3= NH4I/HI/MAXCOMP(GNH3,TINY)/A4 - ONE 
      FD3= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
C      WRITE(*,*) 'fd3 ',fd3
C      WRITE(*,*) 'molal(3)',molal(3),'molal(1)',molal(1)
C      WRITE(*,*) 'GNH3',GNH3
C      WRITE(*,*) 'a4',a4
      RETURN
C
C *** END OF FUNCTION FUNCD3 ********************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCD3P
C *** CASE D3 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
C
C=======================================================================
C
      SUBROUTINE FUNCD3P (P4, Y1)
      INCLUDE 'isrpia_b.inc'
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
C     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      DOUBLE PRECISION P4, Y1, PARM, X
      DOUBLE PRECISION X1, X2, XT, Y1D, Y2, XTD
      DOUBLE PRECISION OMPS, DIAK, DELTA
      CHARACTER(40) ERRINF
      INTEGER   ERRSTKI(25)
      LOGICAL   DEXS,       IEXS,       EOF
      DOUBLE PRECISION ZE, FEPS
      CHARACTER(40) ERRMSGI(25)
      COMMON /CD1A/ OM, PS, ZE
C
C *** SETUP PARAMETERS ************************************************
C
C      W = WP
C
C *** SETUP PARAMETERS ************************************************
C
      FEPS    = 1.D-5
      PARM    = XK10/(R*TEMP)/(R*TEMP)
C
C *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
C
      CNH42S4 = W(2)
      X = MIN(W(3)-2.d0*W(2), W(4))
      IF ((X) .GT. ZERO) THEN
         IF (((W(3) - 2.0*W(2))) .LT. (W(4))) THEN
            PS = ZERO
            OM = W(4) - W(3) + 2.0*W(2)
            IF ((OM) .LT. TINY) THEN
               OM = ZERO
            ENDIF
         ELSE 
            PS = W(3) - W(4) - 2.0*W(2)
            IF ((PS) .LT. TINY) THEN
               PS = ZERO
            ENDIF
            OM = ZERO
         ENDIF
      ELSE
         X  = ZERO
         PS = MAX(W(3) - 2.d0*W(2), ZERO) 
         IF ((PS) .LT. TINY) THEN
            PS = ZERO
         ENDIF
         OM = W(4)
      ENDIF
C
      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
      ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA
C
C *** SPECIATION *******************************************************
C
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3
C
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3
C
      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CHI2
      PSI3 = ZERO
      PSI4 = P4
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
C
C *** NEWTON-RAPHSON DETERMINATION OF ROOT **********************
C
C      WRITE(*,*) 'Before FUNCD3B_DNRD, xt: ',PSI4
      XT  = PSI4
      XTD = 1.D0
C      WRITE(*,*) 'PSI4 ',PSI4
C      CALL FUNCD3B(XT,Y1)
      CALL FUNCD3B_DNRD(XT, XTD, Y1, Y1D)
      X2 = XT - (Y1/(Y1D*1.d0))
      CALL FUNCD3B(X2,Y2)
      IF (abs(Y2).GT. 10.d0*FEPS) THEN
C         WRITE(*,*) '104,D3,',Y2
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCD3 (',(Y2),')'
         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
         RETURN
      ENDIF
C
      IF ((MOLAL(1)).GT.TINY .AND. (MOLAL(5)).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C      PAUSE
C
C         GAS(1) = GNH3                ! Gaseous aerosol species
C         GAS(2) = GHNO3
C         GAS(3) = GHCL
CC
C         DO I=1,NIONS              ! Liquid aerosol species
C            AERLIQ(I) = MOLAL(I)
C         ENDDO 
C 
      RETURN
C
C *** END OF FUNCTION FUNCD3P *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCD3
C *** CASE D3 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
C
C=======================================================================
C
      SUBROUTINE FUNCD3B (P4,FD3B)
      INCLUDE 'isrpia_b.inc'
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
C     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      DOUBLE PRECISION WP(5), MOLALP(7)
      DOUBLE PRECISION P4, BB,  DENM, AHI, AML5, FD3B
      CHARACTER(40)  ERRINF
      INTEGER   ERRSTKI(25), K, J
      LOGICAL   DEXS,       IEXS,       EOF
      CHARACTER(40) ERRMSGI(25)
      LOGICAL TST
C
C *** SETUP PARAMETERS ************************************************
C
      PSI4   = P4
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1,3
C
         A2   = XK7*(WATER/GAMA(4))**3.0
         A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7   = XKW *RH*WATER*WATER
C
         PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         PSI3 = MIN(MAX(PSI3, ZERO), CHI3)
C
         BB   = PSI4 - PSI3
         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF ((DENM).LE.TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = BB + ABB + 2.0*A7/ABB - (2.0*A7*A7)/ABB**3.0 ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.0*A7/DENM
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = AHI                             ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         ! NH4I
         MOLAL (5) = PSI2                            ! SO4I
         MOLAL (6) = ZERO                            ! HSO4I
         MOLAL (7) = PSI3 + PSI1                     ! NO3I
         CNH42S4   = CHI2 - PSI2                     ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                            ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                     ! Gas HNO3
         GNH3      = CHI4 - PSI4                     ! Gas NH3
C
C         CALL CALCMR                                ! Water content
C
         MOLALR(4) = MOLAL(5) + MOLAL(6)             ! (NH4)2SO4
         AML5      = MOLAL(3)-2.D0*MOLALR(4)         ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,MOLAL(7)), ZERO)   ! NH4NO3 = MIN("free", NO3)
C
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C
         CALL CALCACT3F
      ENDDO
C
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
CCC      FUNCD3= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE 
      FD3B = MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
      RETURN
C
C *** END OF FUNCTION FUNCD3P *******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD1A
C *** CASE D1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
C
C     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
C     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
C     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
C     THE SOLID PHASE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCD1A
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: PARM, PS, OM, OMPS, DIAK, ZE, X
C
C *** SETUP PARAMETERS ************************************************
C
      PARM    = XK10/(R*TEMP)/(R*TEMP)
C
C *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
C
      CNH42S4 = W(2)                                    
      X       = MAX(MIN(W(3)-2.d0*CNH42S4, W(4)), ZERO)  ! MAX NH4NO3
      PS      = MAX(W(3) - X - 2.d0*CNH42S4, ZERO)
      OM      = MAX(W(4) - X, ZERO)
C
      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
      ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA
C
C *** SPECIATION *******************************************************
C
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3
C
C      PAUSE
      RETURN
C
C *** END OF SUBROUTINE CALCD1A *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD1A
C *** CASE D1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
C
C     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
C     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
C     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
C     THE SOLID PHASE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCD1AL
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: PARM, PS, OM, OMPS, DIAK, ZE, X
C
C *** SETUP PARAMETERS ************************************************
C
      PARM    = XK10/(R*TEMP)/(R*TEMP)
C
C *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
C
      CNH42S4 = W(2)
      X = MIN(W(3)-2.0*W(2), W(4))
      IF ((X) .GT. ZERO) THEN
         IF ((W(3) - 2.0*W(2)) .LT. (W(4))) THEN
            PS = ZERO
            OM = W(4) - W(3) + 2.0*W(2)
            IF ((OM) .LT. TINY) THEN
               OM = ZERO
            ENDIF
         ELSE 
            PS = W(3) - W(4) - 2.0*W(2)
            IF ((PS) .LT. TINY) THEN
               PS = ZERO
            ENDIF
            OM = ZERO
         ENDIF
      ELSE
         X  = ZERO
         PS = MAX(W(3) - 2.0*W(2), ZERO) 
         IF ((PS) .LT. TINY) THEN
            PS = ZERO
         ENDIF
         OM = W(4)
      ENDIF
C
      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
      ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA
C
C *** SPECIATION *******************************************************
C
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3
C
      RETURN
C
C *** END OF SUBROUTINE CALCD1A *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCG5
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCG5
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST1, TST2
      INTEGER I
      CHARACTER(40)  ERRINF
      INTEGER ERRSTKI(25)
      CHARACTER(40) ERRMSGI(25)
      DOUBLE PRECISION AERLIQ(12), GAS(3), WP(5)
      DOUBLE PRECISION LAMDA
C
C *** SETUP PARAMETERS ************************************************
C
C      WRITE(*,*) 'Beginning CALCG5'
      CALAOU = .TRUE.   
      CHI1   = 0.5d0*W(1)
      CHI2   = MAX (W(2)-0.5d0*W(1), ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI1   = CHI1
      PSI2   = CHI2
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
      WATER  = CHI2/M0(4) + CHI1/M0(2)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      CALL FUNCG5A (X1, Y1)
      IF (CHI6.LE.TINY) THEN 
         X3 = X1
         Y3 = Y1
         GOTO 50  
      ENDIF
ccc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
ccc      IF (WATER .LE. TINY) RETURN                    ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      X2 = X1 
      Y2 = Y1
C      DO I=1,NDIV
C      WRITE(*,*) 'NDIV',NDIV
      I = 1
      TST1 = .TRUE.
      DO WHILE ((I.LE.NDIV) .AND. TST1)
         X1 = X2
         Y1 = Y2
         X2 = X1+DX 
         CALL FUNCG5A (X2, Y2)
C         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         IF ((Y1 .LT. ZERO) .AND. (Y2 .GT. ZERO)) THEN
             TST1 = .FALSE.! (Y1*Y2.LT.ZERO)
         ENDIF
         I = I+1
      ENDDO
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
C      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG5A (PSI6LO)
C      GOTO 50
      IF ((ABS(Y2).GT.EPS).AND.TST1.AND.(I.GT.NDIV+1)) THEN
         CALL RSTGAMP
         CALL FUNCG5A (PSI6LO, Y3)
         X3 = PSI6LO
         CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
C         WRITE(*,*) 'No subdivision with solution found'      
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      I    = 1
      TST2 = .TRUE.
      FEPS = 1.D-5
      DO WHILE ((I.LE.MAXIT).AND.TST2)
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP
         CALL FUNCG5A (X3, Y3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1 .AND. (ABS(Y3).LT.FEPS)) THEN
            TST2 = .FALSE.   !GOTO 40
         ENDIF
         I = I+1
      ENDDO
      IF ((I.GT.(MAXIT+1)) .AND. TST2) THEN
         CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
      ENDIF
C
C *** CONVERGED ; RETURN **********************************************
C
      X3 = 0.5*(X1+X2)
      CALL RSTGAMP
      CALL FUNCG5A (X3, Y3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
C *** Execute differentiable Newton's function once ***********************
C      IF (ABS(Y3).GT. FEPS) THEN
C         WRITE(*, '(A,E12.5,A,E12.5)') 'Error D3: X3',X3, 'Y3',Y3
C         WRITE(ERRINF, '(A,E12.5,A)') 'CALCG5 (',Y3,')'
C         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
C      ENDIF
C
C      WP = W
      CALL FUNCG5AP(X3)!,WP,GAS,AERLIQ)
C
      CALL ISERRINF (ERRSTKI, ERRMSGI, NOFER, STKOFL) ! Obtain error stack
      IF (NOFER.EQ.0) THEN ! No errors
          NONPHYS = .FALSE.
      ELSE
          NONPHYS = .TRUE.
C          WRITE(*,*) 'Setting NONPHYS to TRUE'
C          PAUSE
C          WRITE(*,*) 'After pause'
      ENDIF
C
C      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
C         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
C         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
C         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
C         MOLAL(6) = DELTA                               ! HSO4 EFFECT
C      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCG5 *******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCG5
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCG5_B(wpb, gasb, aerliqb)
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST1, TST2
      INTEGER I
      CHARACTER(40)  ERRINF
      INTEGER ERRSTKI(25)
      CHARACTER(40) ERRMSGI(25)
      DOUBLE PRECISION AERLIQ(12), GAS(3), WP(5)
      DOUBLE PRECISION LAMDA, wpb(5)
      DOUBLE PRECISION gasb(ngasaq), aerliqb(nions+ngasaq+2)
C
C *** SETUP PARAMETERS ************************************************
C
C      WRITE(*,*) 'Beginning CALCG5'
      FEPS   = 1.d-5
      CALAOU = .TRUE.   
      CHI1   = 0.5d0*W(1)
      CHI2   = MAX (W(2)-0.5d0*W(1), ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI1   = CHI1
      PSI2   = CHI2
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
      WATER  = CHI2/M0(4) + CHI1/M0(2)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      CALL FUNCG5A (X1, Y1)
      IF (CHI6.LE.TINY) THEN 
         X3 = X1
         Y3 = Y1
         GOTO 50  
      ENDIF
ccc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
ccc      IF (WATER .LE. TINY) RETURN                    ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      X2 = X1 
      Y2 = Y1
C      DO I=1,NDIV
C      WRITE(*,*) 'NDIV',NDIV
      I = 1
      TST1 = .TRUE.
      DO WHILE ((I.LE.NDIV) .AND. TST1)
         X1 = X2
         Y1 = Y2
         X2 = X1+DX 
         CALL FUNCG5A (X2, Y2)
C         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         IF ((Y1 .LT. ZERO) .AND. (Y2 .GT. ZERO)) THEN
             TST1 = .FALSE.! (Y1*Y2.LT.ZERO)
         ENDIF
         I = I+1
      ENDDO
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
C      IF (ABS(Y2) .GT. EPS) Y2 = FUNCG5A (PSI6LO)
C      GOTO 50
      IF ((ABS(Y2).GT.EPS).AND.TST1.AND.(I.GT.NDIV+1)) THEN
         CALL RSTGAMP
         CALL FUNCG5A (PSI6LO, Y3)
         X3 = PSI6LO
         CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
C         WRITE(*,*) 'No subdivision with solution found'      
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      I    = 1
      TST2 = .TRUE.
      FEPS = 1.D-5
      DO WHILE ((I.LE.MAXIT).AND.TST2)
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP
         CALL FUNCG5A (X3, Y3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1 .AND. (ABS(Y3).LT.FEPS)) THEN
            TST2 = .FALSE.   !GOTO 40
         ENDIF
         I = I+1
      ENDDO
      IF ((I.GT.(MAXIT+1)) .AND. TST2) THEN
         CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
      ENDIF
C
C *** CONVERGED ; RETURN **********************************************
C
      X3 = 0.5*(X1+X2)
      CALL RSTGAMP
      CALL FUNCG5A (X3, Y3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE

C      slc.debug
C      WRITE(*,*) '---- Before FUNCG5AP_GB ----'
C      WRITE(*,*) 'aerliqb: ',aerliqb
C      WRITE(*,*) 'gasb: ',gasb
C *** Execute differentiable Newton's function once ***********************
C      IF (ABS(Y3).GT. FEPS) THEN
C         WRITE(*,*) '104 G5: ',Y3
C         WRITE(*, '(A,E12.5,A,E12.5)') 'Error G5: X3',X3, 'Y3',Y3
C         WRITE(ERRINF, '(A,E12.5,A)') 'CALCG5 (',Y3,')'
C         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
C      ENDIF
C
C      WP = W
C      CALL FUNCG5AP(X3,WP,GAS,AERLIQ)

      wpb = 0.d0
      CALL FUNCG5AP_GB(x3, wpb, gasb, aerliqb)
C
C      ncase = 7
      CALL ISERRINF (ERRSTK, ERRMSG, NOFER, STKOFL) ! Obtain error stack
C      WRITE(*,*) 'Writing error code'
      IF (NOFER.EQ.0) THEN ! No errors
          NONPHYS = .FALSE.
          npflag = 0
C          WRITE(*,*) 'Setting NONPHYS to FALSE'
      ELSE
          NONPHYS = .TRUE.
          npflag = 1
C          WRITE(*,*) 'Setting NONPHYS to TRUE'
C          PAUSE
C          WRITE(*,*) 'After pause'
      ENDIF
C      WRITE(*,*) 'wpb',wpb
C      OPEN (199, FILE='adj_sens.csv',STATUS='UNKNOWN',
C     &             POSITION='APPEND')
C      WRITE(199,888) w,rh,temp,wpb,npflag,ncase
C      CLOSE (199, STATUS='KEEP')
C 888  FORMAT (12(1PE11.4,","),I2,",",I2)
C
C      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
C         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
C         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
C         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
C         MOLAL(6) = DELTA                               ! HSO4 EFFECT
C      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCG5 *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCG5A
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE FUNCG5A (X, FG5A)
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST
      INTEGER SO4FLG
      DOUBLE PRECISION LAMDA, FG5A
C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      I      = 1
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      TST = .TRUE.
      DO WHILE ((I.LE.NSWEEP).AND. TST)
C      WRITE(*,*) 'NSWEEP ',NSWEEP,'WATER',WATER
C
C      WRITE(*,*) 'GAMA ', GAMA
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      AKK = A4*A6
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF
C
CCC      IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = W(1)                                ! NAI
C      MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
      MOLAL (4) = PSI6                                ! CLI
      IF (W(2)-0.5d0*W(1) .GT. ZERO) THEN
         MOLAL(3) = 2.d0*W(2) - W(1) + PSI4
         MOLAL(5) = W(2)                             ! SO4I
      ELSE
         MOLAL(3) = PSI4
         MOLAL(5) = 0.5d0*W(1)                       ! SO4I
      ENDIF
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5                                ! NO3I
C
C      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      SMIN     = PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
      GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
      GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
C
      CNH42S4   = ZERO                                ! Solid (NH4)2SO4
      CNH4NO3   = ZERO                                ! Solid NH4NO3
      CNH4CL    = ZERO                                ! Solid NH4Cl
C
C      CALL CALCMR                                     ! Water content
C
C      WRITE(*,*) 'MOLAL ',MOLAL
         MOLALR(2) = 0.5*W(1)                          ! NA2SO4
         IF ((W(2)-0.5d0*W(1)) .GT. ZERO) THEN
            TOTS4     = W(2)                         ! Total SO4
            MOLALR(4) = W(2)-0.5d0*W(1)              ! (NH4)2SO4
            FRNH4     = MAX(PSI4, ZERO)
         ELSE
            TOTS4     = 0.5d0*W(1)                   ! Total SO4
            MOLALR(4) = ZERO          ! (NH4)2SO4
            FRNH4     = MAX(2.d0*W(2)-W(1) + PSI4, ZERO)
         ENDIF
         IF ((PSI5) .LT. (FRNH4)) THEN
            MOLALR(5) = PSI5
            FRNH4     = MAX(FRNH4 - PSI5, ZERO)
         ELSE 
            MOLALR(5) = FRNH4
            FRNH4     = ZERO
         ENDIF
         MOLALR(6) = MIN(PSI6, FRNH4)                  ! NH4CL
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C         WRITE(*,*) 'After CALCMR: WATER ',WATER
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         TST = .TRUE.
      ELSE 
         TST = .FALSE.
      ENDIF
      CALL CALCACT3
      I = I + 1
C
      ENDDO
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    FG5A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
CCC         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCG5A *******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCG5AB
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE FUNCG5AB (X, FG5AB)
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST
      INTEGER SO4FLG
      DOUBLE PRECISION LAMDA, FG5AB
C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
C      WRITE(*,*) 'NSWEEP ',NSWEEP,'WATER',WATER
      DO I = 1,2
C
C      IF (I.GT.1) CALL CALCACT3
C      WRITE(*,*) 'GAMA ', GAMA
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      AKK = A4*A6
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF
C
CCC      IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = W(1)                                ! NAI
C      MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
      MOLAL (4) = PSI6                                ! CLI
      IF (W(2)-0.5d0*W(1) .GT. ZERO) THEN
         MOLAL(3) = 2.d0*W(2) - W(1) + PSI4
         MOLAL(5) = W(2)                             ! SO4I
      ELSE
         MOLAL(3) = PSI4
         MOLAL(5) = 0.5d0*W(1)                       ! SO4I
      ENDIF
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5                                ! NO3I
C
C      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      SMIN     = PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
      GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
      GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
C
      CNH42S4   = ZERO                                ! Solid (NH4)2SO4
      CNH4NO3   = ZERO                                ! Solid NH4NO3
      CNH4CL    = ZERO                                ! Solid NH4Cl
C
C      CALL CALCMR                                     ! Water content
C
C      WRITE(*,*) 'MOLAL ',MOLAL
         MOLALR(2) = 0.5*W(1)                          ! NA2SO4
         IF (W(2)-0.5d0*W(1) .GT. ZERO) THEN
            TOTS4     = W(2)                         ! Total SO4
            MOLALR(4) = W(2)-0.5d0*W(1)              ! (NH4)2SO4
            FRNH4     = MAX(PSI4, ZERO)
         ELSE
            TOTS4     = 0.5d0*W(1)                   ! Total SO4
            MOLALR(4) = ZERO          ! (NH4)2SO4
            FRNH4     = MAX(2.d0*W(2)-W(1) + PSI4, ZERO)
         ENDIF
         IF (PSI5 .LT. FRNH4) THEN
            MOLALR(5) = PSI5
            FRNH4     = MAX(FRNH4 - PSI5, ZERO)
         ELSE 
            MOLALR(5) = FRNH4
            FRNH4     = ZERO
         ENDIF
         MOLALR(6) = MIN(PSI6, FRNH4)                  ! NH4CL
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C         WRITE(*,*) 'After CALCMR: WATER ',WATER
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         CALL CALCACT3F
      ENDDO
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    FG5AB = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
CCC         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCG5A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCG5A
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE FUNCG5AP (X1) !, WP, GAS, AERLIQ)
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST
      CHARACTER(40)  ERRINF
      DOUBLE PRECISION LAMDA, FEPS
      DOUBLE PRECISION WP(5),  X1
      DOUBLE PRECISION AERLIQ(12), GAS(3)
C
C *** SETUP PARAMETERS ************************************************
C
C      DO I = 1,5
C         W(I) = WP(I)
C      ENDDO
      FEPS   = 1.d-5
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-0.5d0*W(1), ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI1   = CHI1
      PSI2   = CHI2
      I      = 1
      PSI6   = X1
      FRST   = .TRUE.
      CALAIN = .TRUE. 
      TST    = .TRUE.
C
C *** NEWTON-RAPHSON DETERMINATION OF ROOT **********************
C
      XT  = X1
      XTD = 1.D0
CCCC$AD NOCHECKPOINT
      CALL FUNCG5AB_GNRD(XT, XTD, Y1, Y1D)
      X2 = XT - (Y1/(Y1D*1.d0))
      CALL FUNCG5AB(X2,Y2)
      IF (abs(Y2).GT. 10.d0*FEPS) THEN
C         WRITE(*,*) '104,G5,',Y2
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCG5 (',(Y1),')'
         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
         RETURN
      ENDIF
C      CALL FUNCG5AB(XT,Y2)
C
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
C      WRITE(*,*) 'MOLAL', MOLAL
C
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQ(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQ(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQ(NIONS+NGASAQ+2) = COH
CC
C      GAS(1) = GNH3                ! Gaseous aerosol species
C      GAS(2) = GHNO3
C      GAS(3) = GHCL
C 
C
C *** END OF FUNCTION FUNCG5A *******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCH6
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCH6
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST1, TST2, TST
      DOUBLE PRECISION FEPS, WP(5), AERLIQ(12), GAS(3)
      CHARACTER(40)  ERRINF
      INTEGER ERRSTKI(25)
      CHARACTER(40) ERRMSGI(25)
C
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
C     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
C     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
C     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
C     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      TST1   = .TRUE.
      TST2   = .TRUE.
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*W(2), ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
C      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
C      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
C      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
      IF (FRNA .LT. W(4)) THEN
         CHI5 = MAX(W(4)-FRNA, ZERO)
         CHI7 = MIN(ZERO,W(5))
         CHI6 = MAX(W(5),ZERO)
      ELSE 
         CHI5 = ZERO
         IF (MAX(FRNA-W(4),ZERO) .LT. W(5)) THEN
            CHI7 = MAX(FRNA-W(4),ZERO)
            CHI6 = MAX(W(5)-CHI7,ZERO)
         ELSE
            CHI7 = W(5)
            CHI6 = ZERO
         ENDIF
      ENDIF
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      CALL FUNCH6A (X1, Y1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) THEN
         X3 = X1
         Y3 = Y1 
         GOTO 50  
      ENDIF
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      I = 1
      X2 = X1
      Y2 = Y1
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO WHILE ((I.LE.NDIV) .AND. TST1)
         X1 = X2
         Y2 = Y2         
         X2 = X1+DX 
         CALL FUNCH6A (X2, Y2)
         IF ((Y1 .LT. ZERO) .AND. (Y2 .GT. ZERO)) THEN
             TST1 = .FALSE.! (Y1*Y2.LT.ZERO)
         ENDIF
         I = I+1
      ENDDO 
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF ((ABS(Y2).GT.EPS).AND.TST1.AND.(I.GT.NDIV+1)) THEN
         CALL RSTGAMP
         CALL FUNCH6A (PSI6LO, Y3)
         X3 = PSI6LO
         CALL PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
C         WRITE(*,*) 'No subdivision with solution found'      
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      I = 1
      TST2 = .TRUE.
      FEPS = 1.D-5
      DO WHILE ((I .LE. MAXIT) .AND. TST2)
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP
         CALL FUNCH6A (X3, Y3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ((ABS(X2-X1) .LE. EPS*X1) .AND. (ABS(Y3).LT.FEPS)) THEN
            TST2 = .FALSE.
         ENDIF
         I = I+1
      ENDDO
      IF ((I.GT.(MAXIT+1)) .AND. TST2) THEN
         CALL PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
      ENDIF
C
C *** CONVERGED ; RETURN **********************************************
C
      X3 = 0.5*(X1+X2)
      CALL RSTGAMP
      CALL FUNCH6A (X3, Y3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
C
C      WP = W
      CALL FUNCH6AP(X3)!, WP, GAS, AERLIQ)
C
CC      CALL ISERRINF (ERRSTKI, ERRMSGI, NOFER, STKOFL) ! Obtain error stack
C      IF (NOFER.EQ.0) THEN ! No errors
C          NONPHYS = .FALSE.
C      ELSE
C          NONPHYS = .TRUE.
CC          WRITE(*,*) 'Setting NONPHYS to TRUE'
CC          PAUSE
CC          WRITE(*,*) 'After pause'
C      ENDIF
C      
C      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
C         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
C         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
C         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
C         MOLAL(6) = DELTA                                ! HSO4 EFFECT
C      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCH6 ******************************************
C
      END

C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCH6
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCH6_B(wpb, gasb, aerliqb)
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST1, TST2, TST
      DOUBLE PRECISION FEPS, WP(5), AERLIQ(12), GAS(3), wpb(5)
      CHARACTER(40)  ERRINF
      INTEGER ERRSTKI(25), npflag
      CHARACTER(40) ERRMSGI(25)
      DOUBLE PRECISION gasb(3), aerliqb(12)
C
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
C     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
C     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
C     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
C     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      FEPS   = 1.D-5
      CALAOU = .TRUE.   
      TST1   = .TRUE.
      TST2   = .TRUE.
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*W(2), ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
C      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
C      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
C      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
      IF (FRNA .LT. W(4)) THEN
         CHI5 = MAX(W(4)-FRNA, ZERO)
         CHI7 = MIN(ZERO,W(5))
         CHI6 = MAX(W(5),ZERO)
      ELSE 
         CHI5 = ZERO
         IF (MAX(FRNA-W(4),ZERO) .LT. W(5)) THEN
            CHI7 = MAX(FRNA-W(4),ZERO)
            CHI6 = MAX(W(5)-CHI7,ZERO)
         ELSE
            CHI7 = W(5)
            CHI6 = ZERO
         ENDIF
      ENDIF
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      CALL FUNCH6A (X1, Y1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) THEN
         X3 = X1
         Y3 = Y1 
         GOTO 50  
      ENDIF
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      I = 1
      X2 = X1
      Y2 = Y1
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO WHILE ((I.LE.NDIV) .AND. TST1)
         X1 = X2
         Y2 = Y2         
         X2 = X1+DX 
         CALL FUNCH6A (X2, Y2)
         IF ((Y1 .LT. ZERO) .AND. (Y2 .GT. ZERO)) THEN
             TST1 = .FALSE.! (Y1*Y2.LT.ZERO)
         ENDIF
         I = I+1
      ENDDO 
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF ((ABS(Y2).GT.EPS).AND.TST1.AND.(I.GT.NDIV+1)) THEN
         CALL RSTGAMP
         CALL FUNCH6A (PSI6LO, Y3)
         X3 = PSI6LO
         CALL PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
C         WRITE(*,*) 'No subdivision with solution found'      
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      I = 1
      TST2 = .TRUE.
      DO WHILE ((I .LE. MAXIT) .AND. TST2)
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP
         CALL FUNCH6A (X3, Y3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF ((ABS(X2-X1) .LE. EPS*X1) .AND. (ABS(Y3).LT.FEPS)) THEN
            TST2 = .FALSE.
         ENDIF
         I = I+1
      ENDDO
      IF ((I.GT.(MAXIT+1)) .AND. TST2) THEN
         CALL PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
      ENDIF
C
C *** CONVERGED ; RETURN **********************************************
C
      X3 = 0.5*(X1+X2)
      CALL RSTGAMP
      CALL FUNCH6A (X3, Y3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
C      slc.debug
C      WRITE(*,*) '---- Before FUNCH6AP_HB ----'
C      WRITE(*,*) 'aerliqb: ',aerliqb
C      WRITE(*,*) 'gasb: ',gasb
C      IF (ABS(Y3).GT. FEPS) THEN
C         WRITE(*, '(A,E12.5,A,E12.5)') 'Error H6: X3',X3, 'Y3',Y3
C         WRITE(ERRINF, '(A,E12.5,A)') 'CALCH6 (',Y3,')'
C         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
C      ENDIF
C
C      WP = W
C      CALL FUNCH6AP(X3, WP, GAS, AERLIQ)

      wpb = 0.d0
      CALL FUNCH6AP_HB(x3, wpb, gasb, aerliqb)
C
C      ncase = 8
      CALL ISERRINF (ERRSTKI, ERRMSGI, NOFER, STKOFL) ! Obtain error stack
      IF (NOFER.EQ.0) THEN ! No errors
          NONPHYS = .FALSE.
C          npflag = 0
      ELSE
          NONPHYS = .TRUE.
C          npflag = 1
      ENDIF
C
C      OPEN (199, FILE='adj_sens.csv',STATUS='UNKNOWN',
C     &             POSITION='APPEND')
C      WRITE(199,888) w,rh,temp,wpb,npflag,ncase
C      CLOSE (199, STATUS='KEEP')
C 888  FORMAT (12(1PE11.4,","),I2,",",I2)
C      
C      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
C         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
C         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
C         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
C         MOLAL(6) = DELTA                                ! HSO4 EFFECT
C      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCH6 ******************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCH6A
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE FUNCH6A (X, FH6A)
      INCLUDE 'isrpia_b.inc'
C
      LOGICAL TST
C
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
C     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
C     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
C     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
C     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      I      = 1
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE.
      TST    = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO WHILE ((I .LE. NSWEEP) .AND. TST)
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
      MOLAL (3) = PSI4                                  ! NH4I
      MOLAL (4) = PSI6 + PSI7                           ! CLI
      MOLAL (5) = PSI2 + PSI1                           ! SO4I
      MOLAL (6) = ZERO                                  ! HSO4I
      MOLAL (7) = PSI5 + PSI8                           ! NO3I
C
C      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      SMIN      = 2.d0*PSI2 + PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C
      CALL CALCMR                                    ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         TST = .TRUE.
      ELSE 
         TST = .FALSE.
      ENDIF
      CALL CALCACT3
      I = I + 1
      ENDDO
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
      FH6A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCH6A *******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCH6AB
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE FUNCH6AP (X1)!, WP, GAS, AERLIQ)
      INCLUDE 'isrpia_b.inc'
C
      CHARACTER(40)  ERRINF
      DOUBLE PRECISION AERLIQ(12), GAS(3), WP(5), FEPS
C
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
C     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
C     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
C     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
C     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      FEPS   = 1.d-5
C      W      = WP
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*W(2), ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
C      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
C      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
C      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
      IF (FRNA .LT. W(4)) THEN
         CHI5 = MAX(W(4)-FRNA, ZERO)
         CHI7 = MIN(ZERO,W(5))
         CHI6 = MAX(W(5),ZERO)
      ELSE 
         CHI5 = ZERO
         IF (MAX(FRNA-W(4),ZERO) .LT. W(5)) THEN
            CHI7 = MAX(FRNA-W(4),ZERO)
            CHI6 = MAX(W(5)-CHI7,ZERO)
         ELSE
            CHI7 = W(5)
            CHI6 = ZERO
         ENDIF
      ENDIF
C
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      PSI6   = X1
C
C *** NEWTON-RAPHSON DETERMINATION OF ROOT **********************
C
      XT  = X1
      XTD = 1.D0
      CALL FUNCH6AB_HNRD(XT, XTD, Y1, Y1D)
      X2 = XT - (Y1/(Y1D*1.d0))
      CALL FUNCH6AB(X2,Y2)
      IF (abs(Y2).GT. 10.d0*FEPS) THEN
C         WRITE(*,*) '104,H6,',Y2
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCH6 (',(Y2),')'
         CALL PUSHERR (0104, ERRINF)    ! WARNING ERROR: NO CONVERGENCE
         RETURN
      ENDIF
C      CALL FUNCH6AB(XT,Y2)
C
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQ(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQ(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQ(NIONS+NGASAQ+2) = COH
CC
C      GAS(1) = GNH3                ! Gaseous aerosol species
C      GAS(2) = GHNO3
C      GAS(3) = GHCL
C
C *** END OF FUNCTION FUNCH6A *******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCH6AB
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE FUNCH6AB (X, FH6AB)
      INCLUDE 'isrpia_b.inc'
C
      INTEGER J
C
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
C     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
C     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
C     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
C     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1,2
C
         A1  = XK5 *(WATER/GAMA(2))**3.0
         A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
         A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
         A7  = XK8 *(WATER/GAMA(1))**2.0
         A8  = XK9 *(WATER/GAMA(3))**2.0
         A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
         PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
         PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
         PSI5 = MAX(PSI5, TINY)
C
         IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
            BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
            CC   = CHI4*(PSI5+PSI6)
            DD   = BB*BB-4.d0*CC
            PSI4 =0.5d0*(-BB - SQRT(DD))
            PSI4 = MIN(PSI4,CHI4)
         ELSE
            PSI4 = TINY
         ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
         MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
         MOLAL (3) = PSI4                                  ! NH4I
         MOLAL (4) = PSI6 + PSI7                           ! CLI
         MOLAL (5) = PSI2 + PSI1                           ! SO4I
         MOLAL (6) = ZERO                                  ! HSO4I
         MOLAL (7) = PSI5 + PSI8                           ! NO3I
C
C      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
         SMIN      = 2.d0*PSI2 + PSI5 + PSI6 - PSI4
         CALL CALCPH (SMIN, HI, OHI)
         MOLAL (1) = HI
C 
         GNH3      = MAX(CHI4 - PSI4, TINY)
         GHNO3     = MAX(CHI5 - PSI5, TINY)
         GHCL      = MAX(CHI6 - PSI6, TINY)
C 
         CNH42S4   = ZERO
         CNH4NO3   = ZERO
         CNACL     = MAX(CHI7 - PSI7, ZERO)
         CNANO3    = MAX(CHI8 - PSI8, ZERO)
         CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C
C      CALL CALCMR                                    ! Water content
C
C *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM RICH CASE
C *** RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
C
         MOLALR(1) = PSI7                                  ! NACL 
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI8                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
C         FRNO3     = MAX(MOLAL(7) - MOLALR(3), ZERO)       ! "FREE" NO3
         FRNO3     = MAX(PSI5, ZERO)
C         FRCL      = MAX(MOLAL(4) - MOLALR(1), ZERO)       ! "FREE" CL
         FRCL      = MAX(PSI6, ZERO)
C         MOLALR(5) = MIN(MOLAL(3),FRNO3)                   ! NH4NO3
C         FRNH4     = MAX(MOLAL(3) - MOLALR(5), ZERO)       ! "FREE" NH3
C         MOLALR(6) = MIN(FRCL, FRNH4)                      ! NH4CL
         IF (PSI4 .LT. FRNO3) THEN
            MOLALR(5) = PSI4
            FRNH4     = ZERO 
            MOLALR(6) = MIN(FRCL, ZERO)
         ELSE 
            MOLALR(5) = FRNO3
            FRNH4     = MAX(PSI4-FRNO3,ZERO)
            MOLALR(6) = MIN(FRCL, FRNH4)
         ENDIF
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO J=1,NPAIR
            WATER = WATER + MOLALR(J)/M0(J)
         ENDDO
         WATER = MAX(WATER, TINY)
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         CALL CALCACT3F
      ENDDO
C
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
      FH6AB = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCH6A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCI6
C *** CASE I6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCI6
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION AERLIQ(12), GAS(3), WP(5)
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
C     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
C
C *** FIND DRY COMPOSITION **********************************************
C
C      CALL CALCI1A
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CNA2SO4 = 0.5D0*W(1)
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      FRSO4   = MAX(W(2)-CNA2SO4, ZERO)
C
      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
C
      IF (FRSO4.LE.TINY) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF (FRNH4.LE.TINY) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
      ENDIF
C
C *** CALCULATE GAS SPECIES *********************************************
C
      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      J = 1
      DO WHILE ((J.LE.NSWEEP).AND.(CALAIN))
C
         A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
         BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
         CC   =-A6*(PSI2 + PSI3 + PSI1)
         DD   = BB*BB - 4.D0*CC
         PSI6 = 0.5D0*(-BB + SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
         MOLAL (1) = PSI6                                    ! HI
         MOLAL (2) = 2.D0*PSI4 + PSI3                        ! NAI
         MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1            ! NH4I
         MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6               ! SO4I
         MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6               ! HSO4I
         CLC       = ZERO
         CNAHSO4   = ZERO
         CNA2SO4   = CHI4 - PSI4
         CNH42S4   = ZERO
         CNH4HS4   = ZERO
C      CALL CALCMR                                         ! Water content
         MOLALR(04) = PSI5                                 ! (NH4)2SO4
         MOLALR(02) = PSI4                                 ! NA2SO4
         MOLALR(09) = PSI1                                 ! NH4HSO4
         MOLALR(12) = PSI3                                 ! NAHSO4
         MOLALR(13) = PSI2                                 ! LC
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO I=1,NPAIR
            WATER = WATER + MOLALR(I)/M0(I)
         ENDDO
         WATER = MAX(WATER, TINY)
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         CALL CALCACT3
         J = J+1
      ENDDO
C      WRITE(*,*) 'I6: J',J,' CALAIN', CALAIN,'NSWEEP',NSWEEP
C      PAUSE
      IF (CALAIN .AND. (J .GT. (NSWEEP+1))) THEN
C         WRITE(*,*) 'Error writing'
         CALL PUSHERR (0001, 'CALCI6')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C 
20    RETURN
C
C *** END OF SUBROUTINE CALCI6 *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCJ3
C *** CASE J3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCJ3
      INCLUDE 'isrpia_b.inc'
C
      DOUBLE PRECISION LAMDA, KAPA
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      CHI1   = W(1)                           ! NA TOTAL as NaHSO4
      CHI2   = W(3)                           ! NH4 TOTAL as NH4HSO4
      PSI1   = CHI1
      PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      J = 1
      DO WHILE ((J.LE.NSWEEP).AND.(CALAIN))
C
         A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
C
C  CALCULATE DISSOCIATION QUANTITIES
C
         BB   = A3+LAMDA                        ! KAPA
         CC   =-A3*(LAMDA + PSI1 + PSI2)
         DD   = BB*BB-4.D0*CC
         KAPA = 0.5D0*(-BB+SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
         MOLAL (1) = LAMDA + KAPA                 ! HI
         MOLAL (2) = PSI1                         ! NAI
         MOLAL (3) = PSI2                         ! NH4I
         MOLAL (4) = ZERO                         ! CLI
         MOLAL (5) = KAPA                         ! SO4I
         MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA   ! HSO4I
         MOLAL (7) = ZERO                         ! NO3I
C
         CNAHSO4   = ZERO
         CNH4HS4   = ZERO
C
C      CALL CALCMR                              ! Water content
C
         MOLALR(09) = MOLAL(3)                             ! NH4HSO4
         MOLALR(12) = MOLAL(2)                             ! NAHSO4
         MOLALR(07) = MOLAL(5)+MOLAL(6)-MOLAL(3)-MOLAL(2)  ! H2SO4
         MOLALR(07) = MAX(MOLALR(07),ZERO)
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
         WATER = ZERO
         DO I=1,NPAIR
            WATER = WATER + MOLALR(I)/M0(I)
         ENDDO
         WATER = MAX(WATER, TINY)
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         CALL CALCACT3
         J = J+1
      ENDDO
C      WRITE(*,*) 'J3: J',J,' CALAIN', CALAIN,'NSWEEP',NSWEEP
C      PAUSE
      IF (CALAIN .AND. (J .GT. (NSWEEP+1))) THEN
C         WRITE(*,*) 'Error writing'
         CALL PUSHERR (0001, 'CALCJ3')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C 
50    RETURN
C
C *** END OF SUBROUTINE CALCJ3 ******************************************
C
      END
C

