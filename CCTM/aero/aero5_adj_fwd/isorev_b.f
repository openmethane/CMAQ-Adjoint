C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP1R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF 
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY 
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP1R (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION WI(NCOMP)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
C
C *** INITIALIZE COMMON BLOCK VARIABLES *********************************
C
      CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
C      IF (RH.GE.DRNH42S4) THEN         ! WET AEROSOL, NEED NH4 AT SRATIO=2.0
         SULRATW = GETASR(WAER(2), RHI)     ! AEROSOL SULFATE RATIO
C      ELSE
C         SULRATW = 2.0D0                    ! DRY AEROSOL SULFATE RATIO
C      ENDIF
      SULRAT  = WAER(3)/WAER(2)         ! SULFATE RATIO
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (SULRATW.LE.SULRAT) THEN
C
         SCASE = 'S2'
C         WAERP = WAER
         CALL CALCS2!(WAERP,GASP,AERLIQP)               ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
C      W(2) = WAER(2)
C      W(3) = WAER(3)
C
         SCASE = 'B4R'
C         WAERP = WAER
         CALL CALCB4R!(WAERP,GASP,AERLIQP)            ! Only liquid (metastable)
         SCASE = 'B4R'
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C      W(2) = WAER(2)
C      W(3) = WAER(3)
C
         SCASE = 'C2R'
C         WAERP = WAER
         CALL CALCC2R!(WAERP,GASP,AERLIQP)                 ! Only liquid (metastable)
         SCASE = 'C2R'
C
      ENDIF
      RETURN
C
C *** END OF SUBROUTINE ISRP1R *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP2R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF 
C     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP2R (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ
      INTEGER   ADJCASE
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DOUBLE PRECISION WTP(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
C      TRYLIQ = .TRUE.             ! Assume liquid phase, sulfate poor limit 
C
      CALL INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRATW = GETASR(WAER(2), RHI)     ! LIMITING SULFATE RATIO
      SULRAT = WAER(3)/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (SULRATW.LE.SULRAT) THEN                
         SCASE = 'N3'
C         WAERP = WAER
         CALL CALCN3!(WAERP, GASP, AERLIQP)      ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
C     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
C     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE 
C     AEROSOL EQUILIBRIUM.
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
C         ADJCASE = 2
C      W(2) = WAER(2)     ! Now included within subroutine
C      W(3) = WAER(3)
C      W(4) = WAER(4)
C
         SCASE = 'E4R'
C         WAERP = WAER
         CALL CALCE4R!(WAERP, GASP, AERLIQP)         ! Only liquid (metastable)
         SCASE = 'E4R'
C
C *** SULFATE RICH (FREE ACID)
C
C     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
C     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE 
C     AEROSOL EQUILIBRIUM.
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C         ADJCASE = 3
C      W(2) = WAER(2)        ! Now included within subroutine
C      W(3) = WAER(3)
C      W(4) = WAER(4)
C
         SCASE = 'F2R'
C         WAERP = WAER
         CALL CALCF2R!(WAERP, GASP, AERLIQP)     ! Only liquid (metastable)
         SCASE = 'F2R'
CC
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE ISRP2R *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP3R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM 
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP3R (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DIMENSION WI(NCOMP)
      DOUBLE PRECISION WAERP(NCOMP), GASP(3), AERLIQP(NIONS+NGASAQ+2)
      DOUBLE PRECISION WTP(NCOMP)
      LOGICAL   TRYLIQ
ccC
ccC *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
ccC
cc      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
cc      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** INITIALIZE ALL VARIABLES ******************************************
C
C      TRYLIQ = .TRUE.             ! Use liquid phase sulfate poor limit 
C
      CALL ISOINIT3 (WI, RHI, TEMPI) ! COMMON block variables
ccC
ccC *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
ccC
cc      REST = 2.D0*WAER(2) + WAER(4) + WAER(5) 
cc      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
cc         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
cc         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
cc      ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
C      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN  ! ** WET AEROSOL
C      IF (RH.GE.DRNH4NO3) THEN  ! ** WET AEROSOL
         FRSO4   = WAER(2) - WAER(1)/2.0D0     ! SULFATE UNBOUND BY SODIUM
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR(FRSO4, RHI)          ! SULFATE RATIO FOR NH4+
         SULRATW = (WAER(1)+FRSO4*SRI)/WAER(2) ! LIMITING SULFATE RATIO
         SULRATW = MIN (SULRATW, 2.0D0)
C      ELSE
C         SULRATW = 2.0D0                     ! ** DRY AEROSOL
C      ENDIF
      SULRAT = (WAER(1)+WAER(3))/WAER(2)
      SODRAT = WAER(1)/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR ; SODIUM POOR
C
C      WRITE(*,*) 'SULRAT ',SULRAT,' SULRATW ',SULRATW
C      WRITE(*,*) 'SODRAT ',SODRAT
C      PAUSE
      IF (SULRATW.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
C
         SCASE = 'Q5'
C         WAERP = WAER
         CALL CALCQ5!(WAERP, GASP, AERLIQP)         ! Only liquid (metastable)
         SCASE = 'Q5'
C
C *** SULFATE POOR ; SODIUM RICH
C
      ELSE IF (SULRAT.GE.SULRATW .AND. SODRAT.GE.2.0) THEN                
C
         SCASE = 'R6'
C         WAERP = WAER
         CALL CALCR6!(WAERP, GASP, AERLIQP)  ! Only liquid (metastable)
         SCASE = 'R6'
C
C *** SULFATE RICH (NO ACID) 
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
C         W = WAER
C
         SCASE = 'I6R'
C         WAERP = WAER
         CALL CALCI6R!(WAERP, GASP, AERLIQP)         ! Only liquid (metastable)
         SCASE = 'I6R'
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C         W = WAER
C
         SCASE = 'J3R'
C         WAERP = WAER
         CALL CALCJ3R!(WAERP, GASP, AERLIQP)       ! Only liquid (metastable)
         SCASE = 'J3R'
C
      ENDIF
C
C *** IF AFTER CALCULATIONS, SULRATW < SULRAT < 2.0  
C                            and WATER = 0          => SULFATE RICH CASE.
C
C      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0  
C     &                      .AND. WATER.LE.TINY) THEN
C          TRYLIQ = .FALSE.
C          GOTO 10
C      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE ISRP3R *****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCS2
C *** CASE S2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCS2!(WAERP,GASP,AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION NH4I, NH3GI, NH3AQ
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DOUBLE PRECISION WTP(NCOMP)
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU   =.TRUE.     ! Outer loop activity calculation flag
      FRST     =.TRUE.
      CALAIN   =.TRUE.
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      WAER = WAERP
C
C *** CALCULATE WATER CONTENT *****************************************
C
      MOLALR(4)= MIN(WAER(2), 0.5d0*WAER(3))
      WATER    = MOLALR(4)/M0(4)  ! ZSR correlation
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      I = 1
      DO WHILE ((I .LE. NSWEEP).AND.(CALAIN))
CC         A21  = XK21*WATER*R*TEMP
         A2   = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
         AKW  = XKW *RH*WATER*WATER
C
         NH4I = WAER(3)
         SO4I = WAER(2)
         HSO4I= ZERO
C
         CALL CALCPH (2.D0*SO4I - NH4I, HI, OHI)    ! Get pH
C
         NH3AQ = ZERO                               ! AMMONIA EQUILIBRIUM
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ (NH4I, OHI, DEL)
            NH4I  = MAX (NH4I-DEL, ZERO) 
            OHI   = MAX (OHI -DEL, TINY)
            NH3AQ = DEL
            HI    = AKW/OHI
         ENDIF
C
         CALL CALCHS4 (HI, SO4I, ZERO, DEL)         ! SULFATE EQUILIBRIUM
         SO4I  = SO4I - DEL
         HI    = HI   - DEL
         HSO4I = DEL
C
         NH3GI = NH4I/HI/A2   !    NH3AQ/A21
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL(1) = HI
         MOLAL(3) = NH4I
         MOLAL(5) = SO4I
         MOLAL(6) = HSO4I
         COH      = OHI
         GASAQ(1) = NH3AQ
         GNH3     = NH3GI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         CALL CALCACT3     
         I = I + 1
      ENDDO 
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCS2')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
C *** END OF SUBROUTINE CALCS2 ****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB4R
C *** CASE B4R
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
      SUBROUTINE CALCB4R!(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION X,Y, SO4I, HSO4I, BB, CC, DD
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DOUBLE PRECISION WTP(NCOMP) 
      DIMENSION WI(NCOMP)
      INTEGER          I
C
C *** SOLVE EQUATIONS **************************************************
C
C      WAER = WAERP
C
      W(2) = WAER(2)
      W(3) = WAER(3)
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
         CALL CALCACT3  !*** slc.11.2009 moved to beginning of loop
         I = I + 1
C20    CONTINUE
      ENDDO
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCB4')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C
      CALL CALCNH3P          ! Compute NH3(g)
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
C30    RETURN
C
C *** END OF SUBROUTINE CALCB4 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCC2R
C *** CASE C2R 
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
      SUBROUTINE CALCC2R!(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION LAMDA, KAPA, PSI, PARM
      DOUBLE PRECISION BB, CC
      DOUBLE PRECISION WP, GAS(3), AERLIQ(7)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DIMENSION WI(NCOMP) 
      INTEGER          I
C
C      WAER = WAERP
C
      W(2) = WAER(2)
      W(3) = WAER(3)
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
C
      CALL CALCNH3P
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C 
C30    RETURN
C    
C *** END OF SUBROUTINE CALCC2R ****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCN3
C *** CASE N3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCN3!(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION AML5, GG, DEL
      DOUBLE PRECISION NH4I, NO3I, SO4I, HSO4I, HI, OHI, NH3AQ, NO3AQ
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DOUBLE PRECISION WTP(NCOMP) 
      LOGICAL TST
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      WAER = WAERP
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU =.TRUE.              ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
      TST    =.TRUE.
C
C *** AEROSOL WATER CONTENT
C
      MOLALR(4) = MIN(WAER(2),0.5d0*WAER(3))       ! (NH4)2SO4
      AML5      = MAX(WAER(3)-2.D0*MOLALR(4),ZERO) ! "free" NH4
      MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO)     ! NH4NO3=MIN("free",NO3)
      WATER     = MOLALR(4)/M0(4) + MOLALR(5)/M0(5)
      WATER     = MAX(WATER, TINY)
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      I = 1
      DO WHILE ((I .LE. NSWEEP) .AND. TST)
C
         IF (I.GT.1) CALL CALCACT3     
C
         A2    = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
CC         A21   = XK21*WATER*R*TEMP
         A3    = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4    = XK7*(WATER/GAMA(4))**3.0
         AKW   = XKW *RH*WATER*WATER
C
C ION CONCENTRATIONS
C
         NH4I  = WAER(3)
         NO3I  = WAER(4)
         SO4I  = WAER(2)
         HSO4I = ZERO
C
         GG    = 2.D0*SO4I + NO3I - NH4I
         CALL CALCPH (GG, HI, OHI)
C
C AMMONIA ASSOCIATION EQUILIBRIUM
C
         NH3AQ = ZERO
         NO3AQ = ZERO
         GG    = 2.D0*SO4I + NO3I - NH4I
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
         ELSE
            HI    = ZERO
            CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) ! HNO3
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
         ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
         MOLAL (1) = HI
         MOLAL (3) = NH4I
         MOLAL (5) = SO4I
         MOLAL (6) = HSO4I
         MOLAL (7) = NO3I
         COH       = OHI
C
         CNH42S4   = ZERO
         CNH4NO3   = ZERO
C
         GASAQ(1)  = NH3AQ
         GASAQ(3)  = NO3AQ
C
         GHNO3     = HI*NO3I/A3
         GNH3      = NH4I/HI/A2   !   NH3AQ/A21 
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            TST = .TRUE.
         ELSE
            TST = .FALSE. 
         ENDIF
         I = I+1
      ENDDO
C      
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCN3')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
C *** RETURN ***********************************************************
C
      RETURN
C
C *** END OF SUBROUTINE CALCN3 *****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB4R
C *** CASE B4R
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
      SUBROUTINE CALCE4R!(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION X,Y, SO4I, HSO4I, BB, CC, DD
      DOUBLE PRECISION GAS(3), AERLIQ(NIONS+NGASAQ+2)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DOUBLE PRECISION WI(NCOMP)
      INTEGER          I
C
C *** SOLVE EQUATIONS **************************************************
C
C      WAER = WAERP
C
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)
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
      DO WHILE ((I.LE.NSWEEP).AND.(CALAIN))
C        
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
C         CALL CALCMR        ! Water content   - now calculated here
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
         I = I + 1
         CALL CALCACT3  
C
      ENDDO
C      
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCE4R')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C
      MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
      MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
      CALL CALCNAP                   ! HNO3, NH3 dissolved
      CALL CALCNH3P
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
C30    RETURN
C
C *** END OF SUBROUTINE CALCB4 ******************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCF2R
C *** CASE C2R 
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
      SUBROUTINE CALCF2R!(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION LAMDA, KAPA, PSI, PARM
      DOUBLE PRECISION BB, CC
      DOUBLE PRECISION GAS(3), AERLIQ(NIONS+NGASAQ+2)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DIMENSION WI(NCOMP)
      INTEGER          I
C
C      WAER = WAERP
C
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)
C
      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** SOLVE EQUATIONS **************************************************
C
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION
C     
      I = 1
      DO WHILE ((I.LE.NSWEEP).AND.(CALAIN))
C         IF (I.GT.1) CALL CALCACT3 
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.d0
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
         MOLALR(9) = MOLAL(3)                              ! NH4HSO4
         MOLALR(7) = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3),ZERO)  ! H2SO4
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
C
      ENDDO
C      WRITE(*,*) 'I',I,'MOLAL(1,5,6)',MOLAL(1),MOLAL(5),MOLAL(6)
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0001, 'CALCF2R')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C
C *** Add the NO3 to the solution now and calculate partitioning.
C
      MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
      MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
C
      CALL CALCNAP                   ! HNO3, NH3 dissolved
      CALL CALCNH3P
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C    
C *** END OF SUBROUTINE CALCC2R ****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCQ5
C *** CASE Q5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
C     2. LIQUID AND SOLID PHASES ARE POSSIBLE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCQ5(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
C
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ
      DOUBLE PRECISION GAS(3), AERLIQ(NIONS+NGASAQ+2)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3)
      DOUBLE PRECISION WAERP(NCOMP)
      LOGICAL TST
C
C *** SETUP PARAMETERS ************************************************
C
C      WAER = WAERP
C
      FRST    =.TRUE.
      CALAIN  =.TRUE. 
      CALAOU  =.TRUE.
      TST     =.TRUE.
C
C *** CALCULATE INITIAL SOLUTION ***************************************
C
C      CALL CALCQ1A - code inserted below
C
C *** CALCULATE SOLIDS **************************************************
C
      CNA2SO4 = 0.5d0*WAER(1)
      FRSO4   = MAX (WAER(2)-CNA2SO4, ZERO)
C
      CNH42S4 = MAX (MIN(FRSO4,0.5d0*WAER(3)), TINY)
      FRNH3   = MAX (WAER(3)-2.D0*CNH42S4, ZERO)
C
      CNH4NO3 = MIN (FRNH3, WAER(4))
CCC      FRNO3   = MAX (WAER(4)-CNH4NO3, ZERO)
      FRNH3   = MAX (FRNH3-CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRNH3, WAER(5))
CCC      FRCL    = MAX (WAER(5)-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      PSI1   = CNA2SO4      ! SALTS DISSOLVED
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
C
C      CALL CALCMR           ! WATER   -  code inserted here
C
      MOLALR(2) = PSI1       ! NA2SO4
      MOLALR(4) = PSI6       ! (NH4)2SO4
      MOLALR(5) = PSI5       ! NH4NO3
      MOLALR(6) = PSI4       ! NH4CL
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
      WATER = ZERO
      DO I=1,NPAIR
         WATER = WATER + MOLALR(I)/M0(I)
      ENDDO
      WATER = MAX(WATER, TINY)
C
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      I = 1
      DO WHILE ((I.LE.NSWEEP) .AND. TST) 
         AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
C
C ION CONCENTRATIONS
C
         NAI    = WAER(1)
         SO4I   = WAER(2)
         NH4I   = WAER(3)
         NO3I   = WAER(4)
         CLI    = WAER(5)
C
C SOLUTION ACIDIC OR BASIC?
C
         GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
         IF (GG.GT.TINY) THEN                        ! H+ in excess
            BB =-GG
            CC =-AKW
            DD = BB*BB - 4.D0*CC
            HI = 0.5D0*(-BB + SQRT(DD))
            OHI= AKW/HI
         ELSE                                        ! OH- in excess
            BB = GG
            CC =-AKW
            DD = BB*BB - 4.D0*CC
            OHI= 0.5D0*(-BB + SQRT(DD))
            HI = AKW/OHI
         ENDIF
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
            HSO4I = ZERO
         ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
               IF (GGNO3.GT.TINY) THEN
                  IF (GGCL.LE.TINY) HI = ZERO
                  CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
               ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
         ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
         MOLAL(1) = HI
         MOLAL(2) = NAI
         MOLAL(3) = NH4I
         MOLAL(4) = CLI
         MOLAL(5) = SO4I
         MOLAL(6) = HSO4I
         MOLAL(7) = NO3I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         CALL CALCACT3
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            TST = .TRUE.
         ELSE
            TST = .FALSE.
         ENDIF
         I = I+1
      ENDDO
C
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0002, 'CALCQ5')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
C
      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4
C
      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ
C
      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = ZERO
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
      RETURN
C
C *** END OF SUBROUTINE CALCQ5 ******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCR6
C *** CASE R6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCR6!(WAERP, GASP, AERLIQP)  
      INCLUDE 'isrpia_b.inc'
C
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ
      DOUBLE PRECISION GAS(3), AERLIQ(NIONS+NGASAQ+2), WTP(NCOMP)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      LOGICAL TST
C
C      WAER = WAERP
C
      FRST   = .TRUE.
      CALAIN = .TRUE. 
      CALAOU = .TRUE. 
      TST    = .TRUE.
C
C *** SETUP PARAMETERS ************************************************
C
C      CALL CALCR1A     ! Code inserted below
C
C *** CALCULATE SOLIDS **************************************************
C
      CNA2SO4 = WAER(2)
      FRNA    = MAX (WAER(1)-2*CNA2SO4, ZERO)
C
      CNH42S4 = ZERO
C
      CNANO3  = MIN (FRNA, WAER(4))
      FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
      FRNA    = MAX (FRNA-CNANO3, ZERO)
C
      CNACL   = MIN (FRNA, WAER(5))
      FRCL    = MAX (WAER(5)-CNACL, ZERO)
      FRNA    = MAX (FRNA-CNACL, ZERO)
C
      CNH4NO3 = MIN (FRNO3, WAER(3))
      FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
      FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRCL, FRNH3)
      FRCL    = MAX (FRCL-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
C
      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
C
C *** CALCULATE WATER **************************************************
C
C      CALL CALCMR    ! Code inserted below
C
      MOLALR(1) = PSI3                                  ! NACL 
      MOLALR(2) = PSI1                                  ! NA2SO4
      MOLALR(3) = PSI2                                  ! NANO3
      MOLALR(4) = ZERO                                  ! (NH4)2SO4
      MOLALR(5) = PSI5                                  ! NH4NO3
      MOLALR(6) = PSI4                                  ! NH4CL
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
      WATER = ZERO
      DO I=1,NPAIR
         WATER = WATER + MOLALR(I)/M0(I)
      ENDDO
      WATER = MAX(WATER, TINY)
C
C *** SETUP LIQUID CONCENTRATIONS **************************************
C
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      I = 1
      DO WHILE ((I.LE.NSWEEP) .AND. TST) 
         AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+      
C
         NAI    = WAER(1)
         SO4I   = WAER(2)
         NH4I   = WAER(3)
         NO3I   = WAER(4)
         CLI    = WAER(5)
C
C SOLUTION ACIDIC OR BASIC?
C
         GG  = 2.D0*WAER(2) + NO3I + CLI - NAI - NH4I
         IF (GG.GT.TINY) THEN                        ! H+ in excess
            BB =-GG
            CC =-AKW
            DD = BB*BB - 4.D0*CC
            HI = 0.5D0*(-BB + SQRT(DD))
            OHI= AKW/HI
         ELSE                                        ! OH- in excess
            BB = GG
            CC =-AKW
            DD = BB*BB - 4.D0*CC
            OHI= 0.5D0*(-BB + SQRT(DD))
            HI = AKW/OHI
         ENDIF
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
         ELSE
            GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
            GGCL  = MAX(GG-GGNO3, ZERO)
            IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
            IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
               CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
            ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
         ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
         MOLAL(1) = HI
         MOLAL(2) = NAI
         MOLAL(3) = NH4I
         MOLAL(4) = CLI
         MOLAL(5) = SO4I
         MOLAL(6) = HSO4I
         MOLAL(7) = NO3I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         CALL CALCACT3
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            TST = .TRUE. 
         ELSE
            TST = .FALSE. 
         ENDIF
         I = I + 1
      ENDDO
C
      IF (CALAIN .AND. (I .GT. (NSWEEP+1))) THEN
         CALL PUSHERR (0002, 'CALCR6')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
C
      GNH3     = NH4I/HI/A2
      GHNO3    = HI*NO3I/A3
      GHCL     = HI*CLI /A4
C      WRITE(*,*) 'A2: ',A2
C      WRITE(*,*) 'A3: ',A3
C      WRITE(*,*) 'A4: ',A4
C
      GASAQ(1) = NH3AQ
      GASAQ(2) = CLAQ
      GASAQ(3) = NO3AQ
C
      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4CL   = ZERO
      CNACL    = ZERO
      CNANO3   = ZERO
      CNA2SO4  = ZERO 
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                ! Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
      RETURN
C
C *** END OF SUBROUTINE CALCR6 ******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCI6R
C *** CASE I6R
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
      SUBROUTINE CALCI6R!(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION GAS(3), AERLIQ(NIONS+NGASAQ+2)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DOUBLE PRECISION WP(5), WTP(NCOMP)
C         
C      WAER = WAERP
      W = WAER
C
C *** FIND DRY COMPOSITION ********************************************
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
         A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.D0
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
         CALL PUSHERR (0001, 'CALCI6R')    ! WARNING ERROR: NO SOLUTION
      ENDIF
C      CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
C
C     Inserting subroutine here (slc.3.2012 - for adjoint development)
C
C *** IS THERE A LIQUID PHASE? ******************************************
C
      IF ((WATER).LE.TINY) RETURN
C
C *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
C
      A3       = XK3*R*TEMP*(WATER/GAMA(11))**2.0
      A4       = XK4*R*TEMP*(WATER/GAMA(10))**2.0
      MOLAL(1) = MOLAL(1) + WAER(4) + WAER(5)  ! H+ increases because NO3, Cl are added.
C
C *** CALCULATE CONCENTRATIONS ******************************************
C *** ASSUME THAT 'DELT' FROM HNO3 >> 'DELT' FROM HCL
C
      CALL CALCNIAQ (WAER(4), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT 
      MOLAL(7) = WAER(4)  - DELT  ! NO3- = Waer(4) minus any turned into (HNO3aq)
      GASAQ(3) = DELT
C
      CALL CALCCLAQ (WAER(5), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(4) = WAER(5)  - DELT  ! Cl- = Waer(4) minus any turned into (HNO3aq)
      GASAQ(2) = DELT
C
      GHNO3    = MOLAL(1)*MOLAL(7)/A4
      GHCL     = MOLAL(1)*MOLAL(4)/A3
C
      CALL CALCNH3P
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = WAER(1)                C Total gas+aerosol phase
C      WTP(2) = WAER(2)
C      WTP(3) = WAER(3) + GNH3 
C      WTP(4) = WAER(4) + GHNO3
C      WTP(5) = WAER(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
CC
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
      RETURN
C
C *** END OF SUBROUTINE CALCI6R *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCJ3R
C *** CASE J3R
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
      SUBROUTINE CALCJ3R!(WAERP, GASP, AERLIQP)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION LAMDA, KAPA
      DOUBLE PRECISION GAS(3), AERLIQ(NIONS+NGASAQ+2)
      DOUBLE PRECISION AERLIQP(NIONS+NGASAQ+2), GASP(3), WAERP(NCOMP)
      DOUBLE PRECISION WP(5), WTP(NCOMP)
C         
C      WAER = WAERP
      W = WAER
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
C      CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
C
C     Inserting subroutine here (slc.3.2012 - for adjoint development)
C
C *** IS THERE A LIQUID PHASE? ******************************************
C
      IF ((WATER).LE.TINY) RETURN
C
C *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
C
      A3       = XK3*R*TEMP*(WATER/GAMA(11))**2.0
      A4       = XK4*R*TEMP*(WATER/GAMA(10))**2.0
      MOLAL(1) = MOLAL(1) + WAER(4) + WAER(5)  ! H+ increases because NO3, Cl are added.
C
C *** CALCULATE CONCENTRATIONS ******************************************
C *** ASSUME THAT 'DELT' FROM HNO3 >> 'DELT' FROM HCL
C
      CALL CALCNIAQ (WAER(4), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT 
      MOLAL(7) = WAER(4)  - DELT  ! NO3- = Waer(4) minus any turned into (HNO3aq)
      GASAQ(3) = DELT
C
      CALL CALCCLAQ (WAER(5), MOLAL(1)+MOLAL(7)+MOLAL(4), DELT)
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(4) = WAER(5)  - DELT  ! Cl- = Waer(4) minus any turned into (HNO3aq)
      GASAQ(2) = DELT
C
      GHNO3    = MOLAL(1)*MOLAL(7)/A4
      GHCL     = MOLAL(1)*MOLAL(4)/A3
C
      CALL CALCNH3P
C
C *** CALCULATE TOTAL AMOUNT OF SPECIES WITH GAS PHASE ****************
C
C      WTP(1) = W(1)                C Total gas+aerosol phase
C      WTP(2) = W(2)
C      WTP(3) = W(3) + GNH3 
C      WTP(4) = W(4) + GHNO3
C      WTP(5) = W(5) + GHCL
C
C *** ADJOINT DEVELOPMENT *********************************************
C
C      GASP(1) = GNH3                ! Gaseous aerosol species
C      GASP(2) = GHNO3
C      GASP(3) = GHCL
C
C      DO I=1,NIONS              ! Liquid aerosol species
C         AERLIQP(I) = MOLAL(I)
C      ENDDO
C      DO I=1,NGASAQ
C         AERLIQP(NIONS+1+I) = GASAQ(I)
C      ENDDO
C      AERLIQP(NIONS+1)        = WATER*1.0D3/18.0D0
C      AERLIQP(NIONS+NGASAQ+2) = COH
C
      RETURN
C
C *** END OF SUBROUTINE CALCJ3R *****************************************
C
      END
C


