C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISOROPIA
C *** THIS SUBROUTINE IS THE MASTER ROUTINE FOR THE ISORROPIA
C     THERMODYNAMIC EQUILIBRIUM AEROSOL MODEL (VERSION 1.1 and above)
C
C ======================== ARGUMENTS / USAGE ===========================
C
C  INPUT:
C  1. [WI] 
C     DOUBLE PRECISION array of length [5].
C     Concentrations, expressed in moles/m3. Depending on the type of
C     problem solved (specified in CNTRL(1)), WI contains either 
C     GAS+AEROSOL or AEROSOL only concentratios.
C     WI(1) - sodium
C     WI(2) - sulfate
C     WI(3) - ammonium
C     WI(4) - nitrate
C     WI(5) - chloride
C
C  2. [RHI] 
C     DOUBLE PRECISION variable.  
C     Ambient relative humidity expressed on a (0,1) scale.
C
C  3. [TEMPI]
C     DOUBLE PRECISION variable. 
C     Ambient temperature expressed in Kelvins. 
C
C  4. [CNTRL]
C     DOUBLE PRECISION array of length [2].
C     Parameters that control the type of problem solved.
C
C     CNTRL(1): Defines the type of problem solved.
C     0 - Forward problem is solved. In this case, array WI contains 
C         GAS and AEROSOL concentrations together.
C     1 - Reverse problem is solved. In this case, array WI contains
C         AEROSOL concentrations only.
C
C     CNTRL(2): Defines the state of the aerosol
C     0 - The aerosol can have both solid+liquid phases (deliquescent)
C     1 - The aerosol is in only liquid state (metastable aerosol)
C
C  OUTPUT:
C  1. [WT] 
C     DOUBLE PRECISION array of length [5].
C     Total concentrations (GAS+AEROSOL) of species, expressed in moles/m3. 
C     If the foreward probelm is solved (CNTRL(1)=0), array WT is 
C     identical to array WI.
C     WT(1) - total sodium
C     WT(2) - total sulfate
C     WT(3) - total ammonium
C     WT(4) - total nitrate
C     WT(5) - total chloride
C
C  2. [GAS]
C     DOUBLE PRECISION array of length [03]. 
C     Gaseous species concentrations, expressed in moles/m3. 
C     GAS(1) - NH3
C     GAS(2) - HNO3
C     GAS(3) - HCl 
C
C  3. [AERLIQ]
C     DOUBLE PRECISION array of length [11]. 
C     Liquid aerosol species concentrations, expressed in moles/m3. 
C     AERLIQ(01) - H+(aq)          
C     AERLIQ(02) - Na+(aq)         
C     AERLIQ(03) - NH4+(aq)
C     AERLIQ(04) - Cl-(aq)         
C     AERLIQ(05) - SO4--(aq)       
C     AERLIQ(06) - HSO4-(aq)       
C     AERLIQ(07) - NO3-(aq)        
C     AERLIQ(08) - H2O             
C     AERLIQ(09) - NH3(aq) (undissociated)
C     AERLIQ(10) - HNCl(aq) (undissociated)
C     AERLIQ(11) - HNO3(aq) (undissociated)
C     AERLIQ(12) - OH-(aq)
C
C  4. [AERSLD]
C     DOUBLE PRECISION array of length [09]. 
C     Solid aerosol species concentrations, expressed in moles/m3. 
C     AERSLD(01) - NaNO3(s)
C     AERSLD(02) - NH4NO3(s)
C     AERSLD(03) - NaCl(s)         
C     AERSLD(04) - NH4Cl(s)
C     AERSLD(05) - Na2SO4(s)       
C     AERSLD(06) - (NH4)2SO4(s)
C     AERSLD(07) - NaHSO4(s)
C     AERSLD(08) - NH4HSO4(s)
C     AERSLD(09) - (NH4)4H(SO4)2(s)
C
C  5. [SCASI]
C     CHARACTER(15) variable.
C     Returns the subcase which the input corresponds to.
C
C  6. [OTHER]
C     DOUBLE PRECISION array of length [6].
C     Returns solution information.
C
C     OTHER(1): Shows if aerosol water exists.
C     0 - Aerosol is WET
C     1 - Aerosol is DRY
C
C     OTHER(2): Aerosol Sulfate ratio, defined as (in moles/m3) :
C               (total ammonia + total Na) / (total sulfate)
C
C     OTHER(3): Sulfate ratio based on aerosol properties that defines 
C               a sulfate poor system:
C               (aerosol ammonia + aerosol Na) / (aerosol sulfate)
C           
C     OTHER(4): Aerosol sodium ratio, defined as (in moles/m3) :
C               (total Na) / (total sulfate)
C      
C     OTHER(5): Ionic strength of the aqueous aerosol (if it exists).
C      
C     OTHER(6): Total number of calls to the activity coefficient 
C               calculation subroutine.
C 
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISOROPIA (WI, RHI, TEMPI,  CNTRL, 
     &                     WT, GAS, AERLIQ, AERSLD, SCASI, OTHER, 
     &                     TRUSTISO)
      INCLUDE 'isrpia_b.inc'
      PARAMETER (NCTRL=2,NOTHER=6)
      CHARACTER SCASI*15
      DOUBLE PRECISION CNTRL, AERSLD, OTHER
      LOGICAL   TRUSTISO
      DOUBLE PRECISION WI, RHI, TEMPI, WT, GAS, AERLIQ
      DIMENSION WI(NCOMP), WT(NCOMP),   GAS(NGASAQ),  AERSLD(NSLDS), 
     &          AERLIQ(NIONS+NGASAQ+2), CNTRL(NCTRL), OTHER(NOTHER)
	  INTEGER   ERRSTKI(25)
      CHARACTER(40) ERRMSGI(25)
C
C *** PROBLEM TYPE (0=FOREWARD, 1=REVERSE) ******************************
C
      IPROB   = NINT(CNTRL(1))
C
C *** AEROSOL STATE (0=SOLID+LIQUID, 1=METASTABLE) **********************
C
      METSTBL = NINT(CNTRL(2))
C
C *** SOLVE FOREWARD PROBLEM ********************************************
C
50    IF (IPROB.EQ.0) THEN
         IF ((WI(1)+WI(2)+WI(3)+WI(4)+WI(5)) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF ((WI(1)+WI(4)+WI(5)) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1F (WI, RHI, TEMPI)
         ELSE IF ((WI(1)+WI(5)) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2F (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3F (WI, RHI, TEMPI)
         ENDIF
C
C *** SOLVE REVERSE PROBLEM *********************************************
C
      ELSE
         IF ((WI(1)+WI(2)+WI(3)+WI(4)+WI(5)) .LE.TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF ((WI(1)+WI(4)+WI(5)) .LE. TINY) THEN        ! Na,Cl,NO3=0
            CALL ISRP1R (WI, RHI, TEMPI)
         ELSE IF ((WI(1)+WI(5)) .LE. TINY) THEN              ! Na,Cl=0
            CALL ISRP2R (WI, RHI, TEMPI)
         ELSE
            CALL ISRP3R (WI, RHI, TEMPI)
         ENDIF
      ENDIF
C
C *** ADJUST MASS BALANCE ***********************************************
C
      IF (NADJ.EQ.1) CALL ADJUST (WI)
ccC
ccC *** IF METASTABLE AND NO WATER - RESOLVE AS NORMAL ********************
ccC
cc      IF (WATER.LE.TINY .AND. METSTBL.EQ.1) THEN
cc         METSTBL = 0
cc         GOTO 50
cc      ENDIF
C
C *** SAVE RESULTS TO ARRAYS (units = mole/m3) ****************************
C
      GAS(1) = GNH3                ! Gaseous aerosol species
      GAS(2) = GHNO3
      GAS(3) = GHCL
C
C      WRITE(*,*) 'MOLAL',MOLAL
      DO I=1,NIONS              ! Liquid aerosol species
         AERLIQ(I) = MOLAL(I)
      ENDDO 
      DO I=1,NGASAQ
         AERLIQ(NIONS+1+I) = GASAQ(I)
      ENDDO 
      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0
      AERLIQ(NIONS+NGASAQ+2) = COH
C
      AERSLD(1) = CNANO3           ! Solid aerosol species
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNACL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNA2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNAHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC
C
      IF(WATER.LE.TINY) THEN       ! Dry flag
        OTHER(1) = 1.d0
      ELSE
        OTHER(1) = 0.d0
      ENDIF
C
      OTHER(2) = SULRAT            ! Other stuff
      OTHER(3) = SULRATW
      OTHER(4) = SODRAT
      OTHER(5) = IONIC
      OTHER(6) = ICLACT
C
      SCASI = SCASE
C
      WT(1) = WI(1)                ! Total gas+aerosol phase
      WT(2) = WI(2)
      WT(3) = WI(3) 
      WT(4) = WI(4)
      WT(5) = WI(5)
      IF (IPROB.GT.0 .AND. WATER.GT.TINY) THEN 
         WT(3) = WT(3) + GNH3 
         WT(4) = WT(4) + GHNO3
         WT(5) = WT(5) + GHCL
      ENDIF
C      WRITE(*,*) ''
C      WRITE(*,*) '******** Results and Sensitivities in CVM ***********'
C      WRITE(*,*) 'CVM: AERLIQ ', AERLIQ
C      WRITE(*,*) 'AERLIQ Perturbation', IMAG(AERLIQ)/1.d-18
C      WRITE(*,*) 'GAS ', GAS
C      WRITE(*,*) 'GAS Perturbation', IMAG(GAS)/1.d-18
C      WRITE(*,*) '***** End of Results and Sensitivities in CVM *******'
C      WRITE(*,*) ''
C
C *** Check for errors ****************************************************
C
      TRUSTISO = .TRUE.
      CALL ISERRINF (ERRSTK, ERRMSG, NOFER, STKOFL) ! Obtain error stack
      IF (NOFER.GT.0) TRUSTISO = .FALSE.   ! Errors found
C
      RETURN
C
C *** END OF SUBROUTINE ISOROPIA ******************************************
C
      END


C=======================================================================
C
C *** ANISORROPIA CODE
C *** SUBROUTINE ISOROPIA_B
C *** THIS SUBROUTINE IS THE MASTER ROUTINE FOR THE ADJOINT OF ISORROPIA
C     THERMODYNAMIC EQUILIBRIUM AEROSOL MODEL (VERSION 1.1 and above)
C
C=======================================================================
C
      SUBROUTINE ISOROPIA_B(WI, WPB, RHI, TEMPI,  CNTRL, 
     &               WT,  GAS, GASb, AERLIQ, AERLIQb, AERSLD, 
     &               SCASI, OTHER, TRUSTISO)
      INCLUDE 'isrpia_b.inc'
      PARAMETER (NCTRL=2,NOTHER=6)
      CHARACTER(15)    SCASI
      LOGICAL          TRUSTISO, TRYLIQ
      DOUBLE PRECISION wp(ncomp), aerliq, gas
      DOUBLE PRECISION wpb(ncomp)
      DOUBLE PRECISION WTORIG, GASORIG, AERLIQORIG
      DOUBLE PRECISION wi, RHI, aerliqb, aerliqallb, gasb
      DOUBLE PRECISION aerliqab(nions), aerliqbb(nions), aerliqcb(nions) 
      DOUBLE PRECISION aerliqdb(nions), aerliqeb(nions), aerliqfb(nions)
C      DIMENSION WTORIG(NCOMP),GASORIG(NGASAQ),AERLIQORIG(NIONS+NGASAQ+2)
      DIMENSION WI(NCOMP), WT(NCOMP),   GAS(NGASAQ),  AERSLD(NSLDS), 
     &          AERLIQ(NIONS+NGASAQ+2), CNTRL(NCTRL), OTHER(NOTHER)
      DIMENSION gasb(ngasaq), aerliqb(nions)
      DIMENSION aerliqallb(nions+ngasaq+2)
      INTEGER   ERRSTKI(25)
      CHARACTER(40) ERRMSGI(25)
C
C *** PROBLEM TYPE (0=FOREWARD, 1=REVERSE) ******************************
C
      IPROB   = NINT(CNTRL(1))
C
C *** AEROSOL STATE (0=SOLID+LIQUID, 1=METASTABLE) **********************
C
      METSTBL = NINT(CNTRL(2))
C
C *** SOLVE FORWARD PROBLEM ********************************************
C
C      WRITE(*,*) 'AERLIQb',AERLIQb     
      IF (IPROB.EQ.0) THEN
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI,RHI,TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
C           CALL ISRP1F (WI, RHI, TEMPI)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
            CALL INIT1 (WI, RHI, TEMPI)
C            WP = W
C
C *** CALCULATE SULFATE RATIO TO SEND TO APPROPRIATE CALC ***************
C
            SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
        IF (2.0.LE.SULRAT) THEN 
           SCASE = 'A2'
           DO j = 1,nions
              aerliqab(j) = aerliqb(j)
           ENDDO
           CALL ISRP1FA_AB(wpb, gasb, aerliqab)
        ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN
           SCASE = 'B4'
           DO j = 1,nions
              aerliqbb(j) = aerliqb(j)
           ENDDO
           CALL CALCB4_BB(wpb, gasb, aerliqbb)
        ELSEIF (SULRAT.LT.1.0) THEN
           SCASE = 'C2'
           DO j = 1,nions
              aerliqcb(j) = aerliqb(j)
           ENDDO
           CALL CALCC2_CB(wpb, gasb, aerliqcb)
        ELSE
           RETURN
        ENDIF
        RETURN
      ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
        CALL INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO TO SEND TO APPROPRIATE CALC ***************
C
        SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
       IF (2.0.LE.SULRAT) THEN 
          SCASE = 'D3'
          DO j = 1,nions
             aerliqdb(j) = aerliqb(j)
          ENDDO
          CALL CALCD3_B(wpb, gasb, aerliqdb)
       ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN
          SCASE = 'E4'
          DO j = 1,nions
             aerliqeb(j) = aerliqb(j)
          ENDDO
          CALL CALCB4E_EB(wpb, gasb, aerliqeb)
       ELSEIF (SULRAT.LT.1.0) THEN
          SCASE = 'F2'
          DO j = 1,nions
             aerliqfb(j) = aerliqb(j)
          ENDDO
          CALL CALCC2F_FB(wpb, gasb, aerliqfb)
       ELSE
          RETURN
       ENDIF
       RETURN
C
      ELSE IF (WI(1)+WI(5) .GT. TINY) THEN              ! Na,Cl>0
C
C *** SULFATE POOR ; SODIUM POOR
C
         REST = 2.D0*WI(2) + WI(4) + WI(5) 
         IF (WI(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
            WI(1) = (ONE-1D-6)*REST         ! Adjust Na amount
            CALL PUSHERR (0050, 'ISRP3F')  ! Warning error: Na adjusted
         ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
         SULRAT = (WI(1)+WI(3))/WI(2)
         SODRAT = WI(1)/WI(2)
C
         IF (2.0.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
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
            SCASE = 'G5'
            CALL CALCG5_B(wpb, gasb, aerliqb)    ! Only liquid (metastable)
C            
         ELSEIF (SULRAT.GE.2.0 .AND. SODRAT.GE.2.0) THEN                
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
            SCASE = 'H6'
            CALL CALCH6_B(wpb, gasb, aerliqb)             ! Only liquid (metastable)
         ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
            CALL ISOINIT3 (WI, RHI, TEMPI)
            SCASE = 'I6'
            CALL ISRP3F_IB(wpb, gasb, aerliqb)     ! Only liquid (metastable)
         ELSEIF (SULRAT.LT.1.0) THEN             
            CALL ISOINIT3 (WI, RHI, TEMPI)
            SCASE = 'J3'
            CALL ISRP3F_JB(wpb, gasb, aerliqb)     ! Only liquid (metastable)
         ENDIF
      ENDIF
C
C *** SOLVE REVERSE PROBLEM *********************************************
C
       ELSE
C
C *** Ammonium and Sulfate **********************************************
C
         IF (WI(1)+WI(2)+WI(3)+WI(4)+WI(5) .LE. TINY) THEN ! Everything=0
            CALL INIT1 (WI, RHI, TEMPI)
         ELSE IF (WI(1)+WI(4)+WI(5) .LE. TINY) THEN        ! Na,Cl,NO3=0
C            CALL ISRP1R (WI, RHI, TEMPI)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
C
C *** INITIALIZE COMMON BLOCK VARIABLES *********************************
C
            CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
            IF (RH.GE.DRNH42S4) THEN         ! WET AEROSOL, NEED NH4 AT SRATIO=2.0
               SULRATW = GETASR(WAER(2), RHI)     ! AEROSOL SULFATE RATIO
            ELSE
               SULRATW = 2.0D0                    ! DRY AEROSOL SULFATE RATIO
            ENDIF
            SULRAT  = WAER(3)/WAER(2)         ! SULFATE RATIO
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C            WRITE(*,*) 'SULRAT',SULRAT, ' SULRATW ',SULRATW
            IF (2.0.LE.SULRAT) THEN 
               SCASE = 'S2'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
C               WRITE(*,*) 'S2: gasb', gasb
C               WRITE(*,*) 'aerliqallb',aerliqallb
               CALL CALCS2_SB(wpb, gasb, aerliqallb)
            ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN
               SCASE = 'B4R'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCB4R_BRB(wpb, gasb, aerliqallb)
            ELSEIF (SULRAT.LT.1.0) THEN
               SCASE = 'C2R'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCC2R_CRB(wpb, gasb, aerliqallb)
            ELSE
               RETURN
            ENDIF
         ELSE IF (WI(1)+WI(5) .LE. TINY) THEN              ! Na,Cl=0
C            CALL ISRP2R (WI, RHI, TEMPI)
C
C *** *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
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
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCN3_NB(wpb, gasb, aerliqallb)      ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
            ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
               SCASE = 'E4R'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCE4R_ERB(wpb, gasb, aerliqallb)    ! Only liquid (metastable)
C
C *** SULFATE RICH (FREE ACID)
C
            ELSEIF (SULRAT.LT.1.0) THEN             
               SCASE = 'F2R'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCF2R_FRB(wpb, gasb, aerliqallb)    ! Only liquid (metastable)
            ENDIF
         ELSE
C
C            CALL ISRP3R (WI, RHI, TEMPI)  ! Code inserted below
C
            CALL ISOINIT3 (WI, RHI, TEMPI) ! COMMON block variables
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
            FRSO4   = WAER(2) - WAER(1)/2.0D0     ! SULFATE UNBOUND BY SODIUM
            FRSO4   = MAX(FRSO4, TINY)
            SRI     = GETASR(FRSO4, RHI)          ! SULFATE RATIO FOR NH4+
            SULRATW = (WAER(1)+FRSO4*SRI)/WAER(2) ! LIMITING SULFATE RATIO
            SULRATW = MIN (SULRATW, 2.0D0)
            SULRAT = (WAER(1)+WAER(3))/WAER(2)
            SODRAT = WAER(1)/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR ; SODIUM POOR
C
            IF (SULRATW.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
               SCASE = 'Q5'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCQ5_QB(wpb, gasb, aerliqallb)      ! Only liquid (metastable)
C
C *** SULFATE POOR ; SODIUM RICH
C
            ELSEIF (SULRAT.GE.SULRATW .AND. SODRAT.GE.2.0) THEN                
               SCASE = 'R6'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCR6_RB(wpb, gasb, aerliqallb)     ! Only liquid (metastable)
            ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
               SCASE = 'I6'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCI6R_IRB(wpb, gasb, aerliqallb)    ! Only liquid (metastable)
C
C *** SULFATE RICH (FREE ACID)
C
            ELSEIF (SULRAT.LT.1.0) THEN             
               SCASE = 'J3R'
               DO i = 8,12
                  aerliqallb(i) = 0.d0
               ENDDO
               DO i = 1,nions
                  aerliqallb(i) = aerliqb(i)
               ENDDO
               CALL CALCJ3R_JRB(wpb, gasb, aerliqallb)       ! Only liquid (metastable)
            ENDIF
         ENDIF
      ENDIF
C
C *** ADJUST MASS BALANCE ***********************************************
C
      IF (NADJ.EQ.1) CALL ADJUST (WI)
ccC
ccC *** IF METASTABLE AND NO WATER - RESOLVE AS NORMAL ********************
ccC
cc      IF (WATER.LE.TINY .AND. METSTBL.EQ.1) THEN
cc         METSTBL = 0
cc         GOTO 50
cc      ENDIF
C
C *** SAVE RESULTS TO ARRAYS (units = mole/m3) ****************************
C
      GAS(1) = GNH3                ! Gaseous aerosol species
      GAS(2) = GHNO3
      GAS(3) = GHCL
C
      DO I=1,NIONS              ! Liquid aerosol species
         AERLIQ(I) = MOLAL(I)
      ENDDO
      DO I=1,NGASAQ
         AERLIQ(NIONS+1+I) = GASAQ(I)
      ENDDO
      AERLIQ(NIONS+1)        = WATER*1.0D3/18.0D0
      AERLIQ(NIONS+NGASAQ+2) = COH
C
      AERSLD(1) = CNANO3           ! Solid aerosol species
      AERSLD(2) = CNH4NO3
      AERSLD(3) = CNACL
      AERSLD(4) = CNH4CL
      AERSLD(5) = CNA2SO4
      AERSLD(6) = CNH42S4
      AERSLD(7) = CNAHSO4
      AERSLD(8) = CNH4HS4
      AERSLD(9) = CLC
C
      IF(WATER.LE.TINY) THEN       ! Dry flag
        OTHER(1) = 1.d0
      ELSE
        OTHER(1) = 0.d0
      ENDIF
C
      OTHER(2) = SULRAT            ! Other stuff
      OTHER(3) = SULRATW
      OTHER(4) = SODRAT
      OTHER(5) = IONIC
      OTHER(6) = ICLACT
C
      SCASI = SCASE
C
      WT(1) = WI(1)                ! Total gas+aerosol phase
      WT(2) = WI(2)
      WT(3) = WI(3) 
      WT(4) = WI(4)
      WT(5) = WI(5)
      IF (IPROB.GT.0 .AND. WATER.GT.TINY) THEN 
         WT(3) = WT(3) + GNH3 
         WT(4) = WT(4) + GHNO3
         WT(5) = WT(5) + GHCL
      ENDIF
C
C *** Check for errors ****************************************************
C
      TRUSTISO = .TRUE.
      CALL ISERRINF (ERRSTKI, ERRMSGI, NOFER, STKOFL) ! Obtain error stack
      IF (NOFER.GT.0) TRUSTISO = .FALSE.   ! Errors found
C
      RETURN
C
C *** END OF SUBROUTINE ISOROPIA ******************************************
C
      END
C
C
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE SETPARM
C *** THIS SUBROUTINE REDEFINES THE SOLUTION PARAMETERS OF ISORROPIA
C
C ======================== ARGUMENTS / USAGE ===========================
C
C *** NOTE: IF NEGATIVE VALUES ARE GIVEN FOR A PARAMETER, IT IS
C     IGNORED AND THE CURRENT VALUE IS USED INSTEAD.
C 
C  INPUT:
C  1. [WFTYPI] 
C     INTEGER variable.
C     Defines the type of weighting algorithm for the solution in Mutual 
C     Deliquescence Regions (MDR's):
C     0 - MDR's are assumed dry. This is equivalent to the approach 
C         used by SEQUILIB.
C     1 - The solution is assumed "half" dry and "half" wet throughout
C         the MDR.
C     2 - The solution is a relative-humidity weighted mean of the
C         dry and wet solutions (as defined in Nenes et al., 1998)
C
C  2. [IACALCI] 
C     INTEGER variable.
C     Method of activity coefficient calculation:
C     0 - Calculate coefficients during runtime
C     1 - Use precalculated tables
C 
C  3. [EPSI] 
C     DOUBLE PRECITION variable.
C     Defines the convergence criterion for all iterative processes
C     in ISORROPIA, except those for activity coefficient calculations
C     (EPSACTI controls that).
C
C  4. [MAXITI]
C     INTEGER variable.
C     Defines the maximum number of iterations for all iterative 
C     processes in ISORROPIA, except for activity coefficient calculations 
C     (NSWEEPI controls that).
C
C  5. [NSWEEPI]
C     INTEGER variable.
C     Defines the maximum number of iterations for activity coefficient 
C     calculations.
C 
C  6. [EPSACTI] 
C     DOUBLE PRECISION variable.
C     Defines the convergence criterion for activity coefficient 
C     calculations.
C 
C  7. [NDIV] 
C     INTEGER variable.
C     Defines the number of subdivisions needed for the initial root
C     tracking for the bisection method. Usually this parameter should 
C     not be altered, but is included for completeness.
C
C  8. [NADJ]
C     INTEGER variable.
C     Forces the solution obtained to satisfy total mass balance
C     to machine precision
C     0 - No adjustment done (default)
C     1 - Do adjustment
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE SETPARM (WFTYPI,  IACALCI, EPSI, MAXITI, NSWEEPI, 
     &                    EPSACTI, NDIVI, NADJI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION EPSI, EPSACTI
      INTEGER  WFTYPI
C
C *** SETUP SOLUTION PARAMETERS *****************************************
C
      IF (WFTYPI .GE. 0)   WFTYP  = WFTYPI
      IF (IACALCI.GE. 0)   IACALC = IACALCI
      IF (EPSI   .GE.ZERO) EPS    = EPSI
      IF (MAXITI .GT. 0)   MAXIT  = MAXITI
      IF (NSWEEPI.GT. 0)   NSWEEP = NSWEEPI
      IF (EPSACTI.GE.ZERO) EPSACT = EPSACTI
      IF (NDIVI  .GT. 0)   NDIV   = NDIVI
      IF (NADJI  .GE. 0)   NADJ   = NADJI
C
C *** END OF SUBROUTINE SETPARM *****************************************
C
      RETURN
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE GETPARM
C *** THIS SUBROUTINE OBTAINS THE CURRENT VAULES OF THE SOLUTION 
C     PARAMETERS OF ISORROPIA
C
C ======================== ARGUMENTS / USAGE ===========================
C
C *** THE PARAMETERS ARE THOSE OF SUBROUTINE SETPARM
C 
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE GETPARM (WFTYPI,  IACALCI, EPSI, MAXITI, NSWEEPI, 
     &                    EPSACTI, NDIVI, NADJI)
      INCLUDE 'isrpia_b.inc'
      INTEGER  WFTYPI
C
C *** GET SOLUTION PARAMETERS *******************************************
C
      WFTYPI  = WFTYP
      IACALCI = IACALC
      EPSI    = EPS
      MAXITI  = MAXIT
      NSWEEPI = NSWEEP
      EPSACTI = EPSACT
      NDIVI   = NDIV
      NADJI   = NADJ
C
C *** END OF SUBROUTINE GETPARM *****************************************
C
      RETURN
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** BLOCK DATA BLKISO
C *** THIS SUBROUTINE PROVIDES INITIAL (DEFAULT) VALUES TO PROGRAM
C     PARAMETERS VIA DATA STATEMENTS
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C *** ZSR RELATIONSHIP PARAMETERS MODIFIED BY DOUGLAS WALDRON
C *** OCTOBER 2003
C *** BASED ON AIM MODEL III (http://mae.ucdavis.edu/wexler/aim)
C
C=======================================================================
C
      BLOCK DATA BLKISO
      INCLUDE 'isrpia_b.inc'
C
C *** DEFAULT VALUES *************************************************
C
      DATA TEMP/298.0/, R/82.0567D-6/, RH/(0.9D0, 0)/, EPS/1D-6/, 
     &     MAXIT/100/, TINY/1D-20/, GREAT/1D10/, ZERO/0.0D0/, 
     &     ONE/1.0D0/,NSWEEP/4/, TINY2/1D-11/,NDIV/5/
C
      DATA MOLAL/NIONS*0.0D0/, MOLALR/NPAIR*0.0D0/, GAMA/NPAIR*0.1D0/,
     &     GAMOU/NPAIR*1D10/,  GAMIN/NPAIR*1D10/,   CALAIN/.TRUE./,
C     &     CALAOU/.TRUE./,     EPSACT/5D-2/,        ICLACT/0/,
     &     CALAOU/.TRUE./,     EPSACT/1D-10/,        ICLACT/0/,
     &     IACALC/0/,          NADJ/0/,             WFTYP/2/
C
      DATA ERRSTK/NERRMX*0/,   ERRMSG/NERRMX*' '/,  NOFER/0/, 
     &     STKOFL/.FALSE./ 
C
      DATA IPROB/0/, METSTBL/0/
C
      DATA VERSION /'1.7 (03/26/07)'/
C
C *** OTHER PARAMETERS ***********************************************
C
      DATA SMW/58.5,142.,85.0,132.,80.0,53.5,98.0,98.0,115.,63.0,
     &         36.5,120.,247./
     &     IMW/ 1.0,23.0,18.0,35.5,96.0,97.0,63.0/,
     &     WMW/23.0,98.0,17.0,63.0,36.5/
C
      DATA Z /1.0D0,1.0D0,1.0D0,1.0D0,2.0D0,1.0D0,1.0D0/
      DATA ZZ/1,2,1,2,1,1,2,1,1,1,1,1,2/
C      DATA ZZ/1,2,1,2,1,1,2,1,1,1,1,1,2/, Z /1,1,1,1,2,1,1/
C
C *** ZSR RELATIONSHIP PARAMETERS **************************************
C
C awas= ammonium sulfate
C
      DATA AWAS/10*187.72,
     & 158.13,134.41,115.37,100.10, 87.86, 78.00, 70.00, 63.45, 58.02,
     &  53.46,
     &  49.59, 46.26, 43.37, 40.84, 38.59, 36.59, 34.79, 33.16, 31.67,
     &  30.31,
     &  29.07, 27.91, 26.84, 25.84, 24.91, 24.03, 23.21, 22.44, 21.70,
     &  21.01,
     &  20.34, 19.71, 19.11, 18.54, 17.99, 17.46, 16.95, 16.46, 15.99,
     &  15.54,
     &  15.10, 14.67, 14.26, 13.86, 13.47, 13.09, 12.72, 12.36, 12.01,
     &  11.67,
     &  11.33, 11.00, 10.68, 10.37, 10.06,  9.75,  9.45,  9.15,  8.86,
     &   8.57,
     &   8.29,  8.01,  7.73,  7.45,  7.18,  6.91,  6.64,  6.37,  6.10,
     &   5.83,
     &   5.56,  5.29,  5.02,  4.74,  4.47,  4.19,  3.91,  3.63,  3.34,
     &   3.05,
     &   2.75,  2.45,  2.14,  1.83,  1.51,  1.19,  0.87,  0.56,  0.26,
     &  0.1/
C
C awsn= sodium nitrate
C
      DATA AWSN/10*394.54,
     & 338.91,293.01,254.73,222.61,195.56,172.76,153.53,137.32,123.65,
     & 112.08,
     & 102.26, 93.88, 86.68, 80.45, 75.02, 70.24, 66.02, 62.26, 58.89,
     &  55.85,
     &  53.09, 50.57, 48.26, 46.14, 44.17, 42.35, 40.65, 39.06, 37.57,
     &  36.17,
     &  34.85, 33.60, 32.42, 31.29, 30.22, 29.20, 28.22, 27.28, 26.39,
     &  25.52,
     &  24.69, 23.89, 23.12, 22.37, 21.65, 20.94, 20.26, 19.60, 18.96,
     &  18.33,
     &  17.72, 17.12, 16.53, 15.96, 15.40, 14.85, 14.31, 13.78, 13.26,
     &  12.75,
     &  12.25, 11.75, 11.26, 10.77, 10.29,  9.82,  9.35,  8.88,  8.42,
     &   7.97,
     &   7.52,  7.07,  6.62,  6.18,  5.75,  5.32,  4.89,  4.47,  4.05,
     &   3.64,
     &   3.24,  2.84,  2.45,  2.07,  1.70,  1.34,  0.99,  0.65,  0.31,
     &  0.1/
C
C awsc= sodium chloride
C
      DATA AWSC/10*28.16,
     &  27.17, 26.27, 25.45, 24.69, 23.98, 23.33, 22.72, 22.14, 21.59,
     &  21.08,
     &  20.58, 20.12, 19.67, 19.24, 18.82, 18.43, 18.04, 17.67, 17.32,
     &  16.97,
     &  16.63, 16.31, 15.99, 15.68, 15.38, 15.08, 14.79, 14.51, 14.24,
     &  13.97,
     &  13.70, 13.44, 13.18, 12.93, 12.68, 12.44, 12.20, 11.96, 11.73,
     &  11.50,
     &  11.27, 11.05, 10.82, 10.60, 10.38, 10.16,  9.95,  9.74,  9.52,
     &   9.31,
     &   9.10,  8.89,  8.69,  8.48,  8.27,  8.07,  7.86,  7.65,  7.45,
     &   7.24,
     &   7.04,  6.83,  6.62,  6.42,  6.21,  6.00,  5.79,  5.58,  5.36,
     &   5.15,
     &   4.93,  4.71,  4.48,  4.26,  4.03,  3.80,  3.56,  3.32,  3.07,
     &   2.82,
     &   2.57,  2.30,  2.04,  1.76,  1.48,  1.20,  0.91,  0.61,  0.30,
     &  0.1/
C
C awac= ammonium chloride
C
      DATA AWAC/10*1209.00,
     & 1067.60,949.27,848.62,761.82,686.04,619.16,559.55,505.92,457.25,
     & 412.69,
     & 371.55,333.21,297.13,262.81,229.78,197.59,165.98,135.49,108.57,
     &  88.29,
     &  74.40, 64.75, 57.69, 52.25, 47.90, 44.30, 41.27, 38.65, 36.36,
     &  34.34,
     &  32.52, 30.88, 29.39, 28.02, 26.76, 25.60, 24.51, 23.50, 22.55,
     &  21.65,
     &  20.80, 20.00, 19.24, 18.52, 17.83, 17.17, 16.54, 15.93, 15.35,
     &  14.79,
     &  14.25, 13.73, 13.22, 12.73, 12.26, 11.80, 11.35, 10.92, 10.49,
     &  10.08,
     &   9.67,  9.28,  8.89,  8.51,  8.14,  7.77,  7.42,  7.06,  6.72,
     &   6.37,
     &   6.03,  5.70,  5.37,  5.05,  4.72,  4.40,  4.08,  3.77,  3.45,
     &   3.14,
     &   2.82,  2.51,  2.20,  1.89,  1.57,  1.26,  0.94,  0.62,  0.31,
     &  0.1/
C
C awss= sodium sulfate
C
      DATA AWSS/10*24.10,
     &  23.17, 22.34, 21.58, 20.90, 20.27, 19.69, 19.15, 18.64, 18.17,
     &  17.72,
     &  17.30, 16.90, 16.52, 16.16, 15.81, 15.48, 15.16, 14.85, 14.55,
     &  14.27,
     &  13.99, 13.73, 13.47, 13.21, 12.97, 12.73, 12.50, 12.27, 12.05,
     &  11.84,
     &  11.62, 11.42, 11.21, 11.01, 10.82, 10.63, 10.44, 10.25, 10.07,
     &   9.89,
     &   9.71,  9.53,  9.36,  9.19,  9.02,  8.85,  8.68,  8.51,  8.35,
     &   8.19,
     &   8.02,  7.86,  7.70,  7.54,  7.38,  7.22,  7.06,  6.90,  6.74,
     &   6.58,
     &   6.42,  6.26,  6.10,  5.94,  5.78,  5.61,  5.45,  5.28,  5.11,
     &   4.93,
     &   4.76,  4.58,  4.39,  4.20,  4.01,  3.81,  3.60,  3.39,  3.16,
     &   2.93,
     &   2.68,  2.41,  2.13,  1.83,  1.52,  1.19,  0.86,  0.54,  0.25,
     &  0.1/
C
C awab= ammonium bisulfate
C
      DATA AWAB/10*312.84,
     & 271.43,237.19,208.52,184.28,163.64,145.97,130.79,117.72,106.42,
     &  96.64,
     &  88.16, 80.77, 74.33, 68.67, 63.70, 59.30, 55.39, 51.89, 48.76,
     &  45.93,
     &  43.38, 41.05, 38.92, 36.97, 35.18, 33.52, 31.98, 30.55, 29.22,
     &  27.98,
     &  26.81, 25.71, 24.67, 23.70, 22.77, 21.90, 21.06, 20.27, 19.52,
     &  18.80,
     &  18.11, 17.45, 16.82, 16.21, 15.63, 15.07, 14.53, 14.01, 13.51,
     &  13.02,
     &  12.56, 12.10, 11.66, 11.24, 10.82, 10.42, 10.04,  9.66,  9.29,
     &   8.93,
     &   8.58,  8.24,  7.91,  7.58,  7.26,  6.95,  6.65,  6.35,  6.05,
     &   5.76,
     &   5.48,  5.20,  4.92,  4.64,  4.37,  4.09,  3.82,  3.54,  3.27,
     &   2.99,
     &   2.70,  2.42,  2.12,  1.83,  1.52,  1.22,  0.90,  0.59,  0.28,
     &  0.1/
C
C awsa= sulfuric acid
C
      DATA AWSA/34.00, 33.56, 29.22, 26.55, 24.61, 23.11, 21.89, 20.87,
     &  19.99, 18.45,
     &  17.83, 17.26, 16.73, 16.25, 15.80, 15.38, 14.98, 14.61, 14.26,
     &  13.93,
     &  13.61, 13.30, 13.01, 12.73, 12.47, 12.21, 11.96, 11.72, 11.49,
     &  11.26,
     &  11.04, 10.83, 10.62, 10.42, 10.23, 10.03,  9.85,  9.67,  9.49,
     &   9.31,
     &   9.14,  8.97,  8.81,  8.65,  8.49,  8.33,  8.18,  8.02,  7.87,
     &   7.73,
     &   7.58,  7.44,  7.29,  7.15,  7.01,  6.88,  6.74,  6.61,  6.47,
     &   6.34,
     &   6.21,  6.07,  5.94,  5.81,  5.68,  5.55,  5.43,  5.30,  5.17,
     &   5.04,
     &   4.91,  4.78,  4.65,  4.52,  4.39,  4.26,  4.13,  4.00,  3.86,
     &   3.73,
     &   3.59,  3.45,  3.31,  3.17,  3.02,  2.87,  2.71,  2.56,  2.39,
     &   2.22,
     &   2.05,  1.87,  1.68,  1.48,  1.27,  1.04,  0.80,  0.55,  0.28,
     &  0.1/
C
C awlc= (NH4)3H(SO4)2
C
      DATA AWLC/10*125.37,
     & 110.10, 97.50, 86.98, 78.08, 70.49, 63.97, 58.33, 53.43, 49.14,
     &  45.36,
     &  42.03, 39.07, 36.44, 34.08, 31.97, 30.06, 28.33, 26.76, 25.32,
     &  24.01,
     &  22.81, 21.70, 20.67, 19.71, 18.83, 18.00, 17.23, 16.50, 15.82,
     &  15.18,
     &  14.58, 14.01, 13.46, 12.95, 12.46, 11.99, 11.55, 11.13, 10.72,
     &  10.33,
     &   9.96,  9.60,  9.26,  8.93,  8.61,  8.30,  8.00,  7.72,  7.44,
     &   7.17,
     &   6.91,  6.66,  6.42,  6.19,  5.96,  5.74,  5.52,  5.31,  5.11,
     &   4.91,
     &   4.71,  4.53,  4.34,  4.16,  3.99,  3.81,  3.64,  3.48,  3.31,
     &   3.15,
     &   2.99,  2.84,  2.68,  2.53,  2.37,  2.22,  2.06,  1.91,  1.75,
     &   1.60,
     &   1.44,  1.28,  1.12,  0.95,  0.79,  0.62,  0.45,  0.29,  0.14,
     &  0.1/
C
C awan= ammonium nitrate
C
      DATA AWAN/10*960.19,
     & 853.15,763.85,688.20,623.27,566.92,517.54,473.91,435.06,400.26,
     & 368.89,
     & 340.48,314.63,291.01,269.36,249.46,231.11,214.17,198.50,184.00,
     & 170.58,
     & 158.15,146.66,136.04,126.25,117.24,108.97,101.39, 94.45, 88.11,
     &  82.33,
     &  77.06, 72.25, 67.85, 63.84, 60.16, 56.78, 53.68, 50.81, 48.17,
     &  45.71,
     &  43.43, 41.31, 39.32, 37.46, 35.71, 34.06, 32.50, 31.03, 29.63,
     &  28.30,
     &  27.03, 25.82, 24.67, 23.56, 22.49, 21.47, 20.48, 19.53, 18.61,
     &  17.72,
     &  16.86, 16.02, 15.20, 14.41, 13.64, 12.89, 12.15, 11.43, 10.73,
     &  10.05,
     &   9.38,  8.73,  8.09,  7.47,  6.86,  6.27,  5.70,  5.15,  4.61,
     &   4.09,
     &   3.60,  3.12,  2.66,  2.23,  1.81,  1.41,  1.03,  0.67,  0.32,
     &  0.1/
C
C awsb= sodium bisulfate
C
      DATA AWSB/10*55.99,
     &  53.79, 51.81, 49.99, 48.31, 46.75, 45.28, 43.91, 42.62, 41.39,
     &  40.22,
     &  39.10, 38.02, 36.99, 36.00, 35.04, 34.11, 33.21, 32.34, 31.49,
     &  30.65,
     &  29.84, 29.04, 28.27, 27.50, 26.75, 26.01, 25.29, 24.57, 23.87,
     &  23.17,
     &  22.49, 21.81, 21.15, 20.49, 19.84, 19.21, 18.58, 17.97, 17.37,
     &  16.77,
     &  16.19, 15.63, 15.08, 14.54, 14.01, 13.51, 13.01, 12.53, 12.07,
     &  11.62,
     &  11.19, 10.77, 10.36,  9.97,  9.59,  9.23,  8.87,  8.53,  8.20,
     &   7.88,
     &   7.57,  7.27,  6.97,  6.69,  6.41,  6.14,  5.88,  5.62,  5.36,
     &   5.11,
     &   4.87,  4.63,  4.39,  4.15,  3.92,  3.68,  3.45,  3.21,  2.98,
     &   2.74,
     &   2.49,  2.24,  1.98,  1.72,  1.44,  1.16,  0.87,  0.57,  0.28,
     &  0.1/
C
C *** ZSR RELATIONSHIP PARAMETERS **************************************
c
C awas= ammonium sulfate
C
C      DATA AWAS/33*100.,30,30,30,29.54,28.25,27.06,25.94,
C     & 24.89,23.90,22.97,22.10,21.27,20.48,19.73,19.02,18.34,17.69,
C     & 17.07,16.48,15.91,15.37,14.85,14.34,13.86,13.39,12.94,12.50,
C     & 12.08,11.67,11.27,10.88,10.51,10.14, 9.79, 9.44, 9.10, 8.78,
C     &  8.45, 8.14, 7.83, 7.53, 7.23, 6.94, 6.65, 6.36, 6.08, 5.81,
C     &  5.53, 5.26, 4.99, 4.72, 4.46, 4.19, 3.92, 3.65, 3.38, 3.11,
C     &  2.83, 2.54, 2.25, 1.95, 1.63, 1.31, 0.97, 0.63, 0.30, 0.001/
C
C awsn= sodium nitrate
C
C      DATA AWSN/ 9*1.e5,685.59,
C     & 451.00,336.46,268.48,223.41,191.28,
C     & 167.20,148.46,133.44,121.12,110.83,
C     & 102.09,94.57,88.03,82.29,77.20,72.65,68.56,64.87,61.51,58.44,
C     & 55.62,53.03,50.63,48.40,46.32,44.39,42.57,40.87,39.27,37.76,
C     & 36.33,34.98,33.70,32.48,31.32,30.21,29.16,28.14,27.18,26.25,
C     & 25.35,24.50,23.67,22.87,22.11,21.36,20.65,19.95,19.28,18.62,
C     & 17.99,17.37,16.77,16.18,15.61,15.05,14.51,13.98,13.45,12.94,
C     & 12.44,11.94,11.46,10.98,10.51,10.04, 9.58, 9.12, 8.67, 8.22,
C     &  7.77, 7.32, 6.88, 6.43, 5.98, 5.53, 5.07, 4.61, 4.15, 3.69,
C     &  3.22, 2.76, 2.31, 1.87, 1.47, 1.10, 0.77, 0.48, 0.23, 0.001/
C
C awsc= sodium chloride
C
C      DATA AWSC/
C     &  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
C     &  100., 100., 100., 100., 100., 100., 100., 100., 100.,16.34,
C     & 16.28,16.22,16.15,16.09,16.02,15.95,15.88,15.80,15.72,15.64,
C     & 15.55,15.45,15.36,15.25,15.14,15.02,14.89,14.75,14.60,14.43,
C     & 14.25,14.04,13.81,13.55,13.25,12.92,12.56,12.19,11.82,11.47,
C     & 11.13,10.82,10.53,10.26,10.00, 9.76, 9.53, 9.30, 9.09, 8.88,
C     &  8.67, 8.48, 8.28, 8.09, 7.90, 7.72, 7.54, 7.36, 7.17, 6.99,
C     &  6.81, 6.63, 6.45, 6.27, 6.09, 5.91, 5.72, 5.53, 5.34, 5.14,
C     &  4.94, 4.74, 4.53, 4.31, 4.09, 3.86, 3.62, 3.37, 3.12, 2.85,
C     &  2.58, 2.30, 2.01, 1.72, 1.44, 1.16, 0.89, 0.64, 0.40, 0.18/
C
C awac= ammonium chloride
C
C      DATA AWAC/
C     &  100., 100., 100., 100., 100., 100., 100., 100., 100., 100.,
C     &  100., 100., 100., 100., 100., 100., 100., 100., 100.,31.45,
C     & 31.30,31.14,30.98,30.82,30.65,30.48,30.30,30.11,29.92,29.71,
C     & 29.50,29.29,29.06,28.82,28.57,28.30,28.03,27.78,27.78,27.77,
C     & 27.77,27.43,27.07,26.67,26.21,25.73,25.18,24.56,23.84,23.01,
C     & 22.05,20.97,19.85,18.77,17.78,16.89,16.10,15.39,14.74,14.14,
C     & 13.59,13.06,12.56,12.09,11.65,11.22,10.81,10.42,10.03, 9.66,
C     &  9.30, 8.94, 8.59, 8.25, 7.92, 7.59, 7.27, 6.95, 6.63, 6.32,
C     &  6.01, 5.70, 5.39, 5.08, 4.78, 4.47, 4.17, 3.86, 3.56, 3.25,
C     &  2.94, 2.62, 2.30, 1.98, 1.65, 1.32, 0.97, 0.62, 0.26, 0.13/
C
C awss= sodium sulfate
C
C      DATA AWSS/34*1.e5,23*14.30,14.21,12.53,11.47,
C     & 10.66,10.01, 9.46, 8.99, 8.57, 8.19, 7.85, 7.54, 7.25, 6.98,
C     &  6.74, 6.50, 6.29, 6.08, 5.88, 5.70, 5.52, 5.36, 5.20, 5.04,
C     &  4.90, 4.75, 4.54, 4.34, 4.14, 3.93, 3.71, 3.49, 3.26, 3.02,
C     &  2.76, 2.49, 2.20, 1.89, 1.55, 1.18, 0.82, 0.49, 0.22, 0.001/
C
C awab= ammonium bisulfate
C
C      DATA AWAB/356.45,296.51,253.21,220.47,194.85,
C     & 174.24,157.31,143.16,131.15,120.82,
C     & 111.86,103.99,97.04,90.86,85.31,80.31,75.78,71.66,67.90,64.44,
C     &  61.25,58.31,55.58,53.04,50.68,48.47,46.40,44.46,42.63,40.91,
C     &  39.29,37.75,36.30,34.92,33.61,32.36,31.18,30.04,28.96,27.93,
C     &  26.94,25.99,25.08,24.21,23.37,22.57,21.79,21.05,20.32,19.63,
C     &  18.96,18.31,17.68,17.07,16.49,15.92,15.36,14.83,14.31,13.80,
C     &  13.31,12.83,12.36,11.91,11.46,11.03,10.61,10.20, 9.80, 9.41,
C     &   9.02, 8.64, 8.28, 7.91, 7.56, 7.21, 6.87, 6.54, 6.21, 5.88,
C     &   5.56, 5.25, 4.94, 4.63, 4.33, 4.03, 3.73, 3.44, 3.14, 2.85,
C     &   2.57, 2.28, 1.99, 1.71, 1.42, 1.14, 0.86, 0.57, 0.29, 0.001/
C
C awsa= sulfuric acid
C
C      DATA AWSA/
C     & 34.0,33.56,29.22,26.55,24.61,23.11,21.89,20.87,19.99,
C     & 19.21,18.51,17.87,17.29,16.76,16.26,15.8,15.37,14.95,14.56,
C     & 14.20,13.85,13.53,13.22,12.93,12.66,12.40,12.14,11.90,11.67,
C     & 11.44,11.22,11.01,10.8,10.60,10.4,10.2,10.01,9.83,9.65,9.47,
C     & 9.3,9.13,8.96,8.81,8.64,8.48,8.33,8.17,8.02,7.87,7.72,7.58,
C     & 7.44,7.30,7.16,7.02,6.88,6.75,6.61,6.48,6.35,6.21,6.08,5.95,
C     & 5.82,5.69,5.56,5.44,5.31,5.18,5.05,4.92,4.79,4.66,4.53,4.40,
C     & 4.27,4.14,4.,3.87,3.73,3.6,3.46,3.31,3.17,3.02,2.87,2.72,
C     & 2.56,2.4,2.23,2.05,1.87,1.68,1.48,1.27,1.05,0.807,0.552,0.281/
C
C awlc= (NH4)3H(SO4)2
C
C      DATA AWLC/34*1.e5,17.0,16.5,15.94,15.31,14.71,14.14,
C     & 13.60,13.08,12.59,12.12,11.68,11.25,10.84,10.44,10.07, 9.71,
C     &  9.36, 9.02, 8.70, 8.39, 8.09, 7.80, 7.52, 7.25, 6.99, 6.73,
C     &  6.49, 6.25, 6.02, 5.79, 5.57, 5.36, 5.15, 4.95, 4.76, 4.56,
C     &  4.38, 4.20, 4.02, 3.84, 3.67, 3.51, 3.34, 3.18, 3.02, 2.87,
C     &  2.72, 2.57, 2.42, 2.28, 2.13, 1.99, 1.85, 1.71, 1.57, 1.43,
C     &  1.30, 1.16, 1.02, 0.89, 0.75, 0.61, 0.46, 0.32, 0.16, 0.001/
C
C awan= ammonium nitrate
C
C      DATA AWAN/31*1.e5,
C     &       97.17,92.28,87.66,83.15,78.87,74.84,70.98,67.46,64.11,
C     & 60.98,58.07,55.37,52.85,50.43,48.24,46.19,44.26,42.40,40.70,
C     & 39.10,37.54,36.10,34.69,33.35,32.11,30.89,29.71,28.58,27.46,
C     & 26.42,25.37,24.33,23.89,22.42,21.48,20.56,19.65,18.76,17.91,
C     & 17.05,16.23,15.40,14.61,13.82,13.03,12.30,11.55,10.83,10.14,
C     &  9.44, 8.79, 8.13, 7.51, 6.91, 6.32, 5.75, 5.18, 4.65, 4.14,
C     &  3.65, 3.16, 2.71, 2.26, 1.83, 1.42, 1.03, 0.66, 0.30, 0.001/
C
C awsb= sodium bisulfate
C
C      DATA AWSB/173.72,156.88,142.80,130.85,120.57,
C     & 111.64,103.80,96.88,90.71,85.18,
C     & 80.20,75.69,71.58,67.82,64.37,61.19,58.26,55.53,53.00,50.64,
C     & 48.44,46.37,44.44,42.61,40.90,39.27,37.74,36.29,34.91,33.61,
C     & 32.36,31.18,30.05,28.97,27.94,26.95,26.00,25.10,24.23,23.39,
C     & 22.59,21.81,21.07,20.35,19.65,18.98,18.34,17.71,17.11,16.52,
C     & 15.95,15.40,14.87,14.35,13.85,13.36,12.88,12.42,11.97,11.53,
C     & 11.10,10.69,10.28, 9.88, 9.49, 9.12, 8.75, 8.38, 8.03, 7.68,
C     &  7.34, 7.01, 6.69, 6.37, 6.06, 5.75, 5.45, 5.15, 4.86, 4.58,
C     &  4.30, 4.02, 3.76, 3.49, 3.23, 2.98, 2.73, 2.48, 2.24, 2.01,
C     &  1.78, 1.56, 1.34, 1.13, 0.92, 0.73, 0.53, 0.35, 0.17, 0.001/
C
C *** END OF BLOCK DATA SUBPROGRAM *************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE INIT1
C *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM     
C     SULFATE AEROSOL SYSTEMS (SUBROUTINE ISRP1)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE INIT1 (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION WI, RHI, TEMPI
      DIMENSION WI(NCOMP)
      DOUBLE PRECISION IC,GII,GI0,XX
      DOUBLE PRECISION ln10
      PARAMETER (LN10=2.30258509299404568402D0)
C
C *** SAVE INPUT VARIABLES IN COMMON BLOCK ******************************
C
      IF (IPROB.EQ.0) THEN                 ! FORWARD CALCULATION
         DO I=1,NCOMP
C            W(I) = (WI(I), TINY)
            IF (TINY .GT. (WI(I))) THEN
               W(I) = TINY
            ELSE
               W(I) = WI(I)
            ENDIF
C            WB(I) = 0.d0
         ENDDO
      ELSE
         DO I=1,NCOMP                      ! REVERSE CALCULATION
C            WAER(I) = MAX(WI(I), TINY)
            IF (TINY .GT. (WI(I))) THEN
               WAER(I) = TINY
            ELSE
               WAER(I) = WI(I)
            ENDIF
            W(I)    = ZERO
         ENDDO
      ENDIF
      RH      = RHI
      TEMP    = TEMPI
C
C *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
C
      XK1  = 1.015e-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639    ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805e-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK7  = 1.817     ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK12 = 1.382e2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268    ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XKW  = 1.010e-14 ! H2O              <==> H(aq)     + OH(aq)
C
      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2 = XK21*XK22       
C
C *** CALCULATE DELIQUESCENCE RELATIVE HUMIDITIES (UNICOMPONENT) ********
C
      DRH2SO4  = 0.0000D0
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRLC     = 0.6900D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15d0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF) 
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF) 
         DRLC     = DRLC    *EXP(186.*TCF) 
      ENDIF
C
C *** CALCULATE MUTUAL DELIQUESCENCE RELATIVE HUMIDITIES ****************
C
      DRMLCAB = 0.3780D0              ! (NH4)3H(SO4)2 & NH4HSO4 
      DRMLCAS = 0.6900D0              ! (NH4)3H(SO4)2 & (NH4)2SO4 
CCC      IF (INT(TEMP) .NE. 298) THEN      ! For the time being.
CCC         T0       = 298.15d0
CCC         TCF      = 1.0/TEMP - 1.0/T0
CCC         DRMLCAB  = DRMLCAB*EXP(507.506*TCF) 
CCC         DRMLCAS  = DRMLCAS*EXP(133.865*TCF) 
CCC      ENDIF
C
C *** LIQUID PHASE ******************************************************
C
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY
C
      DO I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1d0
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
      ENDDO
C
      DO I=1,NPAIR
         GAMA(I) = 0.1d0
      ENDDO
C
      DO I=1,NIONS
C         MOLALB(I) = 0.d0
         MOLAL(I)=ZERO
      ENDDO
      COH = ZERO
C
      DO I=1,NGASAQ
         GASAQ(I)=ZERO
      ENDDO
C
C *** SOLID PHASE *******************************************************
C
      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO
C
C *** GAS PHASE *********************************************************
C
      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO
C
C *** CALCULATE ZSR PARAMETERS ******************************************
C
      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  ! Position in ZSR arrays
      IRH    = MAX (IRH, 1)
C
C      M0(01) = AWSC(IRH)      ! NACl
C      IF (M0(01) .LT. 100.0) THEN
C         IC = M0(01)
C         CALL KMTAB(IC,298.0,     GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(02) = AWSS(IRH)      ! (NA)2SO4
C      IF (M0(02) .LT. 100.0) THEN
C         IC = 3.0*M0(02)
C         CALL KMTAB(IC,298.0,     XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(03) = AWSN(IRH)      ! NANO3
C      IF (M0(03) .LT. 100.0) THEN
C         IC = M0(03)
C         CALL KMTAB(IC,298.0,     XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(04) = AWAS(IRH)      ! (NH4)2SO4
C      IF (M0(04) .LT. 100.0) THEN
C         IC = 3.0*M0(04)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(05) = AWAN(IRH)      ! NH4NO3
C      IF (M0(05) .LT. 100.0) THEN
C         IC     = M0(05)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
C         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(06) = AWAC(IRH)      ! NH4CL
C      IF (M0(06) .LT. 100.0) THEN
C         IC = M0(06)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
C         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(07) = AWSA(IRH)      ! 2H-SO4
C      IF (M0(07) .LT. 100.0) THEN
C         IC = 3.0*M0(07)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
C         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(08) = AWSA(IRH)      ! H-HSO4
CCC      IF (M0(08) .LT. 100.0) THEN     ! These are redundant, because M0(8) is not used
CCC         IC = M0(08)
CCC         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCCCCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
CCC         M0(08) = M0(08)*EXP(LN10*(GI0-GII))
CCC      ENDIF
C
      M0(09) = AWAB(IRH)      ! NH4HSO4
C      IF (M0(09) .LT. 100.0) THEN
C         IC = M0(09)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
C         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(12) = AWSB(IRH)      ! NAHSO4
C      IF (M0(12) .LT. 100.0) THEN
C         IC = M0(12)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
C         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(13) = AWLC(IRH)      ! (NH4)3H(SO4)2
C      IF (M0(13) .LT. 100.0) THEN
C         IC     = 4.0*M0(13)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
C         G130   = 0.2*(3.0*GI0+2.0*GII)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
C         G13I   = 0.2*(3.0*GI0+2.0*GII)
C         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
C      ENDIF
C
C *** OTHER INITIALIZATIONS *********************************************
C
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = 'XX'
      SULRATW = 2.D0
      SODRAT  = ZERO
      NOFER   = 0
      STKOFL  =.FALSE.
      DO I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
      ENDDO
C
C *** END OF SUBROUTINE INIT1 *******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE INIT2
C *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM,
C     NITRATE, SULFATE AEROSOL SYSTEMS (SUBROUTINE ISRP2)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE INIT2 (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION WI, RHI, TEMPI
      DIMENSION WI(NCOMP)
      LOGICAL FLAGNP
      DOUBLE PRECISION IC,GII,GI0,XX
      DOUBLE PRECISION ln10
      PARAMETER (LN10=2.30258509299404568402D0)
C
C *** SAVE INPUT VARIABLES IN COMMON BLOCK ******************************
C
      IF (IPROB.EQ.0) THEN                 ! FORWARD CALCULATION
         DO I=1,NCOMP
C            W(I) = (WI(I), TINY)
            IF (TINY .GT. (WI(I))) THEN
               W(I) = ZERO! TINY
            ELSE
               W(I) = WI(I)
            ENDIF
C            WB(I) = ZERO
         ENDDO
      ELSE
         DO I=1,NCOMP                   ! REVERSE CALCULATION
C            WAER(I) = MAX(WI(I), TINY)
            IF (TINY .GT. (WI(I))) THEN
               WAER(I) = TINY
            ELSE
               WAER(I) = WI(I)
            ENDIF
            W(I)    = ZERO
         ENDDO
      ENDIF
      RH      = RHI
      TEMP    = TEMPI
C
C *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
C
      XK1  = 1.015e-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639    ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805e-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK4  = 2.511e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
CCC      XK4  = 3.638e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! SEQUIL
      XK41 = 2.100e5   ! HNO3(g)          <==> HNO3(aq)
      XK7  = 1.817     ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK10 = 5.746e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR
CCC      XK10 = 2.985e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! SEQUIL
      XK12 = 1.382e2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268    ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XKW  = 1.010e-14 ! H2O              <==> H(aq)     + OH(aq)
C
      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) !ISORR
CCC         XK4 = XK4 *EXP( 29.47*(T0T-1.0) + 16.840*COEF) ! SEQUIL
         XK41= XK41*EXP( 29.17*(T0T-1.0) + 16.830*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) ! ISORR
CCC         XK10= XK10*EXP(-75.11*(T0T-1.0) + 13.460*COEF) ! SEQUIL
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2  = XK21*XK22       
      XK42 = XK4/XK41
C
C *** CALCULATE DELIQUESCENCE RELATIVE HUMIDITIES (UNICOMPONENT) ********
C
      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRNH4NO3 = 0.6183D0
      DRLC     = 0.6900D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNH4NO3 = DRNH4NO3*EXP(852.*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF) 
         DRLC     = DRLC    *EXP(186.*TCF) 
         DRNH4NO3 = MIN ((DRNH4NO3),(DRNH42S4)) ! ADJUST FOR DRH CROSSOVER AT T<271K
      ENDIF
C
C *** CALCULATE MUTUAL DELIQUESCENCE RELATIVE HUMIDITIES ****************
C
      DRMLCAB = 0.3780D0              ! (NH4)3H(SO4)2 & NH4HSO4 
      DRMLCAS = 0.6900D0              ! (NH4)3H(SO4)2 & (NH4)2SO4 
      DRMASAN = 0.6000D0              ! (NH4)2SO4     & NH4NO3
CCC      IF (INT(TEMP) .NE. 298) THEN    ! For the time being
CCC         T0       = 298.15d0
CCC         TCF      = 1.0/TEMP - 1.0/T0
CCC         DRMLCAB  = DRMLCAB*EXP( 507.506*TCF) 
CCC         DRMLCAS  = DRMLCAS*EXP( 133.865*TCF) 
CCC         DRMASAN  = DRMASAN*EXP(1269.068*TCF)
CCC      ENDIF
C
C *** LIQUID PHASE ******************************************************
C
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY
C
      DO I=1,NPAIR
         MOLALR(I)=ZERO
         GAMA(I)  =0.1D0
         GAMIN(I) = 1.D10 ! GREAT
         GAMOU(I) = 1.D10 !GREAT
         M0(I)    =1d5
      ENDDO
C
      DO I=1,NPAIR
         GAMA(I) = 0.1d0
      ENDDO
C
      DO I=1,NIONS
         MOLAL(I)=ZERO
      ENDDO
      COH = ZERO
C
      DO I=1,NGASAQ
         GASAQ(I)=ZERO
      ENDDO
C
C *** SOLID PHASE *******************************************************
C
      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO
C
C *** GAS PHASE *********************************************************
C
      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO
C
C *** CALCULATE ZSR PARAMETERS ******************************************
C
      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  ! Position in ZSR arrays
      IRH    = MAX (IRH, 1)
C
C      M0(01) = AWSC(IRH)      ! NACl
C      IF (M0(01) .LT. 100.0) THEN
C         IC = M0(01)
C         CALL KMTAB(IC,298.0,     GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(02) = AWSS(IRH)      ! (NA)2SO4
C      IF (M0(02) .LT. 100.0) THEN
C         IC = 3.0*M0(02)
C         CALL KMTAB(IC,298.0,     XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(03) = AWSN(IRH)      ! NANO3
C      IF (M0(03) .LT. 100.0) THEN
C         IC = M0(03)
C         CALL KMTAB(IC,298.0,     XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(04) = AWAS(IRH)      ! (NH4)2SO4
C      IF (M0(04) .LT. 100.0) THEN
C         IC = 3.0*M0(04)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(05) = AWAN(IRH)      ! NH4NO3
C      IF (M0(05) .LT. 100.0) THEN
C         IC     = M0(05)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
C         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(06) = AWAC(IRH)      ! NH4CL
C      IF (M0(06) .LT. 100.0) THEN
C         IC = M0(06)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
C         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(07) = AWSA(IRH)      ! 2H-SO4
C      IF (M0(07) .LT. 100.0) THEN
C         IC = 3.0*M0(07)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
C         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(08) = AWSA(IRH)      ! H-HSO4
CCC      IF (M0(08) .LT. 100.0) THEN     ! These are redundant, because M0(8) is not used
CCC         IC = M0(08)
CCC         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCCCCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
CCC         M0(08) = M0(08)*EXP(LN10*(GI0-GII))
CCC      ENDIF
C
      M0(09) = AWAB(IRH)      ! NH4HSO4
C      IF (M0(09) .LT. 100.0) THEN
C         IC = M0(09)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
C         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
C      ENDIF
C
C      M0(12) = AWSB(IRH)      ! NAHSO4
C      IF (M0(12) .LT. 100.0) THEN
C         IC = M0(12)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
C         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(13) = AWLC(IRH)      ! (NH4)3H(SO4)2
C      IF (M0(13) .LT. 100.0) THEN
C         IC     = 4.0*M0(13)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
C         G130   = 0.2*(3.0*GI0+2.0*GII)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
C         G13I   = 0.2*(3.0*GI0+2.0*GII)
C         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
C      ENDIF
C
C *** OTHER INITIALIZATIONS *********************************************
C
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      FLAGNP  = .FALSE.
      NONPHYS = .FALSE.
      SCASE   = 'XX'
      SULRATW = 2.D0
      SODRAT  = ZERO
      NOFER   = 0
      STKOFL  =.FALSE.
      DO 60 I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
   60 CONTINUE
C
C *** END OF SUBROUTINE INIT2 *******************************************
C
      END





C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISOINIT3
C *** THIS SUBROUTINE INITIALIZES ALL GLOBAL VARIABLES FOR AMMONIUM,
C     SODIUM, CHLORIDE, NITRATE, SULFATE AEROSOL SYSTEMS (SUBROUTINE 
C     ISRP3)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISOINIT3 (WI, RHI, TEMPI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: WI(NCOMP), RHI, TEMPI
      DOUBLE PRECISION      IC,GII,GI0,XX
      DOUBLE PRECISION ln10
      
      PARAMETER (LN10=2.30258509299404568402D0)
C
C *** SAVE INPUT VARIABLES IN COMMON BLOCK ******************************
C
      IF (IPROB.EQ.0) THEN                 ! FORWARD CALCULATION
         DO I=1,NCOMP
C            W(I) = (WI(I), TINY)
            IF (TINY .GT. (WI(I))) THEN
               W(I) = TINY
            ELSE
               W(I) = WI(I)
            ENDIF
         ENDDO
      ELSE
         DO I=1,NCOMP                      ! REVERSE CALCULATION
C            WAER(I) = MAX(WI(I), TINY)
            IF (TINY .GT. (WI(I))) THEN
               WAER(I) = TINY
            ELSE
               WAER(I) = WI(I)
            ENDIF
            W(I)    = ZERO
         ENDDO
      ENDIF
      RH      = RHI
      TEMP    = TEMPI
C
C *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
C
      XK1  = 1.015D-2  ! HSO4(aq)         <==> H(aq)     + SO4(aq)
      XK21 = 57.639D0  ! NH3(g)           <==> NH3(aq)
      XK22 = 1.805D-5  ! NH3(aq)          <==> NH4(aq)   + OH(aq)
      XK3  = 1.971D6   ! HCL(g)           <==> H(aq)     + CL(aq)
      XK31 = 2.500e3   ! HCL(g)           <==> HCL(aq)
      XK4  = 2.511e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! ISORR
CCC      XK4  = 3.638e6   ! HNO3(g)          <==> H(aq)     + NO3(aq) ! SEQUIL
      XK41 = 2.100e5   ! HNO3(g)          <==> HNO3(aq)
      XK5  = 0.4799D0  ! NA2SO4(s)        <==> 2*NA(aq)  + SO4(aq)
      XK6  = 1.086D-16 ! NH4CL(s)         <==> NH3(g)    + HCL(g)
      XK7  = 1.817D0   ! (NH4)2SO4(s)     <==> 2*NH4(aq) + SO4(aq)
      XK8  = 37.661D0  ! NACL(s)          <==> NA(aq)    + CL(aq)
      XK10 = 5.746D-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! ISORR
CCC      XK10 = 2.985e-17 ! NH4NO3(s)        <==> NH3(g)    + HNO3(g) ! SEQUIL
      XK11 = 2.413D4   ! NAHSO4(s)        <==> NA(aq)    + HSO4(aq)
      XK12 = 1.382D2   ! NH4HSO4(s)       <==> NH4(aq)   + HSO4(aq)
      XK13 = 29.268D0  ! (NH4)3H(SO4)2(s) <==> 3*NH4(aq) + HSO4(aq) + SO4(aq)
      XK14 = 22.05D0   ! NH4CL(s)         <==> NH4(aq)   + CL(aq)
      XKW  = 1.010D-14 ! H2O              <==> H(aq)     + OH(aq)
      XK9  = 11.977D0  ! NANO3(s)         <==> NA(aq)    + NO3(aq)
C
      IF (INT(TEMP) .NE. 298) THEN   ! FOR T != 298K or 298.15K
         T0  = 298.15D0
         T0T = T0/TEMP
         COEF= 1.0+LOG(T0T)-T0T
         XK1 = XK1 *EXP(  8.85*(T0T-1.0) + 25.140*COEF)
         XK21= XK21*EXP( 13.79*(T0T-1.0) -  5.393*COEF)
         XK22= XK22*EXP( -1.50*(T0T-1.0) + 26.920*COEF)
         XK3 = XK3 *EXP( 30.20*(T0T-1.0) + 19.910*COEF)
         XK31= XK31*EXP( 30.20*(T0T-1.0) + 19.910*COEF)
         XK4 = XK4 *EXP( 29.17*(T0T-1.0) + 16.830*COEF) !ISORR
CCC         XK4 = XK4 *EXP( 29.47*(T0T-1.0) + 16.840*COEF) ! SEQUIL
         XK41= XK41*EXP( 29.17*(T0T-1.0) + 16.830*COEF)
         XK5 = XK5 *EXP(  0.98*(T0T-1.0) + 39.500*COEF)
         XK6 = XK6 *EXP(-71.00*(T0T-1.0) +  2.400*COEF)
         XK7 = XK7 *EXP( -2.65*(T0T-1.0) + 38.570*COEF)
         XK8 = XK8 *EXP( -1.56*(T0T-1.0) + 16.900*COEF)
         XK9 = XK9 *EXP( -8.22*(T0T-1.0) + 16.010*COEF)
         XK10= XK10*EXP(-74.38*(T0T-1.0) +  6.120*COEF) ! ISORR
CCC         XK10= XK10*EXP(-75.11*(T0T-1.0) + 13.460*COEF) ! SEQUIL
         XK11= XK11*EXP(  0.79*(T0T-1.0) + 14.746*COEF)
         XK12= XK12*EXP( -2.87*(T0T-1.0) + 15.830*COEF)
         XK13= XK13*EXP( -5.19*(T0T-1.0) + 54.400*COEF)
         XK14= XK14*EXP( 24.55*(T0T-1.0) + 16.900*COEF)
         XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)
      ENDIF
      XK2  = XK21*XK22       
      XK42 = XK4/XK41
      XK32 = XK3/XK31
C
C *** CALCULATE DELIQUESCENCE RELATIVE HUMIDITIES (UNICOMPONENT) ********
C
      DRH2SO4  = ZERO
      DRNH42S4 = 0.7997D0
      DRNH4HS4 = 0.4000D0
      DRLC     = 0.6900D0
      DRNACL   = 0.7528D0
      DRNANO3  = 0.7379D0
      DRNH4CL  = 0.7710D0
      DRNH4NO3 = 0.6183D0
      DRNA2SO4 = 0.9300D0
      DRNAHSO4 = 0.5200D0
      IF (INT(TEMP) .NE. 298) THEN
         T0       = 298.15D0
         TCF      = 1.0/TEMP - 1.0/T0
         DRNACL   = DRNACL  *EXP( 25.*TCF)
         DRNANO3  = DRNANO3 *EXP(304.*TCF)
         DRNA2SO4 = DRNA2SO4*EXP( 80.*TCF)
         DRNH4NO3 = DRNH4NO3*EXP(852.*TCF)
         DRNH42S4 = DRNH42S4*EXP( 80.*TCF)
         DRNH4HS4 = DRNH4HS4*EXP(384.*TCF) 
         DRLC     = DRLC    *EXP(186.*TCF)
         DRNH4CL  = DRNH4Cl *EXP(239.*TCF)
         DRNAHSO4 = DRNAHSO4*EXP(-45.*TCF) 
C
C *** ADJUST FOR DRH "CROSSOVER" AT LOW TEMPERATURES
C
C         DRNH4NO3  = MIN (DRNH4NO3, DRNH4CL, DRNH42S4, DRNANO3, DRNACL) slc.1.2011 due to TAPENADE's FORTRAN parser
         DRNH42S4  = MIN (DRNH42S4, DRNANO3, DRNACL)           ! slc.1.2011 due to TAPENADE's FORTRAN parser
         DRNH4NO3  = MIN (DRNH4NO3, DRNH4CL, DRNH42S4)         ! slc.1.2011 due to TAPENADE's FORTRAN parser
         DRNANO3   = MIN (DRNANO3, DRNACL)
         DRNH4CL   = MIN (DRNH4Cl, DRNH42S4)
C
      ENDIF
C
C *** CALCULATE MUTUAL DELIQUESCENCE RELATIVE HUMIDITIES ****************
C
      DRMLCAB = 0.378D0    ! (NH4)3H(SO4)2 & NH4HSO4 
      DRMLCAS = 0.690D0    ! (NH4)3H(SO4)2 & (NH4)2SO4 
      DRMASAN = 0.600D0    ! (NH4)2SO4     & NH4NO3
      DRMG1   = 0.460D0    ! (NH4)2SO4, NH4NO3, NA2SO4, NH4CL
      DRMG2   = 0.691D0    ! (NH4)2SO4, NA2SO4, NH4CL
      DRMG3   = 0.697D0    ! (NH4)2SO4, NA2SO4
      DRMH1   = 0.240D0    ! NA2SO4, NANO3, NACL, NH4NO3, NH4CL
      DRMH2   = 0.596D0    ! NA2SO4, NANO3, NACL, NH4CL
      DRMI1   = 0.240D0    ! LC, NAHSO4, NH4HSO4, NA2SO4, (NH4)2SO4
      DRMI2   = 0.363D0    ! LC, NAHSO4, NA2SO4, (NH4)2SO4  - NO DATA -
      DRMI3   = 0.610D0    ! LC, NA2SO4, (NH4)2SO4 
      DRMQ1   = 0.494D0    ! (NH4)2SO4, NH4NO3, NA2SO4
      DRMR1   = 0.663D0    ! NA2SO4, NANO3, NACL
      DRMR2   = 0.735D0    ! NA2SO4, NACL
      DRMR3   = 0.673D0    ! NANO3, NACL
      DRMR4   = 0.694D0    ! NA2SO4, NACL, NH4CL
      DRMR5   = 0.731D0    ! NA2SO4, NH4CL
      DRMR6   = 0.596D0    ! NA2SO4, NANO3, NH4CL
      DRMR7   = 0.380D0    ! NA2SO4, NANO3, NACL, NH4NO3
      DRMR8   = 0.380D0    ! NA2SO4, NACL, NH4NO3
      DRMR9   = 0.494D0    ! NA2SO4, NH4NO3
      DRMR10  = 0.476D0    ! NA2SO4, NANO3, NH4NO3
      DRMR11  = 0.340D0    ! NA2SO4, NACL, NH4NO3, NH4CL
      DRMR12  = 0.460D0    ! NA2SO4, NH4NO3, NH4CL
      DRMR13  = 0.438D0    ! NA2SO4, NANO3, NH4NO3, NH4CL
CCC      IF (INT(TEMP) .NE. 298) THEN
CCC         T0       = 298.15d0
CCC         TCF      = 1.0/TEMP - 1.0/T0
CCC         DRMLCAB  = DRMLCAB*EXP( 507.506*TCF) 
CCC         DRMLCAS  = DRMLCAS*EXP( 133.865*TCF) 
CCC         DRMASAN  = DRMASAN*EXP(1269.068*TCF)
CCC         DRMG1    = DRMG1  *EXP( 572.207*TCF)
CCC         DRMG2    = DRMG2  *EXP(  58.166*TCF)
CCC         DRMG3    = DRMG3  *EXP(  22.253*TCF)
CCC         DRMH1    = DRMH1  *EXP(2116.542*TCF)
CCC         DRMH2    = DRMH2  *EXP( 650.549*TCF)
CCC         DRMI1    = DRMI1  *EXP( 565.743*TCF)
CCC         DRMI2    = DRMI2  *EXP(  91.745*TCF)
CCC         DRMI3    = DRMI3  *EXP( 161.272*TCF)
CCC         DRMQ1    = DRMQ1  *EXP(1616.621*TCF)
CCC         DRMR1    = DRMR1  *EXP( 292.564*TCF)
CCC         DRMR2    = DRMR2  *EXP(  14.587*TCF)
CCC         DRMR3    = DRMR3  *EXP( 307.907*TCF)
CCC         DRMR4    = DRMR4  *EXP(  97.605*TCF)
CCC         DRMR5    = DRMR5  *EXP(  98.523*TCF)
CCC         DRMR6    = DRMR6  *EXP( 465.500*TCF)
CCC         DRMR7    = DRMR7  *EXP( 324.425*TCF)
CCC         DRMR8    = DRMR8  *EXP(2660.184*TCF)
CCC         DRMR9    = DRMR9  *EXP(1617.178*TCF)
CCC         DRMR10   = DRMR10 *EXP(1745.226*TCF)
CCC         DRMR11   = DRMR11 *EXP(3691.328*TCF)
CCC         DRMR12   = DRMR12 *EXP(1836.842*TCF)
CCC         DRMR13   = DRMR13 *EXP(1967.938*TCF)
CCC      ENDIF
C
C *** LIQUID PHASE ******************************************************
C
      CHNO3  = ZERO
      CHCL   = ZERO
      CH2SO4 = ZERO
      COH    = ZERO
      WATER  = TINY
C
      DO I=1,NPAIR
         MOLALR(I)=ZERO
C         MOLALRB(I) = ZERO
         GAMA(I)  =0.1
         GAMIN(I) =GREAT
         GAMOU(I) =GREAT
         M0(I)    =1d5
      ENDDO
C
      DO I=1,NPAIR
         GAMA(I) = 0.1d0
      ENDDO
C
      DO I=1,NIONS
         MOLAL(I)=ZERO
      ENDDO
      COH = ZERO
C
      DO I=1,NGASAQ
         GASAQ(I)=ZERO
      ENDDO
C
C *** SOLID PHASE *******************************************************
C
      CNH42S4= ZERO
      CNH4HS4= ZERO
      CNACL  = ZERO
      CNA2SO4= ZERO
      CNANO3 = ZERO
      CNH4NO3= ZERO
      CNH4CL = ZERO
      CNAHSO4= ZERO
      CLC    = ZERO
C
C *** GAS PHASE *********************************************************
C
      GNH3   = ZERO
      GHNO3  = ZERO
      GHCL   = ZERO
C
C *** CALCULATE ZSR PARAMETERS ******************************************
C
      IRH    = MIN (INT(RH*NZSR+0.5),NZSR)  ! Position in ZSR arrays
      IRH    = MAX (IRH, 1)
C
      M0(01) = AWSC(IRH)      ! NACl
C      IF (M0(01) .LT. 100.0) THEN
C         IC = M0(01)
C         CALL KMTAB(IC,298.0,     GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(01) = M0(01)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(02) = AWSS(IRH)      ! (NA)2SO4
C      IF (M0(02) .LT. 100.0) THEN
C         IC = 3.0*M0(02)
C         CALL KMTAB(IC,298.0,     XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(02) = M0(02)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(03) = AWSN(IRH)      ! NANO3
C      IF (M0(03) .LT. 100.0) THEN
C         IC = M0(03)
C         CALL KMTAB(IC,298.0,     XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(03) = M0(03)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(04) = AWAS(IRH)      ! (NH4)2SO4
C      IF (M0(04) .LT. 100.0) THEN
C         IC = 3.0*M0(04)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX,XX)
C         M0(04) = M0(04)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(05) = AWAN(IRH)      ! NH4NO3
C      IF (M0(05) .LT. 100.0) THEN
C        IC     = M0(05)
C        CALL KMTAB(IC,298.0,     XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX,XX)
C        CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX,XX)
C         M0(05) = M0(05)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(06) = AWAC(IRH)      ! NH4CL
C      IF (M0(06) .LT. 100.0) THEN
C         IC = M0(06)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX,XX)
C         M0(06) = M0(06)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(07) = AWSA(IRH)      ! 2H-SO4
C      IF (M0(07) .LT. 100.0) THEN
C         IC = 3.0*M0(07)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX,XX)
C         M0(07) = M0(07)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(08) = AWSA(IRH)      ! H-HSO4
CCC      IF (M0(08) .LT. 100.0) THEN     ! These are redundant, because M0(8) is not used
CCC         IC = M0(08)
CCC         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX,XX)
CCCCCC         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX,XX)
CCC         M0(08) = M0(08)*EXP(LN10*(GI0-GII))
CCC      ENDIF
C
      M0(09) = AWAB(IRH)      ! NH4HSO4
C      IF (M0(09) .LT. 100.0) THEN
C         IC = M0(09)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,GI0,XX,XX,XX)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,GII,XX,XX,XX)
C         M0(09) = M0(09)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(12) = AWSB(IRH)      ! NAHSO4
C      IF (M0(12) .LT. 100.0) THEN
C         IC = M0(12)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GI0)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,XX,GII)
C         M0(12) = M0(12)*EXP(LN10*(GI0-GII))
C      ENDIF
C
      M0(13) = AWLC(IRH)      ! (NH4)3H(SO4)2
C      IF (M0(13) .LT. 100.0) THEN
C         IC     = 4.0*M0(13)
C         CALL KMTAB(IC,298.0,     XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
C         G130   = 0.2*(3.0*GI0+2.0*GII)
C         CALL KMTAB(IC,SNGL(TEMP),XX,XX,XX,GI0,XX,XX,XX,XX,GII,XX,XX,XX)
C         G13I   = 0.2*(3.0*GI0+2.0*GII)
C         M0(13) = M0(13)*EXP(LN10*SNGL(G130-G13I))
C      ENDIF
C
C *** OTHER INITIALIZATIONS *********************************************
C
      ICLACT  = 0
      CALAOU  = .TRUE.
      CALAIN  = .TRUE.
      FRST    = .TRUE.
      SCASE   = 'XX'
      SULRATW = 2.D0
      NOFER   = 0
      STKOFL  =.FALSE.
      DO I=1,NERRMX
         ERRSTK(I) =-999
         ERRMSG(I) = 'MESSAGE N/A'
      ENDDO
C
C *** END OF SUBROUTINE ISOINIT3 *******************************************
C
      END
      
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ADJUST
C *** ADJUSTS FOR MASS BALANCE BETWEEN VOLATILE SPECIES AND SULFATE
C     FIRST CALCULATE THE EXCESS OF EACH PRECURSOR, AND IF IT EXISTS, THEN
C     ADJUST SEQUENTIALY AEROSOL PHASE SPECIES WHICH CONTAIN THE EXCESS
C     PRECURSOR.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ADJUST (WI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION WI(*)
C
C *** FOR AMMONIUM *****************************************************
C
      IF (IPROB.EQ.0) THEN         ! Calculate excess (solution - input)
         EXNH4 = GNH3 + MOLAL(3) + CNH4CL + CNH4NO3 + CNH4HS4
     &                + 2D0*CNH42S4       + 3D0*CLC
     &          -WI(3)
      ELSE
         EXNH4 = MOLAL(3) + CNH4CL + CNH4NO3 + CNH4HS4 + 2D0*CNH42S4
     &                    + 3D0*CLC
     &          -WI(3)

      ENDIF
      EXNH4 = MAX(EXNH4,ZERO)
      IF ((EXNH4).LT.TINY) GOTO 20    ! No excess NH4, go to next precursor
C
      IF ((MOLAL(3)).GT.(EXNH4)) THEN   ! Adjust aqueous phase NH4
         MOLAL(3) = MOLAL(3) - EXNH4
         GOTO 20
      ELSE
         EXNH4    = EXNH4 - MOLAL(3)
         MOLAL(3) = ZERO
      ENDIF
C
      IF ((CNH4CL).GT.(EXNH4)) THEN     ! Adjust NH4Cl(s)
         CNH4CL   = CNH4CL - EXNH4  ! more solid than excess
         GHCL     = GHCL   + EXNH4  ! evaporate Cl to gas phase
         GOTO 20
      ELSE                          ! less solid than excess
         GHCL     = GHCL   + CNH4CL ! evaporate into gas phase
         EXNH4    = EXNH4  - CNH4CL ! reduce excess
         CNH4CL   = ZERO            ! zero salt concentration
      ENDIF
C
      IF ((CNH4NO3).GT.(EXNH4)) THEN    ! Adjust NH4NO3(s)
         CNH4NO3  = CNH4NO3- EXNH4  ! more solid than excess
         GHNO3    = GHNO3  + EXNH4  ! evaporate NO3 to gas phase
         GOTO 20
      ELSE                          ! less solid than excess
         GHNO3    = GHNO3  + CNH4NO3! evaporate into gas phase
         EXNH4    = EXNH4  - CNH4NO3! reduce excess
         CNH4NO3  = ZERO            ! zero salt concentration
      ENDIF
C
      IF ((CLC).GT.3d0*(EXNH4)) THEN    ! Adjust (NH4)3H(SO4)2(s)
         CLC      = CLC - EXNH4/3d0 ! more solid than excess
         GOTO 20
      ELSE                          ! less solid than excess
         EXNH4    = EXNH4 - 3d0*CLC ! reduce excess
         CLC      = ZERO            ! zero salt concentration
      ENDIF
C
      IF ((CNH4HS4).GT.(EXNH4)) THEN    ! Adjust NH4HSO4(s)
         CNH4HS4  = CNH4HS4- EXNH4  ! more solid than excess
         GOTO 20
      ELSE                          ! less solid than excess
         EXNH4    = EXNH4  - CNH4HS4! reduce excess
         CNH4HS4  = ZERO            ! zero salt concentration
      ENDIF
C
      IF ((CNH42S4).GT.(EXNH4)) THEN    ! Adjust (NH4)2SO4(s)
         CNH42S4  = CNH42S4- EXNH4  ! more solid than excess
         GOTO 20
      ELSE                          ! less solid than excess
         EXNH4    = EXNH4  - CNH42S4! reduce excess
         CNH42S4  = ZERO            ! zero salt concentration
      ENDIF
C
C *** FOR NITRATE ******************************************************
C
 20   IF (IPROB.EQ.0) THEN         ! Calculate excess (solution - input)
         EXNO3 = GHNO3 + MOLAL(7) + CNH4NO3
     &          -WI(4)
      ELSE
         EXNO3 = MOLAL(7) + CNH4NO3
     &          -WI(4)
      ENDIF
      EXNO3 = MAX(EXNO3,ZERO)
      IF ((EXNO3).LT.TINY) GOTO 30    ! No excess NO3, go to next precursor
C
      IF ((MOLAL(7)).GT.(EXNO3)) THEN   ! Adjust aqueous phase NO3
         MOLAL(7) = MOLAL(7) - EXNO3
         GOTO 30
      ELSE
         EXNO3    = EXNO3 - MOLAL(7)
         MOLAL(7) = ZERO
      ENDIF
C
      IF ((CNH4NO3).GT.(EXNO3)) THEN    ! Adjust NH4NO3(s)
         CNH4NO3  = CNH4NO3- EXNO3  ! more solid than excess
         GNH3     = GNH3   + EXNO3  ! evaporate NO3 to gas phase
         GOTO 30
      ELSE                          ! less solid than excess
         GNH3     = GNH3   + CNH4NO3! evaporate into gas phase
         EXNO3    = EXNO3  - CNH4NO3! reduce excess
         CNH4NO3  = ZERO            ! zero salt concentration
      ENDIF
C
C *** FOR CHLORIDE *****************************************************
C
 30   IF (IPROB.EQ.0) THEN         ! Calculate excess (solution - input)
         EXCl = GHCL + MOLAL(4) + CNH4CL
     &         -WI(5)
      ELSE
         EXCl = MOLAL(4) + CNH4CL
     &         -WI(5)
      ENDIF
      EXCl = MAX(EXCl,ZERO)
      IF ((EXCl).LT.TINY) GOTO 40    ! No excess Cl, go to next precursor
C
      IF ((MOLAL(4)).GT.(EXCL)) THEN   ! Adjust aqueous phase Cl
         MOLAL(4) = MOLAL(4) - EXCL
         GOTO 40
      ELSE
         EXCL     = EXCL - MOLAL(4)
         MOLAL(4) = ZERO
      ENDIF
C
      IF ((CNH4CL).GT.(EXCL)) THEN      ! Adjust NH4Cl(s)
         CNH4CL   = CNH4CL - EXCL   ! more solid than excess
         GHCL     = GHCL   + EXCL   ! evaporate Cl to gas phase
         GOTO 40
      ELSE                          ! less solid than excess
         GHCL     = GHCL   + CNH4CL ! evaporate into gas phase
         EXCL     = EXCL   - CNH4CL ! reduce excess
         CNH4CL   = ZERO            ! zero salt concentration
      ENDIF
C
C *** FOR SULFATE ******************************************************
C
 40   EXS4 = MOLAL(5) + MOLAL(6) + 2.d0*CLC + CNH42S4 + CNH4HS4 +
     &       CNA2SO4  + CNAHSO4 - WI(2)
      EXS4 = MAX(EXS4,ZERO)        ! Calculate excess (solution - input)
      IF ((EXS4).LT.TINY) GOTO 50    ! No excess SO4, return
C
      IF ((MOLAL(6)).GT.(EXS4)) THEN   ! Adjust aqueous phase HSO4
         MOLAL(6) = MOLAL(6) - EXS4
         GOTO 50
      ELSE
         EXS4     = EXS4 - MOLAL(6)
         MOLAL(6) = ZERO
      ENDIF
C
      IF ((MOLAL(5)).GT.(EXS4)) THEN   ! Adjust aqueous phase SO4
         MOLAL(5) = MOLAL(5) - EXS4
         GOTO 50
      ELSE
         EXS4     = EXS4 - MOLAL(5)
         MOLAL(5) = ZERO
      ENDIF
C
      IF ((CLC).GT.2d0*(EXS4)) THEN     ! Adjust (NH4)3H(SO4)2(s)
         CLC      = CLC - EXS4/2d0  ! more solid than excess
         GNH3     = GNH3 +1.5d0*EXS4! evaporate NH3 to gas phase
         GOTO 50
      ELSE                          ! less solid than excess
         GNH3     = GNH3 + 1.5d0*CLC! evaporate NH3 to gas phase
         EXS4     = EXS4 - 2d0*CLC  ! reduce excess
         CLC      = ZERO            ! zero salt concentration
      ENDIF
C
      IF ((CNH4HS4).GT.(EXS4)) THEN     ! Adjust NH4HSO4(s)
         CNH4HS4  = CNH4HS4 - EXS4  ! more solid than excess
         GNH3     = GNH3 + EXS4     ! evaporate NH3 to gas phase
         GOTO 50
      ELSE                          ! less solid than excess
         GNH3     = GNH3 + CNH4HS4  ! evaporate NH3 to gas phase
         EXS4     = EXS4  - CNH4HS4 ! reduce excess
         CNH4HS4  = ZERO            ! zero salt concentration
      ENDIF
C
      IF ((CNH42S4).GT.(EXS4)) THEN     ! Adjust (NH4)2SO4(s)
         CNH42S4  = CNH42S4- EXS4   ! more solid than excess
         GNH3     = GNH3 + 2.d0*EXS4! evaporate NH3 to gas phase
         GOTO 50
      ELSE                          ! less solid than excess
         GNH3     = GNH3+2.d0*CNH42S4 ! evaporate NH3 to gas phase
         EXS4     = EXS4  - CNH42S4 ! reduce excess
         CNH42S4  = ZERO            ! zero salt concentration
      ENDIF
C
C *** RETURN **********************************************************
C
 50   RETURN
      END
      
C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION GETASR
C *** CALCULATES THE LIMITING NH4+/SO4 RATIO OF A SULFATE POOR SYSTEM
C     (i.e. SULFATE RATIO = 2.0) FOR GIVEN SO4 LEVEL AND RH
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION GETASR (SO4I, RHI)
      IMPLICIT NONE
      INTEGER NSO4S, NRHS, NASRD
      PARAMETER (NSO4S=14, NRHS=20, NASRD=NSO4S*NRHS)
      DOUBLE PRECISION WF, ASRAT, ASSO4
      COMMON /ASRC/ ASRAT(NASRD), ASSO4(NSO4S)
      DOUBLE PRECISION SO4I, RHI, RAT
      INTEGER IA1, A1, INDS, INDR, INDSL, INDSH, IPOSL, IPOSH
CCC
CCC *** SOLVE USING FULL COMPUTATIONS, NOT LOOK-UP TABLES **************
CCC
CCC         W(2) = WAER(2)
CCC         W(3) = WAER(2)*2.0001D0
CCC         CALL CALCA2
CCC         SULRATW = MOLAL(3)/WAER(2)
CCC         CALL INIT1 (WI, RHI, TEMPI)   ! Re-initialize COMMON BLOCK
C
C *** CALCULATE INDICES ************************************************
C
      RAT    = SO4I/1.E-9    
      A1     = INT(LOG10(RAT))                   ! Magnitude of RAT
      IA1    = INT(RAT/2.5/10.0**A1)
C
      INDS   = INT(4.0*A1 + MIN(IA1,4))
      INDS   = MIN(MAX(0, INDS), NSO4S-1) + 1     ! SO4 component of IPOS
C
      INDR   = INT(99.0-RHI*100.0) + 1
      INDR   = MIN(MAX(1, INDR), NRHS)            ! RH component of IPOS
C
C *** GET VALUE AND RETURN *********************************************
C
      INDSL  = INDS
      INDSH  = MIN(INDSL+1, NSO4S)
      IPOSL  = (INDSL-1)*NRHS + INDR              ! Low position in array
      IPOSH  = (INDSH-1)*NRHS + INDR              ! High position in array
C
      WF     = (SO4I-ASSO4(INDSL))/(ASSO4(INDSH)-ASSO4(INDSL) + 1.D-7)
C      WF     = MIN(MAX((WF), 0.0), 1.0)
      IF ((WF) .LT. 0.D0) THEN
         WF = 0.d0
      ELSEIF ((WF) .GT. 1.D0) THEN
         WF = 1.d0
      ENDIF
C
      GETASR = WF*ASRAT(IPOSH) + (1.0-WF)*ASRAT(IPOSL)
C
C *** END OF FUNCTION GETASR *******************************************
C
      RETURN
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** BLOCK DATA AERSR
C *** CONTAINS DATA FOR AEROSOL SULFATE RATIO ARRAY NEEDED IN FUNCTION 
C     GETASR
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      BLOCK DATA AERSR
      IMPLICIT NONE
      INTEGER NSO4S, NRHS, NASRD, I
      PARAMETER (NSO4S=14, NRHS=20, NASRD=NSO4S*NRHS)
      DOUBLE PRECISION ASRAT, ASSO4
      COMMON /ASRC/ ASRAT(NASRD), ASSO4(NSO4S)
C
      DATA ASSO4/1.0E-9, 2.5E-9, 5.0E-9, 7.5E-9, 1.0E-8,
     &           2.5E-8, 5.0E-8, 7.5E-8, 1.0E-7, 2.5E-7, 
     &           5.0E-7, 7.5E-7, 1.0E-6, 5.0E-6/
C
      DATA (ASRAT(I), I=1,280)/
     & 1.020464D0, 0.9998130D0,0.9960167D0,0.9984423D0, 1.004004D0,
     & 1.010885D0,  1.018356D0, 1.026726D0, 1.034268D0, 1.043846D0,
     & 1.052933D0,  1.062230D0, 1.062213D0, 1.080050D0, 1.088350D0,
     & 1.096603D0,  1.104289D0, 1.111745D0, 1.094662D0, 1.121594D0,
     & 1.268909D0,  1.242444D0, 1.233815D0, 1.232088D0, 1.234020D0,
     & 1.238068D0,  1.243455D0, 1.250636D0, 1.258734D0, 1.267543D0,
     & 1.276948D0,  1.286642D0, 1.293337D0, 1.305592D0, 1.314726D0,
     & 1.323463D0,  1.333258D0, 1.343604D0, 1.344793D0, 1.355571D0,
     & 1.431463D0,  1.405204D0, 1.395791D0, 1.393190D0, 1.394403D0,
     & 1.398107D0,  1.403811D0, 1.411744D0, 1.420560D0, 1.429990D0,
     & 1.439742D0,  1.449507D0, 1.458986D0, 1.468403D0, 1.477394D0,
     & 1.487373D0,  1.495385D0, 1.503854D0, 1.512281D0, 1.520394D0,
     & 1.514464D0,  1.489699D0, 1.480686D0, 1.478187D0, 1.479446D0,
     & 1.483310D0,  1.489316D0, 1.497517D0, 1.506501D0, 1.515816D0,
     & 1.524724D0,  1.533950D0, 1.542758D0, 1.551730D0, 1.559587D0,
     & 1.568343D0,  1.575610D0, 1.583140D0, 1.590440D0, 1.596481D0,
     & 1.567743D0,  1.544426D0, 1.535928D0, 1.533645D0, 1.535016D0,
     & 1.539003D0,  1.545124D0, 1.553283D0, 1.561886D0, 1.570530D0,
     & 1.579234D0,  1.587813D0, 1.595956D0, 1.603901D0, 1.611349D0,
     & 1.618833D0,  1.625819D0, 1.632543D0, 1.639032D0, 1.645276D0,
     & 1.707390D0,  1.689553D0, 1.683198D0, 1.681810D0, 1.683490D0,
     & 1.687477D0,  1.693148D0, 1.700084D0, 1.706917D0, 1.713507D0,
     & 1.719952D0,  1.726190D0, 1.731985D0, 1.737544D0, 1.742673D0,
     & 1.747756D0,  1.752431D0, 1.756890D0, 1.761141D0, 1.765190D0,
     & 1.785657D0,  1.771851D0, 1.767063D0, 1.766229D0, 1.767901D0,
     & 1.771455D0,  1.776223D0, 1.781769D0, 1.787065D0, 1.792081D0,
     & 1.796922D0,  1.801561D0, 1.805832D0, 1.809896D0, 1.813622D0,
     & 1.817292D0,  1.820651D0, 1.823841D0, 1.826871D0, 1.829745D0,
     & 1.822215D0,  1.810497D0, 1.806496D0, 1.805898D0, 1.807480D0,
     & 1.810684D0,  1.814860D0, 1.819613D0, 1.824093D0, 1.828306D0,
     & 1.832352D0,  1.836209D0, 1.839748D0, 1.843105D0, 1.846175D0,
     & 1.849192D0,  1.851948D0, 1.854574D0, 1.857038D0, 1.859387D0,
     & 1.844588D0,  1.834208D0, 1.830701D0, 1.830233D0, 1.831727D0,
     & 1.834665D0,  1.838429D0, 1.842658D0, 1.846615D0, 1.850321D0,
     & 1.853869D0,  1.857243D0, 1.860332D0, 1.863257D0, 1.865928D0,
     & 1.868550D0,  1.870942D0, 1.873208D0, 1.875355D0, 1.877389D0,
     & 1.899556D0,  1.892637D0, 1.890367D0, 1.890165D0, 1.891317D0,
     & 1.893436D0,  1.896036D0, 1.898872D0, 1.901485D0, 1.903908D0,
     & 1.906212D0,  1.908391D0, 1.910375D0, 1.912248D0, 1.913952D0,
     & 1.915621D0,  1.917140D0, 1.918576D0, 1.919934D0, 1.921220D0,
     & 1.928264D0,  1.923245D0, 1.921625D0, 1.921523D0, 1.922421D0,
     & 1.924016D0,  1.925931D0, 1.927991D0, 1.929875D0, 1.931614D0,
     & 1.933262D0,  1.934816D0, 1.936229D0, 1.937560D0, 1.938769D0,
     & 1.939951D0,  1.941026D0, 1.942042D0, 1.943003D0, 1.943911D0,
     & 1.941205D0,  1.937060D0, 1.935734D0, 1.935666D0, 1.936430D0,
     & 1.937769D0,  1.939359D0, 1.941061D0, 1.942612D0, 1.944041D0,
     & 1.945393D0,  1.946666D0, 1.947823D0, 1.948911D0, 1.949900D0,
     & 1.950866D0,  1.951744D0, 1.952574D0, 1.953358D0, 1.954099D0,
     & 1.948985D0,  1.945372D0, 1.944221D0, 1.944171D0, 1.944850D0,
     & 1.946027D0,  1.947419D0, 1.948902D0, 1.950251D0, 1.951494D0,
     & 1.952668D0,  1.953773D0, 1.954776D0, 1.955719D0, 1.956576D0,
     & 1.957413D0,  1.958174D0, 1.958892D0, 1.959571D0, 1.960213D0,
     & 1.977193D0,  1.975540D0, 1.975023D0, 1.975015D0, 1.975346D0,
     & 1.975903D0,  1.976547D0, 1.977225D0, 1.977838D0, 1.978401D0,
     & 1.978930D0,  1.979428D0, 1.979879D0, 1.980302D0, 1.980686D0,
     & 1.981060D0,  1.981401D0, 1.981722D0, 1.982025D0, 1.982312D0/
C
C *** END OF BLOCK DATA AERSR ******************************************
C
       END
      
       
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCHA
C *** CALCULATES CHLORIDES SPECIATION
C
C     HYDROCHLORIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES,  
C     AND DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE 
C     HYDROCHLORIC ACID DISSOLVED IS CALCULATED FROM THE 
C     HCL(G) <-> (H+) + (CL-) 
C     EQUILIBRIUM, USING THE (H+) FROM THE SULFATES.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCHA
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION KAPA, X, DELT, ALFA, DIAK
CC      CHARACTER*40 errinf
C
C *** CALCULATE HCL DISSOLUTION *****************************************
C
      X    = W(5) 
      DELT = 0.0d0
      IF ((WATER).GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK3*R*TEMP*(WATER/GAMA(11))**2.0
         DIAK = SQRT( (KAPA+ALFA)**2.0 + 4.0*ALFA*X)
         DELT = 0.5*(-(KAPA+ALFA) + DIAK)
CC         IF (DELT/KAPA.GT.0.1d0) THEN
CC            WRITE (ERRINF,'(1PE10.3)') DELT/KAPA*100.0
CC            CALL PUSHERR (0033, ERRINF)    
CC         ENDIF
      ENDIF
C
C *** CALCULATE HCL SPECIATION IN THE GAS PHASE *************************
C
      GHCL     = MAX(X-DELT, 0.0d0)  ! GAS HCL
C
C *** CALCULATE HCL SPECIATION IN THE LIQUID PHASE **********************
C
      MOLAL(4) = DELT                ! CL-
      MOLAL(1) = MOLAL(1) + DELT     ! H+ 
C 
      RETURN
C
C *** END OF SUBROUTINE CALCHA ******************************************
C
      END





C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCHAP
C *** CALCULATES CHLORIDES SPECIATION
C
C     HYDROCHLORIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES, 
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. 
C     THE HYDROCHLORIC ACID DISSOLVED IS CALCULATED FROM THE 
C     HCL(G) -> HCL(AQ)   AND  HCL(AQ) ->  (H+) + (CL-) 
C     EQUILIBRIA, USING (H+) FROM THE SULFATES.
C
C     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOVER
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCHAP
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: ALFA, DELT
C
C *** IS THERE A LIQUID PHASE? ******************************************
C
      IF ((WATER).LE.TINY) RETURN
C
C *** CALCULATE HCL SPECIATION IN THE GAS PHASE *************************
C
      CALL CALCCLAQ (MOLAL(4), MOLAL(1), DELT)
      ALFA     = XK3*R*TEMP*(WATER/GAMA(11))**2.0
      GASAQ(3) = DELT
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(4) = MOLAL(4) - DELT
      GHCL     = MOLAL(1)*MOLAL(4)/ALFA
C 
      RETURN
C
C *** END OF SUBROUTINE CALCHAP *****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNA
C *** CALCULATES NITRATES SPECIATION
C
C     NITRIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES, THAT 
C     DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE NITRIC
C     ACID DISSOLVED IS CALCULATED FROM THE HNO3(G) -> (H+) + (NO3-) 
C     EQUILIBRIUM, USING THE (H+) FROM THE SULFATES.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNA
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: ALFA, DELT, KAPA, DIAK
CC      CHARACTER*40 errinf
C
C *** CALCULATE HNO3 DISSOLUTION ****************************************
C
      X    = W(4) 
      DELT = 0.0d0
      IF ((WATER).GT.TINY) THEN
         KAPA = MOLAL(1)
         ALFA = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         DIAK = SQRT( (KAPA+ALFA)**2.0 + 4.0*ALFA*X)
         DELT = 0.5*(-(KAPA+ALFA) + DIAK)
CC         IF (DELT/KAPA.GT.0.1d0) THEN
CC            WRITE (ERRINF,'(1PE10.3)') DELT/KAPA*100.0
CC            CALL PUSHERR (0019, ERRINF)    ! WARNING ERROR: NO SOLUTION
CC         ENDIF
      ENDIF
C
C *** CALCULATE HNO3 SPECIATION IN THE GAS PHASE ************************
C
      GHNO3    = MAX(X-DELT, 0.0d0)  ! GAS HNO3
C
C *** CALCULATE HNO3 SPECIATION IN THE LIQUID PHASE *********************
C
      MOLAL(7) = DELT                ! NO3-
      MOLAL(1) = MOLAL(1) + DELT     ! H+ 
C 
      RETURN
C
C *** END OF SUBROUTINE CALCNA ******************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNAP
C *** CALCULATES NITRATES SPECIATION
C
C     NITRIC ACID IN THE LIQUID PHASE IS ASSUMED A MINOR SPECIES, THAT 
C     DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM. THE NITRIC
C     ACID DISSOLVED IS CALCULATED FROM THE HNO3(G) -> HNO3(AQ) AND
C     HNO3(AQ) -> (H+) + (CL-) EQUILIBRIA, USING (H+) FROM THE SULFATES.
C
C     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOVER
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNAP
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: ALFA, DELT
C
C *** IS THERE A LIQUID PHASE? ******************************************
C
      IF ((WATER).LE.TINY) RETURN
C
C *** CALCULATE HNO3 SPECIATION IN THE GAS PHASE ************************
C
      CALL CALCNIAQ (MOLAL(7), MOLAL(1), DELT)
      ALFA     = XK4*R*TEMP*(WATER/GAMA(10))**2.0
      GASAQ(3) = DELT
      MOLAL(1) = MOLAL(1) - DELT
      MOLAL(7) = MOLAL(7) - DELT
      GHNO3    = MOLAL(1)*MOLAL(7)/ALFA
      
C      write (*,*) ALFA, MOLAL(1), MOLAL(7), GHNO3, DELT
C 
      RETURN
C
C *** END OF SUBROUTINE CALCNAP *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNH3
C *** CALCULATES AMMONIA IN GAS PHASE
C
C     AMMONIA IN THE GAS PHASE IS ASSUMED A MINOR SPECIES, THAT 
C     DOES NOT SIGNIFICANTLY PERTURB THE AEROSOL EQUILIBRIUM. 
C     AMMONIA GAS IS CALCULATED FROM THE NH3(g) + (H+)(l) <==> (NH4+)(l)
C     EQUILIBRIUM, USING (H+), (NH4+) FROM THE AEROSOL SOLUTION.
C
C     THIS IS THE VERSION USED BY THE DIRECT PROBLEM
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNH3
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: BB, CC, DIAK, PSI
C
C *** IS THERE A LIQUID PHASE? ******************************************
C
      IF ((WATER).LE.TINY) RETURN
C
C *** CALCULATE NH3 SUBLIMATION *****************************************
C
      A1   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      CHI1 = MOLAL(3)
      CHI2 = MOLAL(1)
C
      BB   =(CHI2 + ONE/A1)          ! a=1; b!=1; c!=1 
      CC   =-CHI1/A1             
      DIAK = SQRT(BB*BB - 4.D0*CC)   ! Always > 0
      PSI  = 0.5*(-BB + DIAK)        ! One positive root
      PSI  = MAX(MIN(PSI,CHI1), TINY)! Constrict in acceptible range
C
C *** CALCULATE NH3 SPECIATION IN THE GAS PHASE *************************
C
      GNH3     = PSI                 ! GAS HNO3
C
C *** CALCULATE NH3 AFFECT IN THE LIQUID PHASE **************************
C
      MOLAL(3) = CHI1 - PSI          ! NH4+
      MOLAL(1) = CHI2 + PSI          ! H+ 
C 
      RETURN
C
C *** END OF SUBROUTINE CALCNH3 *****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNH3P
C *** CALCULATES AMMONIA IN GAS PHASE
C
C     AMMONIA GAS IS CALCULATED FROM THE NH3(g) + (H+)(l) <==> (NH4+)(l)
C     EQUILIBRIUM, USING (H+), (NH4+) FROM THE AEROSOL SOLUTION.
C
C     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOLVER
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNH3P
      INCLUDE 'isrpia_b.inc'
C
C *** IS THERE A LIQUID PHASE? ******************************************
C
      IF ((WATER).LE.TINY) RETURN
C
C *** CALCULATE NH3 GAS PHASE CONCENTRATION *****************************
C
      A1   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      GNH3 = MOLAL(3)/MOLAL(1)/A1
C 
      RETURN
C
C *** END OF SUBROUTINE CALCNH3P ****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNHA
C
C     THIS SUBROUTINE CALCULATES THE DISSOLUTION OF HCL, HNO3 AT
C     THE PRESENCE OF (H,SO4). HCL, HNO3 ARE CONSIDERED MINOR SPECIES,
C     THAT DO NOT SIGNIFICANTLY AFFECT THE EQUILIBRIUM POINT.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNHA
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION M1, M2, M3, DELCL, DELNO, OMEGA
      CHARACTER(40) errinf
C
C *** SPECIAL CASE; WATER=ZERO ******************************************
C
      IF ((WATER).LE.TINY) THEN
c wz
         GOTO 55
C
C *** SPECIAL CASE; HCL=HNO3=ZERO ***************************************
C
      ELSEIF ((W(5)).LE.TINY .AND. (W(4)).LE.TINY) THEN
         GOTO 60
C
C *** SPECIAL CASE; HCL=ZERO ********************************************
C
      ELSE IF ((W(5)).LE.TINY) THEN
         CALL CALCNA              ! CALL HNO3 DISSOLUTION ROUTINE
         GOTO 60
C
C *** SPECIAL CASE; HNO3=ZERO *******************************************
C
      ELSE IF ((W(4)).LE.TINY) THEN
         CALL CALCHA              ! CALL HCL DISSOLUTION ROUTINE
         GOTO 60
      ENDIF
C
C *** CALCULATE EQUILIBRIUM CONSTANTS ***********************************
C
      A3 = XK4*R*TEMP*(WATER/GAMA(10))**2.0   ! HNO3
      A4 = XK3*R*TEMP*(WATER/GAMA(11))**2.0   ! HCL
C
C *** CALCULATE CUBIC EQUATION COEFFICIENTS *****************************
C
      DELCL = ZERO
      DELNO = ZERO
C
      OMEGA = MOLAL(1)       ! H+
      CHI3  = W(4)           ! HNO3
      CHI4  = W(5)           ! HCL
C
      C1    = A3*CHI3
      C2    = A4*CHI4
      C3    = A3 - A4
C
      M1    = (C1 + C2 + (OMEGA+A4)*C3)/C3
      M2    = ((OMEGA+A4)*C2 - A4*C3*CHI4)/C3
      M3    =-A4*C2*CHI4/C3
C
C *** CALCULATE ROOTS ***************************************************
C
      CALL POLY3 (M1, M2, M3, DELCL, ISLV) ! HCL DISSOLUTION
      IF (ISLV.NE.0) THEN
         DELCL = TINY       ! TINY AMOUNTS OF HCL ASSUMED WHEN NO ROOT 
         WRITE (ERRINF,'(1PE10.1)') TINY
         CALL PUSHERR (0022, ERRINF)    ! WARNING ERROR: NO SOLUTION
      ENDIF
      DELCL = MIN(DELCL, CHI4)
C
      DELNO = C1*DELCL/(C2 + C3*DELCL)  
      DELNO = MIN(DELNO, CHI3)
C
      IF ((DELCL).LT.ZERO .OR. (DELNO).LT.ZERO .OR.
     &   (DELCL).GT.(CHI4) .OR. (DELNO).GT.(CHI3)) THEN
         DELCL = TINY  ! TINY AMOUNTS OF HCL ASSUMED WHEN NO ROOT 
         DELNO = TINY
         WRITE (ERRINF,'(1PE10.1)') TINY
         CALL PUSHERR (0022, ERRINF)    ! WARNING ERROR: NO SOLUTION
      ENDIF
CCC
CCC *** COMPARE DELTA TO TOTAL H+ ; ESTIMATE EFFECT TO HSO4 ***************
CCC
CC      IF ((DELCL+DELNO)/MOLAL(1).GT.0.1d0) THEN
CC         WRITE (ERRINF,'(1PE10.3)') (DELCL+DELNO)/MOLAL(1)*100.0
CC         CALL PUSHERR (0021, ERRINF)   
CC      ENDIF
C
C *** EFFECT ON LIQUID PHASE ********************************************
C
50    MOLAL(1) = MOLAL(1) + (DELNO+DELCL)  ! H+   CHANGE
      MOLAL(4) = MOLAL(4) + DELCL          ! CL-  CHANGE
      MOLAL(7) = MOLAL(7) + DELNO          ! NO3- CHANGE
C
C *** EFFECT ON GAS PHASE ***********************************************
C
55    GHCL     = MAX((W(5) - MOLAL(4)), TINY)
      GHNO3    = MAX((W(4) - MOLAL(7)), TINY)
C
60    RETURN
C
C *** END OF SUBROUTINE CALCNHA *****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNHP
C
C     THIS SUBROUTINE CALCULATES THE GAS PHASE NITRIC AND HYDROCHLORIC
C     ACID. CONCENTRATIONS ARE CALCULATED FROM THE DISSOLUTION 
C     EQUILIBRIA, USING (H+), (Cl-), (NO3-) IN THE AEROSOL PHASE.
C
C     THIS IS THE VERSION USED BY THE INVERSE PROBLEM SOLVER
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNHP
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION :: DELT
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
      RETURN
C
C *** END OF SUBROUTINE CALCNHP *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCAMAQ
C *** THIS SUBROUTINE CALCULATES THE NH3(aq) GENERATED FROM (H,NH4+).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCAMAQ (NH4I, OHI, DELT)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION NH4I, OHI, OM1, OM2, BB, CC, DD, DEL1, DEL2
CC      CHARACTER*40 errinf
C
C *** EQUILIBRIUM CONSTANTS
C
      A22  = XK22/XKW/WATER*(GAMA(8)/GAMA(9))**2. ! GAMA(NH3) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
C
C *** FIND ROOT
C
      OM1  = NH4I          
      OM2  = OHI
      BB   =-(OM1+OM2+A22*AKW)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)
C
C *** GET APPROPRIATE ROOT.
C
      IF ((DEL1).LT.ZERO) THEN                 
         IF ((DEL2).GT.(NH4I) .OR. (DEL2).GT.(OHI)) THEN
            DELT = ZERO
         ELSE
            DELT = DEL2
         ENDIF
      ELSE
         DELT = DEL1
      ENDIF
CC
CC *** COMPARE DELTA TO TOTAL NH4+ ; ESTIMATE EFFECT *********************
CC
CC      IF (DELTA/HYD.GT.0.1d0) THEN
CC         WRITE (ERRINF,'(1PE10.3)') DELTA/HYD*100.0
CC         CALL PUSHERR (0020, ERRINF)
CC      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCAMAQ ****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCAMAQ2
C
C     THIS SUBROUTINE CALCULATES THE NH3(aq) GENERATED FROM (H,NH4+).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCAMAQ2 (GGNH3, NH4I, OHI, NH3AQ)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION NH4I, NH3AQ, ALF1, ALF2, BB, CC, DEL, OHI
C
C *** EQUILIBRIUM CONSTANTS
C
      A22  = XK22/XKW/WATER*(GAMA(8)/GAMA(9))**2. ! GAMA(NH3) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
C
C *** FIND ROOT
C
      ALF1 = NH4I - GGNH3
      ALF2 = GGNH3
      BB   = ALF1 + A22*AKW
      CC   =-A22*AKW*ALF2
      DEL  = 0.5D0*(-BB + SQRT(BB*BB-4.D0*CC))
C
C *** ADJUST CONCENTRATIONS
C
      NH4I  = ALF1 + DEL
      OHI   = DEL
      IF ((OHI).LE.TINY) OHI = SQRT(AKW)   ! If solution is neutral.
      NH3AQ = ALF2 - DEL 
C
      RETURN
C
C *** END OF SUBROUTINE CALCAMAQ2 ****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCCLAQ
C
C     THIS SUBROUTINE CALCULATES THE HCL(aq) GENERATED FROM (H+,CL-).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCCLAQ (CLI, HI, DELT)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION CLI, HI, OM1, OM2, BB, CC, DD, DELT
C
C *** EQUILIBRIUM CONSTANTS
C
      A32  = XK32*WATER/(GAMA(11))**2. ! GAMA(HCL) ASSUMED 1
C
C *** FIND ROOT
C
      OM1  = CLI          
      OM2  = HI
      BB   =-(OM1+OM2+A32)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)
C
C *** GET APPROPRIATE ROOT.
C
      IF ((DEL1).LT.ZERO) THEN                 
         IF ((DEL2).LT.ZERO .OR. (DEL2).GT.(CLI) .OR.
     +       (DEL2).GT.(HI)) THEN
            DELT = ZERO
         ELSE
            DELT = DEL2
         ENDIF
      ELSE
         DELT = DEL1
      ENDIF
C      WRITE(*,*) 'DELT: ',DELT
C      PAUSE
C
      RETURN
C
C *** END OF SUBROUTINE CALCCLAQ ****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCCLAQ2
C
C     THIS SUBROUTINE CALCULATES THE HCL(aq) GENERATED FROM (H+,CL-).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCCLAQ2 (GGCL, CLI, HI, CLAQ)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION GGCL, CLI, HI, CLAQ, ALF1, ALF2, DEL1
C
C *** EQUILIBRIUM CONSTANTS
C
      A32  = XK32*WATER/(GAMA(11))**2. ! GAMA(HCL) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
C
C *** FIND ROOT
C
      ALF1  = CLI - GGCL
      ALF2  = GGCL
      COEF  = (ALF1+A32)
      DEL1  = 0.5*(-COEF + SQRT(COEF*COEF+4.D0*A32*ALF2))
C
C *** CORRECT CONCENTRATIONS
C
      CLI  = ALF1 + DEL1
      HI   = DEL1
      IF ((HI).LE.TINY) HI = SQRT(AKW)   ! If solution is neutral.
      CLAQ = ALF2 - DEL1
C
      RETURN
C
C *** END OF SUBROUTINE CALCCLAQ2 ****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNIAQ
C
C     THIS SUBROUTINE CALCULATES THE HNO3(aq) GENERATED FROM (H,NO3-).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNIAQ (NO3I, HI, DELT)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION NO3I, HI, DELT
      DOUBLE PRECISION OM1, OM2, BB, CC, DD, DEL1, DEL2
C
C *** EQUILIBRIUM CONSTANTS
C
      A42  = XK42*WATER/(GAMA(10))**2. ! GAMA(HNO3) ASSUMED 1
C
C *** FIND ROOT
C
      OM1  = NO3I          
      OM2  = HI
      BB   =-(OM1+OM2+A42)
      CC   = OM1*OM2
      DD   = SQRT(BB*BB-4.D0*CC)

      DEL1 = 0.5D0*(-BB - DD)
      DEL2 = 0.5D0*(-BB + DD)
C
C *** GET APPROPRIATE ROOT.
C
      IF ((DEL1).LT.ZERO .OR. (DEL1).GT.(HI) .OR.
     &    (DEL1).GT.(NO3I)) THEN
         DELT = ZERO
      ELSE
         DELT = DEL1
         RETURN
      ENDIF
C
      IF ((DEL2).LT.ZERO .OR. (DEL2).GT.(NO3I) .OR. 
     &    (DEL2).GT.(HI)) THEN
         DELT = ZERO
      ELSE
         DELT = DEL2
      ENDIF
C      WRITE(*,*) 'DELT: ',DELT
C      PAUSE
C
      RETURN
C
C *** END OF SUBROUTINE CALCNIAQ ****************************************
C
      END



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCNIAQ2
C
C     THIS SUBROUTINE CALCULATES THE UNDISSOCIATED HNO3(aq)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION GGNO3, NO3I, HI, NO3AQ
      DOUBLE PRECISION OM1, OM2, BB, CC, DD, ALF1, ALF2, ALF3, DEL1
C
C *** EQUILIBRIUM CONSTANTS
C
      A42  = XK42*WATER/(GAMA(10))**2. ! GAMA(HNO3) ASSUMED 1
      AKW  = XKW *RH*WATER*WATER
C
C *** FIND ROOT
C
      ALF1  = NO3I - GGNO3
      ALF2  = GGNO3
      ALF3  = HI
C
      BB    = ALF3 + ALF1 + A42
      CC    = ALF3*ALF1 - A42*ALF2
      DEL1  = 0.5*(-BB + SQRT(BB*BB-4.D0*CC))
C
C *** CORRECT CONCENTRATIONS
C
      NO3I  = ALF1 + DEL1
      HI    = ALF3 + DEL1
      IF ((HI).LE.TINY) HI = SQRT(AKW)   ! If solution is neutral.
      NO3AQ = ALF2 - DEL1
C
      RETURN
C
C *** END OF SUBROUTINE CALCNIAQ2 ****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCMR
C *** THIS SUBROUTINE CALCULATES:
C     1. ION PAIR CONCENTRATIONS (FROM [MOLAR] ARRAY)
C     2. WATER CONTENT OF LIQUID AEROSOL PHASE (FROM ZSR CORRELATION)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCMR
      INCLUDE 'isrpia_b.inc'
C      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
C     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7, PSI8,
C     &               A1,   A2,   A3,   A4,   A5,   A6,   A7,   A8
      CHARACTER SC*1
      DOUBLE PRECISION HSO4I, SO4I, AML5, TOTS4, FRNO3, FRCL, FRNH4
C
C *** CALCULATE ION PAIR CONCENTRATIONS ACCORDING TO SPECIFIC CASE ****
C 
      SC =SCASE(1:1)                   ! SULRAT & SODRAT case
C
C *** NH4-SO4 SYSTEM ; SULFATE POOR CASE
C
      IF (SC.EQ.'A') THEN      
         MOLALR(4) = MOLAL(5)+MOLAL(6) ! (NH4)2SO4 - CORRECT FOR SO4 TO HSO4
C
C *** NH4-SO4 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
C
      ELSE IF (SC.EQ.'B') THEN
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
C *** NH4-SO4 SYSTEM ; SULFATE RICH CASE ; FREE ACID 
C
      ELSE IF (SC.EQ.'C') THEN
         MOLALR(9) = MOLAL(3)                     ! NH4HSO4
         MOLALR(7) = MAX(W(2)-W(3), ZERO)         ! H2SO4
C
C *** NH4-SO4-NO3 SYSTEM ; SULFATE POOR CASE
C
      ELSE IF (SC.EQ.'D') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6)          ! (NH4)2SO4
         AML5      = MOLAL(3)-2.D0*MOLALR(4)      ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,MOLAL(7)), ZERO)! NH4NO3 = MIN("free", NO3)
C
C *** NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
C
      ELSE IF (SC.EQ.'E') THEN      
         SO4I  = MAX(MOLAL(5)-MOLAL(1),ZERO)      ! FROM HSO4 DISSOCIATION 
         HSO4I = MOLAL(6)+MOLAL(1)              
         IF ((SO4I).LT.(HSO4I)) THEN                
            MOLALR(13) = SO4I                     ! [LC] = [SO4] 
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)    ! NH4HSO4
         ELSE                                   
            MOLALR(13) = HSO4I                    ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)    ! (NH4)2SO4
         ENDIF
C
C *** NH4-SO4-NO3 SYSTEM ; SULFATE RICH CASE ; FREE ACID
C
      ELSE IF (SC.EQ.'F') THEN      
         MOLALR(9) = MOLAL(3)                              ! NH4HSO4
         MOLALR(7) = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3),ZERO)  ! H2SO4
C
C *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM POOR CASE
C
      ELSE IF (SC.EQ.'G') THEN      
         MOLALR(2) = 0.5*MOLAL(2)                          ! NA2SO4
         TOTS4     = MOLAL(5)+MOLAL(6)                     ! Total SO4
         MOLALR(4) = MAX(TOTS4 - MOLALR(2), ZERO)          ! (NH4)2SO4
         FRNH4     = MAX(MOLAL(3) - 2.D0*MOLALR(4), ZERO)
         MOLALR(5) = MIN(MOLAL(7),FRNH4)                   ! NH4NO3
         FRNH4     = MAX(FRNH4 - MOLALR(5), ZERO)
         MOLALR(6) = MIN(MOLAL(4), FRNH4)                  ! NH4CL
C
C *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE POOR ; SODIUM RICH CASE
C *** RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
C
      ELSE IF (SC.EQ.'H') THEN      
         MOLALR(1) = PSI7                                  ! NACL 
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI8                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         FRNO3     = MAX(MOLAL(7) - MOLALR(3), ZERO)       ! "FREE" NO3
         FRCL      = MAX(MOLAL(4) - MOLALR(1), ZERO)       ! "FREE" CL
         MOLALR(5) = MIN(MOLAL(3),FRNO3)                   ! NH4NO3
         FRNH4     = MAX(MOLAL(3) - MOLALR(5), ZERO)       ! "FREE" NH3
         MOLALR(6) = MIN(FRCL, FRNH4)                      ! NH4CL
C
C *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; NO FREE ACID
C *** RETREIVE DISSOLVED SALTS DIRECTLY FROM COMMON BLOCK /SOLUT/
C
      ELSE IF (SC.EQ.'I') THEN      
         MOLALR(04) = PSI5                                 ! (NH4)2SO4
         MOLALR(02) = PSI4                                 ! NA2SO4
         MOLALR(09) = PSI1                                 ! NH4HSO4
         MOLALR(12) = PSI3                                 ! NAHSO4
         MOLALR(13) = PSI2                                 ! LC
C
C *** NA-NH4-SO4-NO3-CL SYSTEM ; SULFATE RICH CASE ; FREE ACID
C
      ELSE IF (SC.EQ.'J') THEN      
         MOLALR(09) = MOLAL(3)                             ! NH4HSO4
         MOLALR(12) = MOLAL(2)                             ! NAHSO4
         MOLALR(07) = MOLAL(5)+MOLAL(6)-MOLAL(3)-MOLAL(2)  ! H2SO4
         MOLALR(07) = MAX(MOLALR(07),ZERO)
C
C ======= REVERSE PROBLEMS ===========================================
C
C *** NH4-SO4-NO3 SYSTEM ; SULFATE POOR CASE
C
      ELSE IF (SC.EQ.'N') THEN      
         MOLALR(4) = MOLAL(5) + MOLAL(6)          ! (NH4)2SO4
         AML5      = WAER(3)-2.D0*MOLALR(4)       ! "free" NH4
         MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO) ! NH4NO3 = MIN("free", NO3)
C
C *** NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM POOR CASE
C
      ELSE IF (SC.EQ.'Q') THEN      
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(4) = PSI6                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL
C
C *** NH4-SO4-NO3-NA-CL SYSTEM ; SULFATE POOR, SODIUM RICH CASE
C
      ELSE IF (SC.EQ.'R') THEN      
         MOLALR(1) = PSI3                                  ! NACL 
         MOLALR(2) = PSI1                                  ! NA2SO4
         MOLALR(3) = PSI2                                  ! NANO3
         MOLALR(4) = ZERO                                  ! (NH4)2SO4
         MOLALR(5) = PSI5                                  ! NH4NO3
         MOLALR(6) = PSI4                                  ! NH4CL
C
C *** UNKNOWN CASE
C
      ELSE
         CALL PUSHERR (1001, ' ') ! FATAL ERROR: CASE NOT SUPPORTED 
      ENDIF
C
C *** CALCULATE WATER CONTENT ; ZSR CORRELATION ***********************
C
      WATER = ZERO
      DO I=1,NPAIR
         WATER = WATER + MOLALR(I)/M0(I)
      ENDDO
      WATER = MAX(WATER, TINY)
C
      RETURN
C
C *** END OF SUBROUTINE CALCMR ******************************************
C
      END
C
CC=======================================================================
CC
CC *** ISORROPIA CODE
CC *** SUBROUTINE CALCMDRH
CC
CC     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
CC     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
CC     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
CC     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
CC
CC *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
CC *** GEORGIA INSTITUTE OF TECHNOLOGY
CC *** WRITTEN BY ATHANASIOS NENES
CC *** UPDATED BY CHRISTOS FOUNTOUKIS
CC
CC=======================================================================
CC
C      SUBROUTINE CALCMDRH (RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE)
C      INCLUDE 'isrpia_b.inc'
C      EXTERNAL DRYCASE, LIQCASE
CC
CC *** FIND WEIGHT FACTOR **********************************************
CC
C      IF (WFTYP.EQ.0) THEN
C         WF = ONE
C      ELSEIF (WFTYP.EQ.1) THEN
C         WF = 0.5D0
C      ELSE
C         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
C      ENDIF
C      ONEMWF  = ONE - WF
CC
CC *** FIND FIRST SECTION ; DRY ONE ************************************
CC
C      CALL DRYCASE
C      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL
CC
C      CNH42SO = CNH42S4                  ! FIRST (DRY) SOLUTION
C      CNH4HSO = CNH4HS4
C      CLCO    = CLC 
C      CNH4N3O = CNH4NO3
C      CNH4CLO = CNH4CL
C      CNA2SO  = CNA2SO4
C      CNAHSO  = CNAHSO4
C      CNANO   = CNANO3
C      CNACLO  = CNACL
C      GNH3O   = GNH3
C      GHNO3O  = GHNO3
C      GHCLO   = GHCL
CC
CC *** FIND SECOND SECTION ; DRY & LIQUID ******************************
CC
C      CNH42S4 = ZERO
C      CNH4HS4 = ZERO
C      CLC     = ZERO
C      CNH4NO3 = ZERO
C      CNH4CL  = ZERO
C      CNA2SO4 = ZERO
C      CNAHSO4 = ZERO
C      CNANO3  = ZERO
C      CNACL   = ZERO
C      GNH3    = ZERO
C      GHNO3   = ZERO
C      GHCL    = ZERO
C      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION
CC
CC *** ADJUST THINGS FOR THE CASE THAT THE LIQUID SUB PREDICTS DRY AEROSOL
CC
C      IF ((WATER).LE.TINY) THEN
C         DO 100 I=1,NIONS
C            MOLAL(I)= ZERO           ! Aqueous phase
C  100    CONTINUE
C         WATER   = ZERO
CC
C         CNH42S4 = CNH42SO           ! Solid phase
C         CNA2SO4 = CNA2SO
C         CNAHSO4 = CNAHSO
C         CNH4HS4 = CNH4HSO
C         CLC     = CLCO
C         CNH4NO3 = CNH4N3O
C         CNANO3  = CNANO
C         CNACL   = CNACLO                                                  
C         CNH4CL  = CNH4CLO 
CC
C         GNH3    = GNH3O             ! Gas phase
C         GHNO3   = GHNO3O
C         GHCL    = GHCLO
CC
C         GOTO 200
C      ENDIF
CC
CC *** FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
CC
C      DAMSUL  = CNH42SO - CNH42S4
C      DSOSUL  = CNA2SO  - CNA2SO4
C      DAMBIS  = CNH4HSO - CNH4HS4
C      DSOBIS  = CNAHSO  - CNAHSO4
C      DLC     = CLCO    - CLC
C      DAMNIT  = CNH4N3O - CNH4NO3
C      DAMCHL  = CNH4CLO - CNH4CL
C      DSONIT  = CNANO   - CNANO3
C      DSOCHL  = CNACLO  - CNACL
CC
CC *** FIND GAS DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
CC
C      DAMG    = GNH3O   - GNH3 
C      DHAG    = GHCLO   - GHCL
C      DNAG    = GHNO3O  - GHNO3
CC
CC *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
CC
CC     LIQUID
CC
C      MOLAL(1)= ONEMWF*MOLAL(1)                                 ! H+
C      MOLAL(2)= ONEMWF*(2.D0*DSOSUL + DSOBIS + DSONIT + DSOCHL) ! NA+
C      MOLAL(3)= ONEMWF*(2.D0*DAMSUL + DAMG   + DAMBIS + DAMCHL +
C     &                  3.D0*DLC    + DAMNIT )                  ! NH4+
C      MOLAL(4)= ONEMWF*(     DAMCHL + DSOCHL + DHAG)            ! CL-
C      MOLAL(5)= ONEMWF*(     DAMSUL + DSOSUL + DLC - MOLAL(6))  ! SO4-- !VB 17 Sept 2001
C      MOLAL(6)= ONEMWF*(   MOLAL(6) + DSOBIS + DAMBIS + DLC)    ! HSO4-
C      MOLAL(7)= ONEMWF*(     DAMNIT + DSONIT + DNAG)            ! NO3-
C      WATER   = ONEMWF*WATER
CC
CC     SOLID
CC
C      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
C      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
C      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
C      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
C      CLC     = WF*CLCO    + ONEMWF*CLC
C      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
C      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
C      CNACL   = WF*CNACLO  + ONEMWF*CNACL
C      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL
CC
CC     GAS
CC
C      GNH3    = WF*GNH3O   + ONEMWF*GNH3
C      GHNO3   = WF*GHNO3O  + ONEMWF*GHNO3
C      GHCL    = WF*GHCLO   + ONEMWF*GHCL
CC
CC *** RETURN POINT
CC
C200   RETURN
CC
CC *** END OF SUBROUTINE CALCMDRH ****************************************
CC
C      END
C
C
C
C
C
C
CC=======================================================================
CC
CC *** ISORROPIA CODE
CC *** SUBROUTINE CALCMDRP
CC
CC     THIS IS THE CASE WHERE THE RELATIVE HUMIDITY IS IN THE MUTUAL
CC     DRH REGION. THE SOLUTION IS ASSUMED TO BE THE SUM OF TWO WEIGHTED
CC     SOLUTIONS ; THE 'DRY' SOLUTION (SUBROUTINE DRYCASE) AND THE
CC     'SATURATED LIQUID' SOLUTION (SUBROUTINE LIQCASE).
CC
CC *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
CC *** GEORGIA INSTITUTE OF TECHNOLOGY
CC *** WRITTEN BY ATHANASIOS NENES
CC *** UPDATED BY CHRISTOS FOUNTOUKIS
CC
CC=======================================================================
CC
C      SUBROUTINE CALCMDRP (RHI, RHDRY, RHLIQ, DRYCASE, LIQCASE)
C      INCLUDE 'isrpia_b.inc'
C      EXTERNAL DRYCASE, LIQCASE
CC
CC *** FIND WEIGHT FACTOR **********************************************
CC
C      IF (WFTYP.EQ.0) THEN
C         WF = ONE
C      ELSEIF (WFTYP.EQ.1) THEN
C         WF = 0.5D0
C      ELSE
C         WF = (RHLIQ-RHI)/(RHLIQ-RHDRY)
C      ENDIF
C      ONEMWF  = ONE - WF
CC
CC *** FIND FIRST SECTION ; DRY ONE ************************************
CC
C      CALL DRYCASE
C      IF (ABS(ONEMWF).LE.1D-5) GOTO 200  ! DRY AEROSOL
CC
C      CNH42SO = CNH42S4              ! FIRST (DRY) SOLUTION
C      CNH4HSO = CNH4HS4
C      CLCO    = CLC 
C      CNH4N3O = CNH4NO3
C      CNH4CLO = CNH4CL
C      CNA2SO  = CNA2SO4
C      CNAHSO  = CNAHSO4
C      CNANO   = CNANO3
C      CNACLO  = CNACL
CC
CC *** FIND SECOND SECTION ; DRY & LIQUID ******************************
CC
C      CNH42S4 = ZERO
C      CNH4HS4 = ZERO
C      CLC     = ZERO
C      CNH4NO3 = ZERO
C      CNH4CL  = ZERO
C      CNA2SO4 = ZERO
C      CNAHSO4 = ZERO
C      CNANO3  = ZERO
C      CNACL   = ZERO
C      GNH3    = ZERO
C      GHNO3   = ZERO
C      GHCL    = ZERO
C      CALL LIQCASE                   ! SECOND (LIQUID) SOLUTION
CC
CC *** ADJUST THINGS FOR THE CASE THAT THE LIQUID SUB PREDICTS DRY AEROSOL
CC
C      IF ((WATER).LE.TINY) THEN
C         WATER = ZERO
C         DO 100 I=1,NIONS
C            MOLAL(I)= ZERO
C 100     CONTINUE
C         CALL DRYCASE
C         GOTO 200
C      ENDIF
CC
CC *** FIND SALT DISSOLUTIONS BETWEEN DRY & LIQUID SOLUTIONS.
CC
C      DAMBIS  = CNH4HSO - CNH4HS4
C      DSOBIS  = CNAHSO  - CNAHSO4
C      DLC     = CLCO    - CLC
CC
CC *** FIND SOLUTION AT MDRH BY WEIGHTING DRY & LIQUID SOLUTIONS.
CC
CC *** SOLID
CC
C      CNH42S4 = WF*CNH42SO + ONEMWF*CNH42S4
C      CNA2SO4 = WF*CNA2SO  + ONEMWF*CNA2SO4
C      CNAHSO4 = WF*CNAHSO  + ONEMWF*CNAHSO4
C      CNH4HS4 = WF*CNH4HSO + ONEMWF*CNH4HS4
C      CLC     = WF*CLCO    + ONEMWF*CLC
C      CNH4NO3 = WF*CNH4N3O + ONEMWF*CNH4NO3
C      CNANO3  = WF*CNANO   + ONEMWF*CNANO3
C      CNACL   = WF*CNACLO  + ONEMWF*CNACL
C      CNH4CL  = WF*CNH4CLO + ONEMWF*CNH4CL
CC
CC *** LIQUID
CC
C      WATER   = ONEMWF*WATER
CC
C      MOLAL(2)= WAER(1) - 2.D0*CNA2SO4 - CNAHSO4 - CNANO3 -     
C     &                         CNACL                            ! NA+
C      MOLAL(3)= WAER(3) - 2.D0*CNH42S4 - CNH4HS4 - CNH4CL - 
C     &                    3.D0*CLC     - CNH4NO3                ! NH4+
C      MOLAL(4)= WAER(5) - CNACL - CNH4CL                        ! CL-
C      MOLAL(7)= WAER(4) - CNANO3 - CNH4NO3                      ! NO3-
C      MOLAL(6)= ONEMWF*(MOLAL(6) + DSOBIS + DAMBIS + DLC)       ! HSO4-
C      MOLAL(5)= WAER(2) - MOLAL(6) - CLC - CNH42S4 - CNA2SO4    ! SO4--
CC
C      A8      = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C      IF ((MOLAL(5)).LE.TINY) THEN
C         HIEQ = SQRT(XKW *RH*WATER*WATER)  ! Neutral solution
C      ELSE
C         HIEQ = A8*MOLAL(6)/MOLAL(5)          
C      ENDIF
C      HIEN    = MOLAL(4) + MOLAL(7) + MOLAL(6) + 2.D0*MOLAL(5) -
C     &          MOLAL(2) - MOLAL(3)
C      MOLAL(1)= MAX (HIEQ, HIEN)                                ! H+
CC
CC *** GAS (ACTIVITY COEFS FROM LIQUID SOLUTION)
CC
C      A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
C      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
C      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
CC
C      GNH3    = MOLAL(3)/MAX((MOLAL(1)),TINY)/A2
C      GHNO3   = MOLAL(1)*MOLAL(7)/A3
C      GHCL    = MOLAL(1)*MOLAL(4)/A4
CC
C200   RETURN
CC
CC *** END OF SUBROUTINE CALCMDRP ****************************************
CC
C      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCHS4
C *** THIS SUBROUTINE CALCULATES THE HSO4 GENERATED FROM (H,SO4).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCHS4 (HI, SO4I, HSO4I, DELTA)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION HI, SO4I, HSO4I, DELTA, BB, CC, DD, SQDD, DELTA1,
     &                 DELTA2
CC      CHARACTER*40 errinf
C
C *** IF TOO LITTLE WATER, DONT SOLVE
C
      IF ((WATER).LE.1d1*TINY) THEN
         DELTA = ZERO 
         RETURN
      ENDIF
C
C *** CALCULATE HSO4 SPECIATION *****************************************
C
      A8 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
      BB =-(HI + SO4I + A8)
      CC = HI*SO4I - HSO4I*A8
      DD = BB*BB - 4.D0*CC
C
      IF ((DD).GE.ZERO) THEN
         SQDD   = SQRT(DD)
         DELTA1 = 0.5*(-BB + SQDD)
         DELTA2 = 0.5*(-BB - SQDD)
         IF ((HSO4I).LE.TINY) THEN
            DELTA = DELTA2
         ELSEIF( (HI*SO4I) .GE. (A8*HSO4I) ) THEN
            DELTA = DELTA2
         ELSEIF( (HI*SO4I) .LT. (A8*HSO4I) ) THEN
            DELTA = DELTA1
         ELSE
            DELTA = ZERO
         ENDIF
      ELSE
         DELTA  = ZERO
      ENDIF
CCC
CCC *** COMPARE DELTA TO TOTAL H+ ; ESTIMATE EFFECT OF HSO4 ***************
CCC
CC      HYD = MAX(HI, MOLAL(1))
CC      IF (HYD.GT.TINY) THEN
CC         IF (DELTA/HYD.GT.0.1d0) THEN
CC            WRITE (ERRINF,'(1PE10.3)') DELTA/HYD*100.0
CC            CALL PUSHERR (0020, ERRINF)
CC         ENDIF
CC      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCHS4 *****************************************
C
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCPH
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCPH (GG, HI, OHI)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION CN, GG, HI, OHI, BB, CC, DD
C
      AKW  = XKW *RH*WATER*WATER
      CN   = SQRT(AKW)
C
C *** GG = (negative charge) - (positive charge)
C
      IF ((GG).GT.TINY) THEN                        ! H+ in excess
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = MAX(0.5D0*(-BB + SQRT(DD)),CN)
         OHI= AKW/HI
      ELSE                                        ! OH- in excess
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= MAX(0.5D0*(-BB + SQRT(DD)),CN)
         HI = AKW/OHI
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCPH ******************************************
C
      END

C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE RSTGAM
C *** RESETS ACTIVITY COEFFICIENT ARRAYS TO DEFAULT VALUE OF 0.1
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE RSTGAM
      INCLUDE 'isrpia_b.inc'
C
      DO I=1, NPAIR
         GAMA(I) = 0.1D0
      ENDDO
C
C *** END OF SUBROUTINE RSTGAM ******************************************
C
      RETURN
      END      
C
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE RSTGAMP
C *** RESETS ACTIVITY COEFFICIENT ARRAYS TO DEFAULT VALUE OF 0.1 IF 
C *** GREATER THAN THE THRESHOLD VALUE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE RSTGAMP
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION GMAX, GTHRESH
      INTEGER I
C
      GTHRESH = 100.D0
      GMAX    = 0.1D0
      DO I=1, NPAIR
         GMAX = MAX(GMAX,GAMA(I))
      ENDDO
      IF ((GMAX) .GT. (GTHRESH)) THEN
         DO I = 1,NPAIR
            GAMA(I)  = 1.D-1
            GAMIN(I) = GREAT
            GAMOU(I) = GREAT
         ENDDO
         CALAOU   = .TRUE.
         FRST     = .TRUE.
      ENDIF
C      
      END SUBROUTINE RSTGAMP
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCACT3
C *** CALCULATES MULTI-COMPONENT ACTIVITY COEFFICIENTS FROM BROMLEYS
C     METHOD FOR AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
C     THE BINARY ACTIVITY COEFFICIENTS ARE CALCULATED BY
C     KUSIK-MEISNER RELATION (SUBROUTINE KMTAB or SUBROUTINE KMFUL3).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCACT3
      INCLUDE 'isrpia_b.inc'
C
      DOUBLE PRECISION G0(6,4),ZPL,ZMI,AGAMA,SION,H,CH,F1(3),F2(4)
      DOUBLE PRECISION MPL, XIJ, YJI, CHECK
      CHARACTER(40) errinf
C
C
C      G(I,J)= (F1(I)/Z(I) + F2(J)/Z(J+3)) / (Z(I)+Z(J+3)) - H
C
C *** SAVE ACTIVITIES IN OLD ARRAY *************************************
C
      IF (FRST) THEN               ! Outer loop
         GAMOU = GAMA
      ENDIF
C
      GAMIN = GAMA
C
C *** CALCULATE IONIC ACTIVITY OF SOLUTION *****************************
C
      IONIC=0.D0
      DO I=1,7
         IONIC=IONIC + MOLAL(I)*Z(I)*Z(I)
      ENDDO
      CHECK = 0.5d0*IONIC/WATER
C      WRITE(*,*) 'Check',CHECK
      IF (CHECK .GT. 200.d0) THEN
C         WRITE(*,*) 'Threshold exceeded in CALCACT: WATER',water,'IONIc'
C     &       ,IONIC
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCACT (',CHECK,')'
         CALL PUSHERR (0102, ERRINF)    ! WARNING ERROR: EXCEED IONIC THRESHOLD
      ELSEIF (CHECK .LT. TINY) THEN
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCACT (',CHECK,')'
         CALL PUSHERR (0102, ERRINF)    ! WARNING ERROR: EXCEED IONIC THRESHOLD
      ENDIF
      IONIC = MAX(MIN(0.5D0*IONIC/WATER,200.d0), TINY)
C
C *** CALCULATE BINARY ACTIVITY COEFFICIENTS ***************************
C
C  G0(1,1)=G11;G0(1,2)=G07;G0(1,3)=G08;G0(1,4)=G10;G0(2,1)=G01;G0(2,2)=G02
C  G0(2,3)=G12;G0(2,4)=G03;G0(3,1)=G06;G0(3,2)=G04;G0(3,3)=G09;G0(3,4)=G05
C
C
      CALL KMFUL3 (IONIC, TEMP,G01,G02,G03,
     &           G04,G05,G06,G07,G08,G09,
     &           G10,G11,G12)
C
      G0(1,1)=G11
      G0(1,2)=G07
      G0(1,3)=G08
      G0(1,4)=G10
      G0(2,1)=G01
      G0(2,2)=G02
      G0(2,3)=G12
      G0(2,4)=G03
      G0(3,1)=G06
      G0(3,2)=G04
      G0(3,3)=G09
      G0(3,4)=G05
C
C *** CALCULATE MULTICOMPONENT ACTIVITY COEFFICIENTS *******************
C
      AGAMA = 0.511D0*(298.D0/TEMP)**1.5D0    ! Debye Huckel const. at T
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1.D0+SION)

C
      DO I=1,3
         F1(I)=0.D0
         F2(I)=0.D0
      ENDDO
      F2(4)=0.D0
C
      DO I=1,3
         ZPL = Z(I)
         MPL = MOLAL(I)/WATER
         DO J=1,4
            ZMI   = Z(J+3)
            CH    = 0.25D0*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*MPL
            YJI   = CH*MOLAL(J+3)/WATER
            F1(I) = F1(I) + (YJI*(G0(I,J) + ZPL*ZMI*H))
            F2(J) = F2(J) + (XIJ*(G0(I,J) + ZPL*ZMI*H))
         ENDDO
      ENDDO
C
C *** LOG10 OF ACTIVITY COEFFICIENTS ***********************************
C
C      GAMA(01) = G(2,1)*ZZ(01)                     ! NACL
      GAMA(01) = ((F1(2)/Z(2) + F2(1)/Z(4)) / (Z(2)+Z(4)) - H)*ZZ(01)  ! NACL
C      GAMA(02) = G(2,2)*ZZ(02)                     ! NA2SO4
      GAMA(02) = ((F1(2)/Z(2) + F2(2)/Z(5)) / (Z(2)+Z(5)) - H)*ZZ(02)  ! NA2SO4
C      GAMA(03) = G(2,4)*ZZ(03)                     ! NANO3
      GAMA(03) = ((F1(2)/Z(2) + F2(4)/Z(7)) / (Z(2)+Z(7)) - H)*ZZ(03)  ! NANO3
C      GAMA(04) = G(3,2)*ZZ(04)                     ! (NH4)2SO4
      GAMA(04) = ((F1(3)/Z(3) + F2(2)/Z(5)) / (Z(3)+Z(5)) - H)*ZZ(04)  ! (NH4)2SO4
C      GAMA(05) = G(3,4)*ZZ(05)                     ! NH4NO3
      GAMA(05) = ((F1(3)/Z(3) + F2(4)/Z(7)) / (Z(3)+Z(7)) - H)*ZZ(05)  ! NH4NO3
C      GAMA(06) = G(3,1)*ZZ(06)                     ! NH4CL
      GAMA(06) = ((F1(3)/Z(3) + F2(1)/Z(4)) / (Z(3)+Z(4)) - H)*ZZ(06)  ! NH4CL
C      GAMA(07) = G(1,2)*ZZ(07)                     ! 2H-SO4
      GAMA(07) = ((F1(1)/Z(1) + F2(2)/Z(5)) / (Z(1)+Z(5)) - H)*ZZ(07)  ! 2H-SO4
C      GAMA(08) = G(1,3)*ZZ(08)                     ! H-HSO4
      GAMA(08) = ((F1(1)/Z(1) + F2(3)/Z(6)) / (Z(1)+Z(6)) - H)*ZZ(08)  ! H-HSO4
C      GAMA(09) = G(3,3)*ZZ(09)                     ! NH4HSO4
      GAMA(09) = ((F1(3)/Z(3) + F2(3)/Z(6)) / (Z(3)+Z(6)) - H)*ZZ(09)  ! NH4HSO4
C      GAMA(10) = G(1,4)*ZZ(10)                     ! HNO3
      GAMA(10) = ((F1(1)/Z(1) + F2(4)/Z(7)) / (Z(1)+Z(7)) - H)*ZZ(10)  ! HNO3
C      GAMA(11) = G(1,1)*ZZ(11)                     ! HCL
      GAMA(11) = ((F1(1)/Z(1) + F2(1)/Z(4)) / (Z(1)+Z(4)) - H)*ZZ(11)  ! HCL
C      GAMA(12) = G(2,3)*ZZ(12)                     ! NAHSO4
      GAMA(12) = ((F1(2)/Z(2) + F2(3)/Z(6)) / (Z(2)+Z(6)) - H)*ZZ(12)  ! NAHSO4
      GAMA(13) = 0.2D0*(3.D0*GAMA(04)+2.D0*GAMA(09))  ! LC ; SCAPE
C
C *** CONVERT LOG (GAMA) COEFFICIENTS TO GAMA **************************
C
      DO I=1,13
         GAMA(I)=MAX(MIN(GAMA(I),5.0d0), -5.0d0) ! F77 LIBRARY ROUTINE
         GAMA(I)=10.D0**GAMA(I)
      ENDDO
C
C *** SETUP ACTIVITY CALCULATION FLAGS *********************************
C
C OUTER CALCULATION LOOP ; ONLY IF FRST=.TRUE.
C
      IF (FRST) THEN
         ERROU = ZERO                    ! CONVERGENCE CRITERION
         DO I=1,13
            ERROU=MAX(ERROU, ((GAMOU(I)-GAMA(I))/GAMOU(I)))
         ENDDO
         CALAOU = (ERROU) .GE. (EPSACT)      ! SETUP FLAGS
         FRST   =.FALSE.
      ENDIF
C
C INNER CALCULATION LOOP ; ALWAYS
C
      ERRIN = ZERO                       ! CONVERGENCE CRITERION
      DO I=1,13
         ERRIN = MAX(ERRIN, ABS((GAMIN(I)-GAMA(I))/GAMIN(I)))
      ENDDO
      CALAIN = (ERRIN) .GE. (EPSACT)
C
      ICLACT = ICLACT + 1                ! Increment ACTIVITY call counter
C
C *** END OF SUBROUTINE ACTIVITY ****************************************
C
      RETURN
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCACT3
C *** CALCULATES MULTI-COMPONENT ACTIVITY COEFFICIENTS FROM BROMLEYS
C     METHOD FOR AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
C     THE BINARY ACTIVITY COEFFICIENTS ARE CALCULATED BY
C     KUSIK-MEISNER RELATION (SUBROUTINE KMTAB or SUBROUTINE KMFUL3).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCACT3P
      INCLUDE 'isrpia_b.inc'
C
      DOUBLE PRECISION G0(6,4),ZPL,ZMI,AGAMA,SION,H,CH,F1(3),F2(4)
      DOUBLE PRECISION MPL, XIJ, YJI, CHECK
      CHARACTER(40) errinf
C
C
C      G(I,J)= (F1(I)/Z(I) + F2(J)/Z(J+3)) / (Z(I)+Z(J+3)) - H
C
C *** CALCULATE IONIC ACTIVITY OF SOLUTION *****************************
C
      IONIC=0.D0
      DO I=1,7
         IONIC=IONIC + MOLAL(I)*Z(I)*Z(I)
      ENDDO
      CHECK = 0.5d0*IONIC/WATER
      IF (CHECK .GT. 200.d0) THEN
C         WRITE(*,*) 'Threshold exceeded in CALCACT: WATER',water,'IONIc'
C     &       ,IONIC
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCACT (',CHECK,')'
         CALL PUSHERR (0102, ERRINF)    ! WARNING ERROR: EXCEED IONIC THRESHOLD
      ELSEIF (CHECK .LT. TINY) THEN
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCACT (',CHECK,')'
         CALL PUSHERR (0102, ERRINF)    ! WARNING ERROR: EXCEED IONIC THRESHOLD
      ENDIF
      IONIC = MAX(MIN(0.5D0*IONIC/WATER,200.d0), TINY)
C
C *** CALCULATE BINARY ACTIVITY COEFFICIENTS ***************************
C
C  G0(1,1)=G11;G0(1,2)=G07;G0(1,3)=G08;G0(1,4)=G10;G0(2,1)=G01;G0(2,2)=G02
C  G0(2,3)=G12;G0(2,4)=G03;G0(3,1)=G06;G0(3,2)=G04;G0(3,3)=G09;G0(3,4)=G05
C
C
      CALL KMFUL3 (IONIC, TEMP,G01,G02,G03,
     &           G04,G05,G06,G07,G08,G09,
     &           G10,G11,G12)
C
      G0(1,1)=G11
      G0(1,2)=G07
      G0(1,3)=G08
      G0(1,4)=G10
      G0(2,1)=G01
      G0(2,2)=G02
      G0(2,3)=G12
      G0(2,4)=G03
      G0(3,1)=G06
      G0(3,2)=G04
      G0(3,3)=G09
      G0(3,4)=G05
C
C *** CALCULATE MULTICOMPONENT ACTIVITY COEFFICIENTS *******************
C
      AGAMA = 0.511D0*(298.D0/TEMP)**1.5D0    ! Debye Huckel const. at T
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1.D0+SION)
C
      DO I=1,3
         F1(I)=0.D0
         F2(I)=0.D0
      ENDDO
      F2(4)=0.D0
C
      DO I=1,3
         ZPL = Z(I)
         MPL = MOLAL(I)/WATER
         DO J=1,4
            ZMI   = Z(J+3)
            CH    = 0.25D0*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*MPL
            YJI   = CH*MOLAL(J+3)/WATER
            F1(I) = F1(I) + (YJI*(G0(I,J) + ZPL*ZMI*H))
            F2(J) = F2(J) + (XIJ*(G0(I,J) + ZPL*ZMI*H))
         ENDDO
      ENDDO
C
C *** LOG10 OF ACTIVITY COEFFICIENTS ***********************************
C
C      GAMA(01) = G(2,1)*ZZ(01)                     ! NACL
      GAMA(01) = ((F1(2)/Z(2) + F2(1)/Z(4)) / (Z(2)+Z(4)) - H)*ZZ(01)  ! NACL
C      GAMA(02) = G(2,2)*ZZ(02)                     ! NA2SO4
      GAMA(02) = ((F1(2)/Z(2) + F2(2)/Z(5)) / (Z(2)+Z(5)) - H)*ZZ(02)  ! NA2SO4
C      GAMA(03) = G(2,4)*ZZ(03)                     ! NANO3
      GAMA(03) = ((F1(2)/Z(2) + F2(4)/Z(7)) / (Z(2)+Z(7)) - H)*ZZ(03)  ! NANO3
C      GAMA(04) = G(3,2)*ZZ(04)                     ! (NH4)2SO4
      GAMA(04) = ((F1(3)/Z(3) + F2(2)/Z(5)) / (Z(3)+Z(5)) - H)*ZZ(04)  ! (NH4)2SO4
C      GAMA(05) = G(3,4)*ZZ(05)                     ! NH4NO3
      GAMA(05) = ((F1(3)/Z(3) + F2(4)/Z(7)) / (Z(3)+Z(7)) - H)*ZZ(05)  ! NH4NO3
C      GAMA(06) = G(3,1)*ZZ(06)                     ! NH4CL
      GAMA(06) = ((F1(3)/Z(3) + F2(1)/Z(4)) / (Z(3)+Z(4)) - H)*ZZ(06)  ! NH4CL
C      GAMA(07) = G(1,2)*ZZ(07)                     ! 2H-SO4
      GAMA(07) = ((F1(1)/Z(1) + F2(2)/Z(5)) / (Z(1)+Z(5)) - H)*ZZ(07)  ! 2H-SO4
C      GAMA(08) = G(1,3)*ZZ(08)                     ! H-HSO4
      GAMA(08) = ((F1(1)/Z(1) + F2(3)/Z(6)) / (Z(1)+Z(6)) - H)*ZZ(08)  ! H-HSO4
C      GAMA(09) = G(3,3)*ZZ(09)                     ! NH4HSO4
      GAMA(09) = ((F1(3)/Z(3) + F2(3)/Z(6)) / (Z(3)+Z(6)) - H)*ZZ(09)  ! NH4HSO4
C      GAMA(10) = G(1,4)*ZZ(10)                     ! HNO3
      GAMA(10) = ((F1(1)/Z(1) + F2(4)/Z(7)) / (Z(1)+Z(7)) - H)*ZZ(10)  ! HNO3
C      GAMA(11) = G(1,1)*ZZ(11)                     ! HCL
      GAMA(11) = ((F1(1)/Z(1) + F2(1)/Z(4)) / (Z(1)+Z(4)) - H)*ZZ(11)  ! HCL
C      GAMA(12) = G(2,3)*ZZ(12)                     ! NAHSO4
      GAMA(12) = ((F1(2)/Z(2) + F2(3)/Z(6)) / (Z(2)+Z(6)) - H)*ZZ(12)  ! NAHSO4
      GAMA(13) = 0.2D0*(3.D0*GAMA(04)+2.D0*GAMA(09))  ! LC ; SCAPE
C
C *** CONVERT LOG (GAMA) COEFFICIENTS TO GAMA **************************
C
      DO I=1,13
         GAMA(I)=MAX(MIN(GAMA(I),5.0d0), -5.0d0) ! F77 LIBRARY ROUTINE
         GAMA(I)=10.D0**GAMA(I)
      ENDDO
C
      ICLACT = ICLACT + 1                ! Increment ACTIVITY call counter
C
C *** END OF SUBROUTINE ACTIVITY ****************************************
C
      RETURN
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCACT3F
C *** CALCULATES MULTI-COMPONENT ACTIVITY COEFFICIENTS FROM BROMLEYS
C     METHOD FOR AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
C     THE BINARY ACTIVITY COEFFICIENTS ARE CALCULATED BY
C     KUSIK-MEISNER RELATION (SUBROUTINE KMTAB or SUBROUTINE KMFUL3).
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCACT3F
      INCLUDE 'isrpia_b.inc'
C
      DOUBLE PRECISION G0(6,4),ZPL,ZMI,AGAMA,SION,H,CH,F1(3),F2(4)
      DOUBLE PRECISION MPL, XIJ, YJI, CHECK
      CHARACTER(40) errinf
C
C      G(I,J)= (F1(I)/Z(I) + F2(J)/Z(J+3)) / (Z(I)+Z(J+3)) - H
C
C *** CALCULATE IONIC ACTIVITY OF SOLUTION *****************************
C
      IONIC=0.D0
      DO I=1,7
         IONIC=IONIC + MOLAL(I)*Z(I)*Z(I)
      ENDDO
      CHECK = 0.5d0*IONIC/WATER
      IF (CHECK .GT. 200.d0) THEN
C         WRITE(*,*) 'Threshold exceeded in CALCACT: WATER',water,'IONIc'
C     &       ,IONIC
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCACT (',CHECK,')'
C  *  slc.debug         
C         WRITE(*,*) '102,ACT3F,',CHECK
C         WRITE(*,*) 'Water: ',WATER ,', IONIC: ',IONIC 
C         WRITE(*,*) 'W: ',W
C         WRITE(*,*) 'RH: ',RH, ', TEMP:',TEMP
         CALL PUSHERR (0102, ERRINF)    ! WARNING ERROR: EXCEED IONIC THRESHOLD
      ELSEIF (CHECK .LT. TINY) THEN
         WRITE(ERRINF, '(A,E12.5,A)') 'CALCACT (',CHECK,')'
C  *  slc.debug         
C         WRITE(*,*) '102,ACT3F,',CHECK
C         WRITE(*,*) 'Water: ',WATER ,', IONIC: ',IONIC 
C         WRITE(*,*) 'W: ',W
C         WRITE(*,*) 'RH: ',RH, ', TEMP:',TEMP
         CALL PUSHERR (0102, ERRINF)    ! WARNING ERROR: EXCEED IONIC THRESHOLD
      ENDIF
      IONIC = MAX(MIN(0.5D0*IONIC/WATER,200.d0), TINY)
C
C *** CALCULATE BINARY ACTIVITY COEFFICIENTS ***************************
C
C  G0(1,1)=G11;G0(1,2)=G07;G0(1,3)=G08;G0(1,4)=G10;G0(2,1)=G01;G0(2,2)=G02
C  G0(2,3)=G12;G0(2,4)=G03;G0(3,1)=G06;G0(3,2)=G04;G0(3,3)=G09;G0(3,4)=G05
C
C
      CALL KMFUL3 (IONIC, TEMP,G01,G02,G03,
     &           G04,G05,G06,G07,G08,G09,
     &           G10,G11,G12)
C
      G0(1,1)=G11
      G0(1,2)=G07
      G0(1,3)=G08
      G0(1,4)=G10
      G0(2,1)=G01
      G0(2,2)=G02
      G0(2,3)=G12
      G0(2,4)=G03
      G0(3,1)=G06
      G0(3,2)=G04
      G0(3,3)=G09
      G0(3,4)=G05
C
C *** CALCULATE MULTICOMPONENT ACTIVITY COEFFICIENTS *******************
C
      AGAMA = 0.511D0*(298.D0/TEMP)**1.5D0    ! Debye Huckel const. at T
      SION  = SQRT(IONIC)
      H     = AGAMA*SION/(1.D0+SION)
C
      DO I=1,3
         F1(I)=0.D0
         F2(I)=0.D0
      ENDDO
      F2(4)=0.D0
C
      DO I=1,3
         ZPL = Z(I)
         MPL = MOLAL(I)/WATER
         DO J=1,4
            ZMI   = Z(J+3)
            CH    = 0.25D0*(ZPL+ZMI)*(ZPL+ZMI)/IONIC
            XIJ   = CH*MPL
            YJI   = CH*MOLAL(J+3)/WATER
            F1(I) = F1(I) + (YJI*(G0(I,J) + ZPL*ZMI*H))
            F2(J) = F2(J) + (XIJ*(G0(I,J) + ZPL*ZMI*H))
         ENDDO
      ENDDO
C
C *** LOG10 OF ACTIVITY COEFFICIENTS ***********************************
C
C      GAMA(01) = G(2,1)*ZZ(01)                     ! NACL
      GAMA(01) = ((F1(2)/Z(2) + F2(1)/Z(4)) / (Z(2)+Z(4)) - H)*ZZ(01)  ! NACL
C      GAMA(02) = G(2,2)*ZZ(02)                     ! NA2SO4
      GAMA(02) = ((F1(2)/Z(2) + F2(2)/Z(5)) / (Z(2)+Z(5)) - H)*ZZ(02)  ! NA2SO4
C      GAMA(03) = G(2,4)*ZZ(03)                     ! NANO3
      GAMA(03) = ((F1(2)/Z(2) + F2(4)/Z(7)) / (Z(2)+Z(7)) - H)*ZZ(03)  ! NANO3
C      GAMA(04) = G(3,2)*ZZ(04)                     ! (NH4)2SO4
      GAMA(04) = ((F1(3)/Z(3) + F2(2)/Z(5)) / (Z(3)+Z(5)) - H)*ZZ(04)  ! (NH4)2SO4
C      GAMA(05) = G(3,4)*ZZ(05)                     ! NH4NO3
      GAMA(05) = ((F1(3)/Z(3) + F2(4)/Z(7)) / (Z(3)+Z(7)) - H)*ZZ(05)  ! NH4NO3
C      GAMA(06) = G(3,1)*ZZ(06)                     ! NH4CL
      GAMA(06) = ((F1(3)/Z(3) + F2(1)/Z(4)) / (Z(3)+Z(4)) - H)*ZZ(06)  ! NH4CL
C      GAMA(07) = G(1,2)*ZZ(07)                     ! 2H-SO4
      GAMA(07) = ((F1(1)/Z(1) + F2(2)/Z(5)) / (Z(1)+Z(5)) - H)*ZZ(07)  ! 2H-SO4
C      GAMA(08) = G(1,3)*ZZ(08)                     ! H-HSO4
      GAMA(08) = ((F1(1)/Z(1) + F2(3)/Z(6)) / (Z(1)+Z(6)) - H)*ZZ(08)  ! H-HSO4
C      GAMA(09) = G(3,3)*ZZ(09)                     ! NH4HSO4
      GAMA(09) = ((F1(3)/Z(3) + F2(3)/Z(6)) / (Z(3)+Z(6)) - H)*ZZ(09)  ! NH4HSO4
C      GAMA(10) = G(1,4)*ZZ(10)                     ! HNO3
      GAMA(10) = ((F1(1)/Z(1) + F2(4)/Z(7)) / (Z(1)+Z(7)) - H)*ZZ(10)  ! HNO3
C      GAMA(11) = G(1,1)*ZZ(11)                     ! HCL
      GAMA(11) = ((F1(1)/Z(1) + F2(1)/Z(4)) / (Z(1)+Z(4)) - H)*ZZ(11)  ! HCL
C      GAMA(12) = G(2,3)*ZZ(12)                     ! NAHSO4
      GAMA(12) = ((F1(2)/Z(2) + F2(3)/Z(6)) / (Z(2)+Z(6)) - H)*ZZ(12)  ! NAHSO4
      GAMA(13) = 0.2D0*(3.D0*GAMA(04)+2.D0*GAMA(09))  ! LC ; SCAPE
C
C *** CONVERT LOG (GAMA) COEFFICIENTS TO GAMA **************************
C
      DO I=1,13
         GAMA(I)=MAX(MIN(GAMA(I),5.0d0), -5.0d0) ! F77 LIBRARY ROUTINE
         GAMA(I)=10.D0**GAMA(I)
      ENDDO
C
      ICLACT = ICLACT + 1                ! Increment ACTIVITY call counter
C
C *** END OF SUBROUTINE ACTIVITY ****************************************
C
      RETURN
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE KMFUL3
C *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD
C     FOR AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE KMFUL3 (IONIC,TEMP,G01,G02,G03,G04,G05,G06,G07,G08,G09,
     &                  G10,G11,G12)
      IMPLICIT NONE
      DOUBLE PRECISION IONIC, TEMP, SION, TI, TC, CF1, CF2
      DOUBLE PRECISION G01,G02,G03,G04,G05,G06,G07,G08,G09,
     &                  G10,G11,G12
      DOUBLE PRECISION Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q10, Q11
      DOUBLE PRECISION Z01,Z02,Z03,Z04,Z05,Z06,Z07,Z08,Z10,Z11
      DATA Z01,Z02,Z03,Z04,Z05,Z06,Z07,Z08,Z10,Z11
     &    /1.d0, 2.d0, 1.d0, 2.d0, 1.d0, 1.d0, 2.d0, 1.d0, 1.d0, 1.d0/
C
      SION = SQRT(IONIC)
C
C *** Coefficients at 25 oC
C
      Q1  = 2.230D0
      Q2  = -0.19D0
      Q3  = -0.39D0
      Q4  = -0.25D0
      Q5  = -1.15D0
      Q6  = 0.820D0
      Q7  = -.100D0
      Q8  = 8.000D0
      Q10 = 2.600D0
      Q11 = 6.000D0
C
      CALL MKBI(Q1 , IONIC, SION, Z01, G01)
      CALL MKBI(Q2 , IONIC, SION, Z02, G02)
      CALL MKBI(Q3 , IONIC, SION, Z03, G03)
      CALL MKBI(Q4 , IONIC, SION, Z04, G04)
      CALL MKBI(Q5 , IONIC, SION, Z05, G05)
      CALL MKBI(Q6 , IONIC, SION, Z06, G06)
      CALL MKBI(Q7 , IONIC, SION, Z07, G07)
      CALL MKBI(Q8 , IONIC, SION, Z08, G08)
      CALL MKBI(Q10, IONIC, SION, Z10, G10)
      CALL MKBI(Q11, IONIC, SION, Z11, G11)
C
C *** Correct for T other than 298 K
C
      TI  = TEMP-273.D0
      TC  = TI-25.D0
      IF (ABS(TC) .GT. 1.D0) THEN
         CF1 = 1.125D0-0.005D0*TI
         CF2 = (0.125D0-0.005D0*TI)*(0.039D0*IONIC**0.92D0-
     &         0.41D0*SION/(1.D0+SION))
         G01 = CF1*G01 - CF2*Z01
         G02 = CF1*G02 - CF2*Z02
         G03 = CF1*G03 - CF2*Z03
         G04 = CF1*G04 - CF2*Z04
         G05 = CF1*G05 - CF2*Z05
         G06 = CF1*G06 - CF2*Z06
         G07 = CF1*G07 - CF2*Z07
         G08 = CF1*G08 - CF2*Z08
         G10 = CF1*G10 - CF2*Z10
         G11 = CF1*G11 - CF2*Z11
      ENDIF
C
      G09 = G06 + G08 - G11
      G12 = G01 + G08 - G11
C
C *** Return point ; End of subroutine
C
      RETURN
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE MKBI
C *** CALCULATES BINARY ACTIVITY COEFFICIENTS BY KUSIK-MEISSNER METHOD. 
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE MKBI(Q,IONIC,SION,ZIP,BI)
C
      IMPLICIT NONE
      DOUBLE PRECISION Q, IONIC, SION, ZIP, BI
      DOUBLE PRECISION B, C, XX
C
      B=.75D0-.065D0*Q
C      C= 1.0
C      IF (IONIC.LT.6.0) C=1.+.055*Q*EXP(-.023*IONIC*IONIC*IONIC)
      C=1.+.055D0*Q*EXP(-.023D0*IONIC*IONIC*IONIC)
      XX=-0.5107D0*SION/(1.D0+C*SION)
      BI=(1.D0+B*(1.D0+.1D0*IONIC)**Q-B)
      BI=ZIP*LOG10(BI) + ZIP*XX
C
      RETURN
      END
C
CC*************************************************************************
CC
CC  TOOLBOX LIBRARY v.1.0 (May 1995)
CC
CC  Program unit   : SUBROUTINE CHRBLN
CC  Purpose        : Position of last non-blank character in a string
CC  Author         : Athanasios Nenes
CC
CC  ======================= ARGUMENTS / USAGE =============================
CC
CC  STR        is the CHARACTER variable containing the string examined
CC  IBLK       is a INTEGER variable containing the position of last non
CC             blank character. If string is all spaces (ie '   '), then
CC             the value returned is 1.
CC
CC  EXAMPLE:
CC             STR = 'TEST1.DAT     '
CC             CALL CHRBLN (STR, IBLK)
CC
CC  after execution of this code segment, "IBLK" has the value "9", which
CC  is the position of the last non-blank character of "STR".
CC
CC***********************************************************************
CC
      SUBROUTINE CHRBLN (STR, IBLK)
CC
CC***********************************************************************
      CHARACTER*(*) STR
C
      IBLK = 1                       ! Substring pointer (default=1)
      ILEN = LEN(STR)                ! Length of string
      DO 10 i=ILEN,1,-1
         IF (STR(i:i).NE.' ' .AND. STR(i:i).NE.CHAR(0)) THEN
            IBLK = i
            RETURN
         ENDIF
10    CONTINUE
      RETURN
C
      END


CC*************************************************************************
CC
CC  TOOLBOX LIBRARY v.1.0 (May 1995)
CC
CC  Program unit   : SUBROUTINE SHFTRGHT
CC  Purpose        : RIGHT-JUSTIFICATION FUNCTION ON A STRING
CC  Author         : Athanasios Nenes
CC
CC  ======================= ARGUMENTS / USAGE =============================
CC
CC  STRING     is the CHARACTER variable with the string to be justified
CC
CC  EXAMPLE:
CC             STRING    = 'AAAA    '
CC             CALL SHFTRGHT (STRING)
CC          
CC  after execution of this code segment, STRING contains the value
CC  '    AAAA'.
CC
CC*************************************************************************
CC
      SUBROUTINE SHFTRGHT (CHR)
CC
CC***********************************************************************
      CHARACTER CHR*(*)
C
      I1  = LEN(CHR)             ! Total length of string
      CALL CHRBLN(CHR,I2)        ! Position of last non-blank character
      IF (I2.EQ.I1) RETURN
C
      DO 10 I=I2,1,-1            ! Shift characters
         CHR(I1+I-I2:I1+I-I2) = CHR(I:I)
         CHR(I:I) = ' '
10    CONTINUE
      RETURN
C
      END




CC*************************************************************************
CC
CC  TOOLBOX LIBRARY v.1.0 (May 1995)
CC
CC  Program unit   : SUBROUTINE RPLSTR
CC  Purpose        : REPLACE CHARACTERS OCCURING IN A STRING
CC  Author         : Athanasios Nenes
CC
CC  ======================= ARGUMENTS / USAGE =============================
CC
CC  STRING     is the CHARACTER variable with the string to be edited
CC  OLD        is the old character which is to be replaced
CC  NEW        is the new character which OLD is to be replaced with
CC  IERR       is 0 if everything went well, is 1 if 'NEW' contains 'OLD'.
CC             In this case, this is invalid, and no change is done.
CC
CC  EXAMPLE:
CC             STRING    = 'AAAA'
CC             OLD       = 'A'
CC             NEW       = 'B' 
CC             CALL RPLSTR (STRING, OLD, NEW)
CC          
CC  after execution of this code segment, STRING contains the value
CC  'BBBB'.
CC
CC*************************************************************************
CC
      SUBROUTINE RPLSTR (STRING, OLD, NEW, IERR)
CC
CC***********************************************************************
      CHARACTER STRING*(*), OLD*(*), NEW*(*)
C
C *** INITIALIZE ********************************************************
C
      ILO = LEN(OLD)
C
C *** CHECK AND SEE IF 'NEW' CONTAINS 'OLD', WHICH CANNOT ***************
C      
      IP = INDEX(NEW,OLD)
      IF (IP.NE.0) THEN
         IERR = 1
         RETURN
      ELSE
         IERR = 0
      ENDIF
C
C *** PROCEED WITH REPLACING *******************************************
C      
10    IP = INDEX(STRING,OLD)      ! SEE IF 'OLD' EXISTS IN 'STRING'
      IF (IP.EQ.0) RETURN         ! 'OLD' DOES NOT EXIST ; RETURN
      STRING(IP:IP+ILO-1) = NEW   ! REPLACE SUBSTRING 'OLD' WITH 'NEW'
      GOTO 10                     ! GO FOR NEW OCCURANCE OF 'OLD'
C
      END
        

CC*************************************************************************
CC
CC  TOOLBOX LIBRARY v.1.0 (May 1995)
CC
CC  Program unit   : SUBROUTINE INPTD
CC  Purpose        : Prompts user for a value (DOUBLE). A default value
CC                   is provided, so if user presses <Enter>, the default
CC                   is used. 
CC  Author         : Athanasios Nenes
CC
CC  ======================= ARGUMENTS / USAGE =============================
CC
CC  VAR        is the DOUBLE PRECISION variable which value is to be saved 
CC  DEF        is a DOUBLE PRECISION variable, with the default value of VAR.        
CC  PROMPT     is a CHARACTER varible containing the prompt string.     
CC  PRFMT      is a CHARACTER variable containing the FORMAT specifier
CC             for the default value DEF.
CC  IERR       is an INTEGER error flag, and has the values:
CC             0 - No error detected.
CC             1 - Invalid FORMAT and/or Invalid default value.
CC             2 - Bad value specified by user
CC
CC  EXAMPLE:
CC             CALL INPTD (VAR, 1.0D0, 'Give value for A ', '*', Ierr)
CC          
CC  after execution of this code segment, the user is prompted for the
CC  value of variable VAR. If <Enter> is pressed (ie no value is specified)
CC  then 1.0 is assigned to VAR. The default value is displayed in free-
CC  format. The error status is specified by variable Ierr
CC
CC***********************************************************************
CC
      SUBROUTINE INPTD (VAR, DEF, PROMPT, PRFMT, IERR)
CC
CC***********************************************************************
      CHARACTER PROMPT*(*), PRFMT*(*), BUFFER*128
      DOUBLE PRECISION DEF, VAR
      INTEGER IERR
C
      IERR = 0
C
C *** WRITE DEFAULT VALUE TO WORK BUFFER *******************************
C
      WRITE (BUFFER, FMT=PRFMT, ERR=10) DEF
      CALL CHRBLN (BUFFER, IEND)
C
C *** PROMPT USER FOR INPUT AND READ IT ********************************
C
C      WRITE (*,*) PROMPT,' [',BUFFER(1:IEND),']: '
C      READ  (*, '(A)', ERR=20, END=20) BUFFER
      CALL CHRBLN (BUFFER,IEND)
C
C *** READ DATA OR SET DEFAULT ? ****************************************
C
      IF (IEND.EQ.1 .AND. BUFFER(1:1).EQ.' ') THEN
         VAR = DEF
      ELSE
         READ (BUFFER, *, ERR=20, END=20) VAR
      ENDIF
C
C *** RETURN POINT ******************************************************
C
30    RETURN
C
C *** ERROR HANDLER *****************************************************
C
10    IERR = 1       ! Bad FORMAT and/or bad default value
      GOTO 30
C
20    IERR = 2       ! Bad number given by user
      GOTO 30
C
      END


CC*************************************************************************
CC
CC  TOOLBOX LIBRARY v.1.0 (May 1995)
CC
CC  Program unit   : SUBROUTINE Pushend 
CC  Purpose        : Positions the pointer of a sequential file at its end
CC                   Simulates the ACCESS='APPEND' clause of a F77L OPEN
CC                   statement with Standard Fortran commands.
CC
CC  ======================= ARGUMENTS / USAGE =============================
CC
CC  Iunit      is a INTEGER variable, the file unit which the file is 
CC             connected to.
CC
CC  EXAMPLE:
CC             CALL PUSHEND (10)
CC          
CC  after execution of this code segment, the pointer of unit 10 is 
CC  pushed to its end.
CC
CC***********************************************************************
CC
      SUBROUTINE Pushend (Iunit)
CC
CC***********************************************************************
C
      LOGICAL OPNED
C
C *** INQUIRE IF Iunit CONNECTED TO FILE ********************************
C
      INQUIRE (UNIT=Iunit, OPENED=OPNED)
      IF (.NOT.OPNED) GOTO 25
C
C *** Iunit CONNECTED, PUSH POINTER TO END ******************************
C
10    READ (Iunit,'()', ERR=20, END=20)
      GOTO 10
C
C *** RETURN POINT ******************************************************
C
20    BACKSPACE (Iunit)
25    RETURN
      END



CC*************************************************************************
CC
CC  TOOLBOX LIBRARY v.1.0 (May 1995)
CC
CC  Program unit   : SUBROUTINE APPENDEXT
CC  Purpose        : Fix extension in file name string
CC
CC  ======================= ARGUMENTS / USAGE =============================
CC
CC  Filename   is the CHARACTER variable with the file name
CC  Defext     is the CHARACTER variable with extension (including '.',
CC             ex. '.DAT')
CC  Overwrite  is a LOGICAL value, .TRUE. overwrites any existing extension
CC             in "Filename" with "Defext", .FALSE. puts "Defext" only if 
CC             there is no extension in "Filename".
CC
CC  EXAMPLE:
CC             FILENAME1 = 'TEST.DAT'
CC             FILENAME2 = 'TEST.DAT'
CC             CALL APPENDEXT (FILENAME1, '.TXT', .FALSE.)
CC             CALL APPENDEXT (FILENAME2, '.TXT', .TRUE. )
CC          
CC  after execution of this code segment, "FILENAME1" has the value 
CC  'TEST.DAT', while "FILENAME2" has the value 'TEST.TXT'
CC
CC***********************************************************************
CC
      SUBROUTINE Appendext (Filename, Defext, Overwrite)
CC
CC***********************************************************************
      CHARACTER*(*) Filename, Defext
      LOGICAL       Overwrite
C
      CALL CHRBLN (Filename, Iend)
      IF (Filename(1:1).EQ.' ' .AND. Iend.EQ.1) RETURN  ! Filename empty
      Idot = INDEX (Filename, '.')                      ! Append extension ?
      IF (Idot.EQ.0) Filename = Filename(1:Iend)//Defext
      IF (Overwrite .AND. Idot.NE.0)
     &              Filename = Filename(:Idot-1)//Defext
      RETURN
      END



C	  SUBROUTINE TEST_QTCRT(degree,a,z)
C	  !-----------------------------------------------------------------------
C      !     Test program written to be compatible with ELF90 by
C      !        Alan Miller
C      !        amiller @ bigpond.net.au
C      !     WWW-page: http://users.bigpond.net.au/amiller
C      !     Latest revision - 27 February 1997
C      !-----------------------------------------------------------------------
C      USE constants_NSWC
C      IMPLICIT NONE
C      
C      INTEGER      :: degree, i
C      REAL (dp)    :: a(0:4)
C      COMPLEX (dp) :: z(4)
C      
C      INTERFACE
C        SUBROUTINE qdcrt (a, z)
C            USE constants_NSWC
C            IMPLICIT NONE
C            REAL (dp), INTENT(IN)     :: a(:)
C            COMPLEX (dp), INTENT(OUT) :: z(:)
C        END SUBROUTINE qdcrt
C      
C        SUBROUTINE cbcrt (a, z)
C            USE constants_NSWC
C            IMPLICIT NONE
C            REAL (dp), INTENT(IN)     :: a(:)
C            COMPLEX (dp), INTENT(OUT) :: z(:)
C        END SUBROUTINE cbcrt
C      
C        SUBROUTINE qtcrt (a, z)
C            USE constants_NSWC
C            IMPLICIT NONE
C            REAL (dp), INTENT(IN)     :: a(:)
C            COMPLEX (dp), INTENT(OUT) :: z(:)
C        END SUBROUTINE qtcrt
C      END INTERFACE
C      
C      WRITE(*, *)  'Solve quadratic, cubic, quartic eq. w/REAL coeffs'
C      WRITE(*, *)
C      
CC      DO
C        WRITE(*, *)'Enter 2, 3, 4 for quadratic, cubic or quartic eqn.:'
CC        READ(*, *) degree
C        SELECT CASE (degree)
C            CASE (2)
C              WRITE(*, *)'Enter a(0), a(1) then a(2): '
C              WRITE(*, *) a(0), a(1), a(2)
C              CALL qdcrt(a, z)
C              WRITE(*, '(a, 2(/2g20.12))') ' Rts: REAL PART  IMAG PART',
C     &                 (DBLE(z(i)), AIMAG(z(i)), i=1,2)
C            CASE (3)
C              WRITE(*, *)'Enter a(0), a(1), a(2) then a(3): '
C              WRITE(*, *) a(0), a(1), a(2), a(3)
C              CALL cbcrt(a, z)
C              WRITE(*, '(a, 3(/2g20.12))') ' Rts: REAL PART  IMAG PART',  
C     &              (DBLE(z(i)), AIMAG(z(i)), i=1,3)
C            CASE (4)
C              WRITE(*, *)'Enter a(0), a(1), a(2), a(3) then a(4): '
C              WRITE(*, *) a(0), a(1), a(2), a(3), a(4)
C              CALL qtcrt(a, z)
C              WRITE(*, '(a, 4(/2g20.12))') ' Rts: REAL PART  IMAG PART',  
C     &              (DBLE(z(i)), AIMAG(z(i)), i=1,4)
C            CASE DEFAULT
C              WRITE(*, *)'*** Try again! ***'
C              WRITE(*, *)'Use Ctrl-C to exit the program'
C        END SELECT
CC      END DO
C      
C      RETURN
C      END SUBROUTINE TEST_QTCRT
      

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE POLY3
C *** FINDS THE REAL ROOTS OF THE THIRD ORDER ALGEBRAIC EQUATION:
C     X**3 + A1*X**2 + A2*X + A3 = 0.0
C     THE EQUATION IS SOLVED ANALYTICALLY.
C
C     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
C     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS 
C     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
C     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
C
C     SOLUTION FORMULA IS FOUND IN PAGE 32 OF:
C     MATHEMATICAL HANDBOOK OF FORMULAS AND TABLES
C     SCHAUM'S OUTLINE SERIES
C     MURRAY SPIEGER, McGRAW-HILL, NEW YORK, 1968
C     (GREEK TRANSLATION: BY SOTIRIOS PERSIDES, ESPI, ATHENS, 1976)
C
C     A SPECIAL CASE IS CONSIDERED SEPERATELY ; WHEN A3 = 0, THEN
C     ONE ROOT IS X=0.0, AND THE OTHER TWO FROM THE SOLUTION OF THE
C     QUADRATIC EQUATION X**2 + A1*X + A2 = 0.0
C     THIS SPECIAL CASE IS CONSIDERED BECAUSE THE ANALYTICAL FORMULA 
C     DOES NOT YIELD ACCURATE RESULTS (DUE TO NUMERICAL ROUNDOFF ERRORS)
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE POLY3 (A1, A2, A3, ROOT, ISLV)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (EXPON=1.D0/3.D0,     ZERO=0.D0, THET1=120.D0/180.D0, 
     &           THET2=240.D0/180.D0, PI=3.14159265358932, EPS=1.D-50)
C      DOUBLE PRECISION X(3)
      DOUBLE PRECISION  X(3), A1, A2, A3, ROOT 
C
C *** SPECIAL CASE : QUADRATIC*X EQUATION *****************************
C
      IF (ABS(A3).LE.EPS) THEN 
         ISLV = 1
         IX   = 1
         X(1) = ZERO
         D    = A1*A1-4.D0*A2
         IF ((D).GE.ZERO) THEN
            IX   = 3
            SQD  = SQRT(D)
            X(2) = 0.5*(-A1+SQD)
            X(3) = 0.5*(-A1-SQD)
         ELSE
C               WRITE(*,*) 'No solution being determined'
            PAUSE
         ENDIF
      ELSE
C
C *** NORMAL CASE : CUBIC EQUATION ************************************
C
C DEFINE PARAMETERS Q, U, S, T, D 
C
         ISLV= 1
         Q   = (3.D0*A2 - A1*A1)/9.D0
         U   = (9.D0*A1*A2 - 27.D0*A3 - 2.D0*A1*A1*A1)/54.D0
         D   = Q*Q*Q + U*U
C
C *** CALCULATE ROOTS *************************************************
C
C  D < 0, THREE REAL ROOTS
C
         IF ((D).LT.-EPS) THEN        ! D < -EPS  : D < ZERO
            IX   = 3
            THET = EXPON*ACOS(U/SQRT(-Q*Q*Q))
            COEF = 2.D0*SQRT(-Q)
            X(1) = COEF*COS(THET)            - EXPON*A1
            X(2) = COEF*COS(THET + THET1*PI) - EXPON*A1
            X(3) = COEF*COS(THET + THET2*PI) - EXPON*A1
C
C  D = 0, THREE REAL (ONE DOUBLE) ROOTS
C
         ELSE IF ((D).LE.EPS) THEN    ! -EPS <= D <= EPS  : D = ZERO
            IX   = 2
            SSIG = SIGN (1.D0, U)
            S    = SSIG*(ABS(U))**EXPON
            X(1) = 2.D0*S  - EXPON*A1
            X(2) =     -S  - EXPON*A1
C
C  D > 0, ONE REAL ROOT
C
         ELSE                       ! D > EPS  : D > ZERO
            IX   = 1
            SQD  = SQRT(D)
            SSIG = SIGN (1.D0, U+SQD)       ! TRANSFER SIGN TO SSIG
            TSIG = SIGN (1.D0, U-SQD)
            S    = SSIG*(ABS(U+SQD))**EXPON ! EXPONENTIATE ABS() 
            T    = TSIG*(ABS(U-SQD))**EXPON
            X(1) = S + T - EXPON*A1
         ENDIF
      ENDIF
C
C *** SELECT APPROPRIATE ROOT *****************************************
C
      ROOT = 1.D30
      DO I=1,IX
         IF ((X(I)).GT.ZERO) THEN
            ROOT = MIN(ROOT, X(I))
            ISLV = 0
         ENDIF
      ENDDO
C
C *** END OF SUBROUTINE POLY3 *****************************************
C
      RETURN
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE POLY3B
C *** FINDS A REAL ROOT OF THE THIRD ORDER ALGEBRAIC EQUATION:
C     X**3 + A1*X**2 + A2*X + A3 = 0.0
C     THE EQUATION IS SOLVED NUMERICALLY (BISECTION).
C
C     PARAMETERS A1, A2, A3 ARE SPECIFIED BY THE USER. THE MINIMUM
C     NONEGATIVE ROOT IS RETURNED IN VARIABLE 'ROOT'. IF NO ROOT IS 
C     FOUND (WHICH IS GREATER THAN ZERO), ROOT HAS THE VALUE 1D30.
C     AND THE FLAG ISLV HAS A VALUE GREATER THAN ZERO.
C
C     RTLW, RTHI DEFINE THE INTERVAL WHICH THE ROOT IS LOOKED FOR.
C
C=======================================================================
C
      SUBROUTINE POLY3B (A1, A2, A3, RTLW, RTHI, ROOT, ISLV)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ZERO=0.D0, EPS=1D-15, MAXIT=100, NDIV=5)
C
      FUNC(X) = X**3.d0 + A1*X**2.0 + A2*X + A3
C
C *** INITIAL VALUES FOR BISECTION *************************************
C
      X1   = RTLW
      Y1   = FUNC(X1)
      IF (ABS(Y1).LE.EPS) THEN     ! Is low a root?
         ROOT = RTLW
         GOTO 50
      ENDIF
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO ***********************
C
      DX = (RTHI-RTLW)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNC (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2) .LT. ZERO) GOTO 20 ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 is a root
         ROOT = X2
      ELSE
         ROOT = 1.d30
         ISLV = 1
      ENDIF
      GOTO 50
C
C *** BISECTION *******************************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         Y3 = FUNC (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
C
C *** CONVERGED ; RETURN ***********************************************
C
40    X3   = 0.5*(X1+X2)
      Y3   = FUNC (X3)
      ROOT = X3
      ISLV = 0
C
50    RETURN
C
C *** END OF SUBROUTINE POLY3B *****************************************
C
      END

ccc      PROGRAM DRIVER
ccc      DOUBLE PRECISION ROOT
cccC
ccc      CALL POLY3 (-1.d0, 1.d0, -1.d0, ROOT, ISLV)
ccc      IF (ISLV.NE.0) STOP 'Error in POLY3'
ccc      WRITE (*,*) 'Root=', ROOT
cccC
ccc      CALL POLY3B (-1.d0, 1.d0, -1.d0, -10.d0, 10.d0, ROOT, ISLV)
ccc      IF (ISLV.NE.0) STOP 'Error in POLY3B'
ccc      WRITE (*,*) 'Root=', ROOT
cccC
ccc      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION EX10
C *** 10^X FUNCTION ; ALTERNATE OF LIBRARY ROUTINE ; USED BECAUSE IT IS
C     MUCH FASTER BUT WITHOUT GREAT LOSS IN ACCURACY. , 
C     MAXIMUM ERROR IS 2%, EXECUTION TIME IS 42% OF THE LIBRARY ROUTINE 
C     (ON A 80286/80287 MACHINE, using Lahey FORTRAN 77 v.3.0).
C
C     EXPONENT RANGE IS BETWEEN -K AND K (K IS THE REAL ARGUMENT 'K')
C     MAX VALUE FOR K: 9.999
C     IF X < -K, X IS SET TO -K, IF X > K, X IS SET TO K
C
C     THE EXPONENT IS CALCULATED BY THE PRODUCT ADEC*AINT, WHERE ADEC
C     IS THE MANTISSA AND AINT IS THE MAGNITUDE (EXPONENT). BOTH 
C     MANTISSA AND MAGNITUDE ARE PRE-CALCULATED AND STORED IN LOOKUP
C     TABLES ; THIS LEADS TO THE INCREASED SPEED.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      FUNCTION EX10(X,K)
      REAL    X, EX10, Y, AINT10, ADEC10, K
      INTEGER K1, K2
      COMMON /EXPNC/ AINT10(20), ADEC10(200)
C
C *** LIMIT X TO [-K, K] RANGE *****************************************
C
      Y    = MAX(-K, MIN(X,K))   ! MIN: -9.999, MAX: 9.999
C
C *** GET INTEGER AND DECIMAL PART *************************************
C
      K1   = INT(Y)
      K2   = INT(100*(Y-K1))
C
C *** CALCULATE EXP FUNCTION *******************************************
C
      EX10 = AINT10(K1+10)*ADEC10(K2+100)
C
C *** END OF EXP FUNCTION **********************************************
C
      RETURN
      END


C=======================================================================
C
C *** ISORROPIA CODE
C *** BLOCK DATA EXPON
C *** CONTAINS DATA FOR EXPONENT ARRAYS NEEDED IN FUNCTION EXP10
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      BLOCK DATA EXPONb
C
C *** Common block definition
C
      REAL AINT10, ADEC10
      COMMON /EXPNC/ AINT10(20), ADEC10(200)
C
C *** Integer part        
C
      DATA AINT10/
     & 0.1000E-08, 0.1000E-07, 0.1000E-06, 0.1000E-05, 0.1000E-04,
     & 0.1000E-03, 0.1000E-02, 0.1000E-01, 0.1000E+00, 0.1000E+01,
     & 0.1000E+02, 0.1000E+03, 0.1000E+04, 0.1000E+05, 0.1000E+06,
     & 0.1000E+07, 0.1000E+08, 0.1000E+09, 0.1000E+10, 0.1000E+11
     & /
C
C *** decimal part        
C
      DATA (ADEC10(I),I=1,200)/
     & 0.1023E+00, 0.1047E+00, 0.1072E+00, 0.1096E+00, 0.1122E+00,
     & 0.1148E+00, 0.1175E+00, 0.1202E+00, 0.1230E+00, 0.1259E+00,
     & 0.1288E+00, 0.1318E+00, 0.1349E+00, 0.1380E+00, 0.1413E+00,
     & 0.1445E+00, 0.1479E+00, 0.1514E+00, 0.1549E+00, 0.1585E+00,
     & 0.1622E+00, 0.1660E+00, 0.1698E+00, 0.1738E+00, 0.1778E+00,
     & 0.1820E+00, 0.1862E+00, 0.1905E+00, 0.1950E+00, 0.1995E+00,
     & 0.2042E+00, 0.2089E+00, 0.2138E+00, 0.2188E+00, 0.2239E+00,
     & 0.2291E+00, 0.2344E+00, 0.2399E+00, 0.2455E+00, 0.2512E+00,
     & 0.2570E+00, 0.2630E+00, 0.2692E+00, 0.2754E+00, 0.2818E+00,
     & 0.2884E+00, 0.2951E+00, 0.3020E+00, 0.3090E+00, 0.3162E+00,
     & 0.3236E+00, 0.3311E+00, 0.3388E+00, 0.3467E+00, 0.3548E+00,
     & 0.3631E+00, 0.3715E+00, 0.3802E+00, 0.3890E+00, 0.3981E+00,
     & 0.4074E+00, 0.4169E+00, 0.4266E+00, 0.4365E+00, 0.4467E+00,
     & 0.4571E+00, 0.4677E+00, 0.4786E+00, 0.4898E+00, 0.5012E+00,
     & 0.5129E+00, 0.5248E+00, 0.5370E+00, 0.5495E+00, 0.5623E+00,
     & 0.5754E+00, 0.5888E+00, 0.6026E+00, 0.6166E+00, 0.6310E+00,
     & 0.6457E+00, 0.6607E+00, 0.6761E+00, 0.6918E+00, 0.7079E+00,
     & 0.7244E+00, 0.7413E+00, 0.7586E+00, 0.7762E+00, 0.7943E+00,
     & 0.8128E+00, 0.8318E+00, 0.8511E+00, 0.8710E+00, 0.8913E+00,
     & 0.9120E+00, 0.9333E+00, 0.9550E+00, 0.9772E+00, 0.1000E+01,
     & 0.1023E+01, 0.1047E+01, 0.1072E+01, 0.1096E+01, 0.1122E+01,
     & 0.1148E+01, 0.1175E+01, 0.1202E+01, 0.1230E+01, 0.1259E+01,
     & 0.1288E+01, 0.1318E+01, 0.1349E+01, 0.1380E+01, 0.1413E+01,
     & 0.1445E+01, 0.1479E+01, 0.1514E+01, 0.1549E+01, 0.1585E+01,
     & 0.1622E+01, 0.1660E+01, 0.1698E+01, 0.1738E+01, 0.1778E+01,
     & 0.1820E+01, 0.1862E+01, 0.1905E+01, 0.1950E+01, 0.1995E+01,
     & 0.2042E+01, 0.2089E+01, 0.2138E+01, 0.2188E+01, 0.2239E+01,
     & 0.2291E+01, 0.2344E+01, 0.2399E+01, 0.2455E+01, 0.2512E+01,
     & 0.2570E+01, 0.2630E+01, 0.2692E+01, 0.2754E+01, 0.2818E+01,
     & 0.2884E+01, 0.2951E+01, 0.3020E+01, 0.3090E+01, 0.3162E+01,
     & 0.3236E+01, 0.3311E+01, 0.3388E+01, 0.3467E+01, 0.3548E+01,
     & 0.3631E+01, 0.3715E+01, 0.3802E+01, 0.3890E+01, 0.3981E+01,
     & 0.4074E+01, 0.4169E+01, 0.4266E+01, 0.4365E+01, 0.4467E+01,
     & 0.4571E+01, 0.4677E+01, 0.4786E+01, 0.4898E+01, 0.5012E+01,
     & 0.5129E+01, 0.5248E+01, 0.5370E+01, 0.5495E+01, 0.5623E+01,
     & 0.5754E+01, 0.5888E+01, 0.6026E+01, 0.6166E+01, 0.6310E+01,
     & 0.6457E+01, 0.6607E+01, 0.6761E+01, 0.6918E+01, 0.7079E+01,
     & 0.7244E+01, 0.7413E+01, 0.7586E+01, 0.7762E+01, 0.7943E+01,
     & 0.8128E+01, 0.8318E+01, 0.8511E+01, 0.8710E+01, 0.8913E+01,
     & 0.9120E+01, 0.9333E+01, 0.9550E+01, 0.9772E+01, 0.1000E+02
     & /
C
C *** END OF BLOCK DATA EXPON ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE PUSHERR
C *** THIS SUBROUTINE SAVES AN ERROR MESSAGE IN THE ERROR STACK
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE PUSHERR (IERR,ERRINF)
      INCLUDE 'isrpia_b.inc'
      CHARACTER(LEN=*) ERRINF
C
C *** SAVE ERROR CODE IF THERE IS ANY SPACE ***************************
C
C      WRITE(*,*) 'Calling Error, IERR: ',IERR,ERRINF
C      PAUSE
      IF (NOFER.LT.NERRMX) THEN   
         NOFER         = NOFER + 1 
         ERRSTK(NOFER) = IERR
         ERRMSG(NOFER) = ERRINF
         STKOFL        =.FALSE.
      ELSE
         STKOFL        =.TRUE.      ! STACK OVERFLOW
      ENDIF
C
C *** END OF SUBROUTINE PUSHERR ****************************************
C
      END
      


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISERRINF
C *** THIS SUBROUTINE OBTAINS A COPY OF THE ERROR STACK (& MESSAGES) 
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISERRINF (ERRSTKI, ERRMSGI, NOFERI, STKOFLI)
      INCLUDE 'isrpia_b.inc'
      CHARACTER(40) ERRMSGI(NERRMX)
      INTEGER   ERRSTKI
      LOGICAL   STKOFLI
      DIMENSION ERRSTKI(NERRMX)
C
C *** OBTAIN WHOLE ERROR STACK ****************************************
C
      DO I=1,NOFER              ! Error messages & codes
        ERRSTKI(I) = ERRSTK(I)
        ERRMSGI(I) = ERRMSG(I)
      ENDDO
C
      STKOFLI = STKOFL
      NOFERI  = NOFER
C
C      WRITE(*,*) 'NOFER', NOFER
      RETURN
C
C *** END OF SUBROUTINE ISERRINF ***************************************
C
      END
      


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ERRSTAT
C *** THIS SUBROUTINE REPORTS ERROR MESSAGES TO UNIT 'IO'
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ERRSTAT (IO,IERR,ERRINF)
      INCLUDE 'isrpia_b.inc'
      CHARACTER CER*4, NCIS*29, NCIF*27, NSIS*26, NSIF*24
      CHARACTER*40 ERRINF(1)
      DATA NCIS /'NO CONVERGENCE IN SUBROUTINE '/,
     &     NCIF /'NO CONVERGENCE IN FUNCTION '  /,
     &     NSIS /'NO SOLUTION IN SUBROUTINE '   /,
     &     NSIF /'NO SOLUTION IN FUNCTION '     /
C
C *** WRITE ERROR IN CHARACTER *****************************************
C
C      WRITE (*,'(I4)') IERR
      WRITE (CER,'(I4)') IERR
      CALL RPLSTR (CER, ' ', '0',IOK)   ! REPLACE BLANKS WITH ZEROS
      CALL CHRBLN (ERRINF, IEND)        ! LAST POSITION OF ERRINF CHAR
C
C *** WRITE ERROR TYPE (FATAL, WARNING ) *******************************
C
      IF (IERR.EQ.0) THEN
         WRITE (IO,1000) 'NO ERRORS DETECTED '
         GOTO 10
C
      ELSE IF (IERR.LT.0) THEN
         WRITE (IO,1000) 'ERROR STACK EXHAUSTED '
         GOTO 10
C
      ELSE IF (IERR.GT.1000) THEN
         WRITE (IO,1100) 'FATAL',CER
C
      ELSE
         WRITE (IO,1100) 'WARNING',CER
      ENDIF
C
C *** WRITE ERROR MESSAGE **********************************************
C
C FATAL MESSAGES
C
      IF (IERR.EQ.1001) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE (IO,1000) 'CASE NOT SUPPORTED IN CALCMR ['//SCASE(1:IEND)
     &                   //']'
C
      ELSEIF (IERR.EQ.1002) THEN 
         CALL CHRBLN (SCASE, IEND)
         WRITE (IO,1000) 'CASE NOT SUPPORTED ['//SCASE(1:IEND)//']'
C
C WARNING MESSAGES
C
      ELSEIF (IERR.EQ.0001) THEN 
         WRITE (IO,1000) NSIS,ERRINF
C
      ELSEIF (IERR.EQ.0002) THEN 
         WRITE (IO,1000) NCIS,ERRINF
C
      ELSEIF (IERR.EQ.0003) THEN 
         WRITE (IO,1000) NSIF,ERRINF
C
      ELSEIF (IERR.EQ.0004) THEN 
         WRITE (IO,1000) NCIF,ERRINF
C
      ELSE IF (IERR.EQ.0019) THEN
         WRITE (IO,1000) 'HNO3(aq) AFFECTS H+, WHICH '//
     &                   'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE (IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'
C
      ELSE IF (IERR.EQ.0020) THEN
         IF ((W(4)).GT.TINY .AND. (W(5)).GT.TINY) THEN
            WRITE (IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT HNO3,'
     &                    //'HCL DISSOLUTION'
         ELSE
            WRITE (IO,1000) 'HSO4-SO4 EQUILIBRIUM MIGHT AFFECT NH3 '
     &                    //'DISSOLUTION'
         ENDIF
         WRITE (IO,1000) 'DIRECT DECREASE IN H+ [',ERRINF(1:IEND),'] %'
C
      ELSE IF (IERR.EQ.0021) THEN
         WRITE (IO,1000) 'HNO3(aq),HCL(aq) AFFECT H+, WHICH '//
     &                   'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE (IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'
C
      ELSE IF (IERR.EQ.0022) THEN
         WRITE (IO,1000) 'HCL(g) EQUILIBRIUM YIELDS NONPHYSICAL '//
     &                   'DISSOLUTION'
         WRITE (IO,1000) 'A TINY AMOUNT [',ERRINF(1:IEND),'] IS '//
     &                   'ASSUMED TO BE DISSOLVED'
C
      ELSEIF (IERR.EQ.0033) THEN
         WRITE (IO,1000) 'HCL(aq) AFFECTS H+, WHICH '//
     &                   'MIGHT AFFECT SO4/HSO4 RATIO'
         WRITE (IO,1000) 'DIRECT INCREASE IN H+ [',ERRINF(1:IEND),'] %'
C
      ELSEIF (IERR.EQ.0050) THEN
         WRITE (IO,1000) 'TOO MUCH SODIUM GIVEN AS INPUT.'
         WRITE (IO,1000) 'REDUCED TO COMPLETELY NEUTRALIZE SO4,Cl,NO3.'
         WRITE (IO,1000) 'EXCESS SODIUM IS IGNORED.'
C
      ELSEIF (IERR.EQ.0100) THEN
C         WRITE(*,*) 'Executing PUSHERR 100'
         WRITE (IO,1000) 'CONVERGENCE TO VALUE OTHER THAN 0 '
         WRITE (IO,1000) 'FUNCTION AND VALUE: ',ERRINF(1:IEND),'.'
C
      ELSEIF (IERR.EQ.0101) THEN
C         WRITE(*,*) 'Executing PUSHERR 101'
         WRITE (IO,1000) 'CONVERGENCE AT INITIAL VALUE.'
         WRITE (IO,1000) 'FUNCTION AND VALUE: ',ERRINF(1:IEND),'.'
C
      ELSEIF (IERR.EQ.0102) THEN
C         WRITE(*,*) 'Executing PUSHERR 102'
         WRITE (IO,1000) 'EXCEEDED THE THRESHOLD VALUE FOR IONIC.'
         WRITE (IO,1000) 'FUNCTION AND VALUE: ',ERRINF(1:IEND),'.'
C
      ELSEIF (IERR.EQ.0103) THEN
C         WRITE(*,*) 'Executing PUSHERR 103'
         WRITE (IO,1000) 'VERY SMALL VALUE FOR TEST VARIABLE.'
         WRITE (IO,1000) 'FUNCTION AND VALUE: ',ERRINF(1:IEND),'.'
C
      ELSEIF (IERR.EQ.0104) THEN
C         WRITE(*,*) 'Executing PUSHERR 104'
         WRITE (IO,1000) 'NEWTON METHOD NOT CONVERGING.'
         WRITE (IO,1000) 'FUNCTION AND Y1-Y2 DIFF: ',ERRINF(1:IEND),'.'
C
      ELSE
         WRITE (IO,1000) 'NO DIAGNOSTIC MESSAGE AVAILABLE'
      ENDIF
C
10    RETURN
C
C *** FORMAT STATEMENTS *************************************
C
1000  FORMAT (1X,A:A:A:A:A)
1100  FORMAT (1X,A,' ERROR [',A4,']:')
C
C *** END OF SUBROUTINE ERRSTAT *****************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISORINF
C *** THIS SUBROUTINE PROVIDES INFORMATION ABOUT ISORROPIA
C
C ======================== ARGUMENTS / USAGE ===========================
C
C  OUTPUT:
C  1. [VERSI]
C     CHARACTER*15 variable. 
C     Contains version-date information of ISORROPIA 
C
C  2. [NCMP]
C     INTEGER variable. 
C     The number of components needed in input array WI
C     (or, the number of major species accounted for by ISORROPIA)
C
C  3. [NION]
C     INTEGER variable
C     The number of ions considered in the aqueous phase
C
C  4. [NAQGAS]
C     INTEGER variable
C     The number of undissociated species found in aqueous aerosol
C     phase
C
C  5. [NSOL]
C     INTEGER variable
C     The number of solids considered in the solid aerosol phase
C
C  6. [NERR]
C     INTEGER variable
C     The size of the error stack (maximum number of errors that can
C     be stored before the stack exhausts).
C
C  7. [TIN]
C     DOUBLE PRECISION variable
C     The value used for a very small number.
C
C  8. [GRT]
C     DOUBLE PRECISION variable
C     The value used for a very large number.
C 
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISORINF (VERSI, NCMP, NION, NAQGAS, NSOL, NERR, TIN,
     &                    GRT)
      INCLUDE 'isrpia_b.inc'
      DOUBLE PRECISION TIN, GRT
      CHARACTER VERSI*(*)
C
C *** ASSIGN INFO *******************************************************
C
      VERSI  = VERSION
      NCMP   = NCOMP
      NION   = NIONS
      NAQGAS = NGASAQ
      NSOL   = NSLDS
      NERR   = NERRMX
      TIN    = TINY
      GRT    = GREAT
C
      RETURN
C
C *** END OF SUBROUTINE ISORINF *******************************************
C
      END
