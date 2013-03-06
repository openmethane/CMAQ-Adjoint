      SUBROUTINE CALC_H2O_ADJ( WI, WI_ADJ, RH, T, H2O_NEW, H2O_NEW_ADJ )

      IMPLICIT NONE

      INTEGER, PARAMETER :: NCMP = 5
!     Arguments:
      REAL( 8 ), INTENT( IN )  :: WI( NCMP )
      REAL( 8 ), INTENT( IN )  :: RH, T
      REAL( 8 )                :: H2O_NEW

!     Parameters:
      INTEGER, PARAMETER :: NPAIR = 13
      REAL( 8 ), PARAMETER :: SMALL = 1.0D-20
      REAL( 8 ), PARAMETER :: Mw = 0.018D0  ! molar mass H2O (kg/mol)

!     Local Variables:
      REAL( 8 ) :: FSO4, FNH4, FNA, FNO3, FCL ! "free" ion amounts
      REAL( 8 ) :: WATER           ! kg of water for new time step
      REAL( 8 ) :: X, Y
      REAL( 8 ) :: CONC ( NCMP )   ! concentration (mol/m^3)
      REAL( 8 ) :: CONCR( NPAIR )  ! concentration (mol/m^3) ion "pairs"
      REAL( 8 ) :: M0I  ( NPAIR )  ! single-solute molalities
      INTEGER :: J
      CHARACTER( 2 ) :: SC         ! sub-case for composition

!     Adjoint variables:
      REAL( 8 ), INTENT ( OUT ) :: WI_ADJ ( NCMP )
      REAL( 8 )                 :: H2O_NEW_ADJ 
      REAL( 8 ) :: FSO4_ADJ
      REAL( 8 ) :: FNH4_ADJ
      REAL( 8 ) :: FNA_ADJ
      REAL( 8 ) :: FNO3_ADJ
      REAL( 8 ) :: FCL_ADJ
      REAL( 8 ) :: WATER_ADJ
      REAL( 8 ) :: X_ADJ
      REAL( 8 ) :: Y_ADJ
      REAL( 8 ) :: CONC_ADJ ( NCMP )
      REAL( 8 ) :: CONCR_ADJ ( NPAIR )

! Adjoint variables
      REAL( 8 ) :: FNH4_CHECK, FNA_CHECK, FNA_CHECK1, WATER_CHK
      REAL( 8 ) :: FNH4_CHK, CONC_N3_CHK( NCMP ), CONCR_N3_CHK
      REAL( 8 ) :: CONC_SC_CHK( NCMP ), CONC_Q5_CHK( NCMP )
      REAL( 8 ) :: CONCR_Q5_5_CHK
      REAL( 8 ) :: CONCR_Q5_4_CHK
      REAL( 8 ) :: CONCR_Q5_2_CHK
      REAL( 8 ) :: CONC_P2_CHK( NCMP )

! Adjoint variable declaration
      WI_ADJ = 0.0
      FSO4_ADJ = 0.0
      FNH4_ADJ = 0.0
      FNA_ADJ = 0.0
      FNO3_ADJ = 0.0
      FCL_ADJ = 0.0
      WATER_ADJ = 0.0
      X_ADJ = 0.0
      Y_ADJ = 0.0
      CONC_ADJ  = 0.0
      CONCR_ADJ = 0.0

!--------------- Begin Execution -----------------

! Forward Code:

C     Return if small concentration
      IF ( WI( 1 ) + WI( 2 ) + WI( 3 ) + WI( 4 ) + WI( 5 ) .LE. SMALL ) 
     &                                                      THEN
         H2O_NEW = SMALL
         H2O_NEW_ADJ = 0.0
         RETURN
      END IF

C     Set component array (mol/m^3) for determining salts
      CONC = WI   ! array assignment

      ! Checkpointing
      CONC_SC_CHK = CONC

C     Get the sub-case to use in determining salts
      CALL GETSC ( CONC, RH, T, SC )

C     Initialize ion "pairs" (i.e., salts) used in ZSR
      CONCR = 0.0D0   ! array assignment

C     Depending on case, determine moles of salts in solution (i.e., CONCR)
C     for ZSR calculation below

      IF ( SC .EQ. 'K2' ) THEN    ! sulfate poor (NH4-SO4 system)
         CONCR( 4 ) = MIN( CONC( 2 ), 0.5D0 * CONC( 3 ) )  ! (NH4)2SO4

      ELSE IF ( SC .EQ. 'L4' .OR. SC .EQ. 'O4' ) THEN  ! sulfate rich (no acid)
         X = 2.0D0 * CONC( 2 ) - CONC( 3 )   ! 2SO4 - NH4
         Y = CONC( 3 ) - CONC( 2 )           ! NH4 - SO4
         IF ( X .LE. Y ) THEN
            CONCR( 13 ) = X                  ! (NH4)3H(SO4)2 is MIN (X,Y)
            CONCR( 4 )  = Y - X              ! (NH4)2SO4
         ELSE
            CONCR( 13 ) = Y                  ! (NH4)3H(SO4)2 is MIN (X,Y)
            CONCR( 9 )  = X - Y              ! NH4HSO4
         END IF

      ELSE IF ( SC .EQ. 'M2' .OR. SC .EQ. 'P2' ) THEN   ! sulfate rich (free acid)
         ! Checkpointing
         CONC_P2_CHK = CONC

         CONCR( 9 ) = CONC( 3 )                                    ! NH4HSO4
         CONCR( 7 ) = MAX( CONC( 2 )-CONC( 3 ), 0.0D0 )            ! H2SO4

      ELSE IF ( SC .EQ. 'N3' ) THEN            ! sulfate poor (NH4-SO4-NO3 system)

         ! Checkpointing
         CONC_N3_CHK = CONC

         CONCR( 4 ) = MIN( CONC( 2 ), 0.5D0 * CONC( 3 ) )          ! (NH4)2SO4
         ! Checkpointing
         CONCR_N3_CHK = CONCR(4)

         FNH4       = MAX( CONC( 3 ) - 2.0D0 * CONCR( 4 ), 0.0D0 ) ! available NH4
         ! Checkpointing
         FNH4_CHK = FNH4

         CONCR( 5 ) = MAX( MIN( FNH4, CONC( 4 ) ), 0.0D0 )         ! NH4NO3=MIN(NH4,NO3)

      ELSE IF ( SC .EQ. 'Q5' ) THEN    ! sulfate poor, sodium poor (NH4-SO4-NO3-Cl-Na)
         ! Checkpointing
         CONC_Q5_CHK = CONC

         CONCR( 2 ) = 0.5D0 * CONC( 1 )                            ! Na2SO4
         ! Checkpointing
         CONCR_Q5_2_CHK = CONCR(2)

         FSO4       = MAX( CONC( 2 ) - CONCR( 2 ), 0.0D0 )         ! available SO4
         CONCR( 4 ) = MAX( MIN( FSO4, 0.5D0 * CONC( 3 ) ), SMALL ) ! NH42S4=MIN(NH4,S4)
         ! Checkpointing
         CONCR_Q5_4_CHK = CONCR(4)

         FNH4       = MAX( CONC( 3 ) - 2.0D0 * CONCR( 4 ), 0.0D0 ) ! available NH4
         ! Checkpointing
         FNH4_CHK = FNH4

         CONCR( 5 ) = MIN( FNH4, CONC( 4 ) )                       ! NH4NO3=MIN(NH4,NO3)
         ! Checkpointing
         CONCR_Q5_5_CHK = CONCR(5)

         FNH4       = MAX( FNH4 - CONCR( 5 ), 0.0D0 )              ! avaialable NH4
         ! Checkpointing
         FNH4_CHECK = FNH4

         CONCR( 6 ) = MIN( FNH4, CONC( 5 ) )                       ! NH4Cl=MIN(NH4,Cl)

      ELSE IF ( SC .EQ. 'R6' ) THEN   ! sulfate poor, sodium rich (NH4-SO4-NO3-Cl-Na)
         CONCR( 2 ) = CONC( 2 )                                    ! Na2SO4
         FNA        = MAX( CONC( 1 ) - 2.0D0 * CONCR( 2 ), 0.0D0 )

         CONCR( 3 ) = MIN( FNA, CONC( 4 ) )                        ! NaNO3
         FNO3       = MAX( CONC( 4 ) - CONCR( 3 ), 0.0D0 )
         FNA        = MAX( FNA-CONCR( 3 ), 0.0D0 )

         CONCR( 1 ) = MIN( FNA, CONC( 5 ) )                        ! NaCl
         FCL        = MAX( CONC( 5 ) - CONCR( 1 ), 0.0D0 )
         FNA        = MAX( FNA-CONCR( 1 ), 0.0D0 )
         ! Checkpointing
         FNA_CHECK1 = FNA

         CONCR( 5 ) = MIN( FNO3, CONC( 3 ) )                       ! NH4NO3
         FNO3       = MAX( FNO3-CONCR( 5 ), 0.0D0 )
         FNH4       = MAX( CONC( 3 ) - CONCR( 5 ), 0.0D0 )
         ! Checkpointing
         FNA_CHECK = FNA

         CONCR( 6 ) = MIN( FCL, FNH4 )                             ! NH4Cl

      ELSE IF ( SC .EQ. 'S6' ) THEN       ! sulfate rich (no acid) (NH4-SO4-NO3-Cl-Na)
         CONCR( 2 )  = 0.5D0 * CONC( 1 )                           ! Na2SO4
         FSO4        = MAX( CONC( 2 ) - CONCR( 2 ), 0.0D0 )
         CONCR( 13 ) = MIN( CONC( 3 ) / 3.0D0, FSO4 / 2.0D0 )      ! (NH4)3H(SO4)2
         FSO4        = MAX( FSO4 - 2.0D0 * CONCR( 13 ), 0.0D0 )
         FNH4        = MAX( CONC( 3 ) - 3.0D0 * CONCR( 13 ), 0.0D0 )

         IF ( FSO4 .LE. SMALL ) THEN          ! reduce (NH4)3H(SO4)2, add (NH4)2SO4
            CONCR( 13 ) = MAX( CONCR( 13 ) - FNH4, 0.0D0 )         ! (NH4)3H(SO4)2
            CONCR( 4 )  = 2.0D0 * FNH4                             ! (NH4)2SO4
         ELSE IF ( FNH4 .LE. SMALL ) THEN     ! reduce (NH4)3H(SO4)2, add NH4HSO4
            CONCR( 9 )  = 3.0D0 * MIN( FSO4, CONCR( 13 ) )         ! NH4HSO4
            CONCR( 13 ) = MAX( CONCR( 13 ) - FSO4, 0.0D0 )
            IF ( CONCR( 2 ) .GT. SMALL ) THEN ! reduce Na2SO4, add NaHSO4
               FSO4      = MAX( FSO4 - CONCR( 9 ) / 3.0D0, 0.0D0 )
               CONCR( 12 ) = 2.0D0 * FSO4                          ! NaHSO4
               CONCR( 2 )  = MAX( CONCR( 2 ) - FSO4, 0.0D0 )       ! Na2SO4
            END IF
         END IF

      ELSE IF ( SC .EQ. 'T3' ) THEN      ! sulfate rich (free acid) (NH4-SO4-NO3-Cl-Na)
         CONCR( 9 )  = CONC( 3 )                                   ! NH4HSO4
         CONCR( 12 ) = CONC( 1 )                                   ! NAHSO4
         CONCR( 7 )  = MAX( CONC( 2 ) - CONC( 3 ) - CONC( 1 ), 0.0D0 ) ! H2SO4

      ELSE
        PRINT*, 'aero_subs.f: case not supported ',
     &          '(metastable reverse only)'
C        STOP
      END IF

C     Get single-solute molalities for ZSR calculation             
      CALL GETM0I ( M0I )
      
C     Calculate H2O with ZSR and determine delta water             
      WATER = 0.0D0 
      DO J = 1, NPAIR
         WATER = WATER + CONCR( J ) / M0I( J )                     
      END DO
      
      ! Checkpointing
      WATER_CHK = WATER
   
      WATER = MAX ( WATER, SMALL )
      H2O_NEW = WATER / Mw 

! Adjoint Code:

      !------
      ! fwd code:
      ! H2O_NEW = WATER / Mw
      ! adj code:
      WATER_ADJ = ( 1.0 / Mw ) * H2O_NEW_ADJ
      H2O_NEW_ADJ = 0.0

      ! Checkpointing
      WATER = WATER_CHK
   
      !------
      ! fwd code:
      ! WATER = MAX ( WATER, SMALL )
      ! adj code:
      IF ( WATER < SMALL ) THEN
         WATER_ADJ = 0.0
      END IF

      DO J = NPAIR, 1, -1
         !------
         ! fwd code:
         ! WATER = WATER + CONCR( J ) / M0I( J )
         ! adj code:
         CONCR_ADJ ( J ) = CONCR_ADJ ( J ) + ( 1.0 / M0I( J ) ) *
     &                     WATER_ADJ
      END DO
 
      IF ( SC .EQ. 'K2' ) THEN    ! sulfate poor (NH4-SO4 system)
         !------
         ! fwd code: 
         ! CONCR( 4 ) = MIN( CONC( 2 ), 0.5D0 * CONC( 3 ) )  ! (NH4)2SO4
         ! adj code:
         IF ( CONC( 2 ) < 0.5D0 * CONC( 3 ) ) THEN
            CONC_ADJ ( 2 ) = CONCR_ADJ ( 4 )
            CONCR_ADJ ( 4 ) = 0.0
         ELSE
            CONC_ADJ ( 3 ) = 0.5D0 * CONCR_ADJ ( 4 )
            CONCR_ADJ ( 4 ) = 0.0
         END IF
         
      ELSE IF ( SC .EQ. 'L4' .OR. SC .EQ. 'O4' ) THEN  ! sulfate rich (no acid)
         IF ( X .LE. Y ) THEN
            !------
            ! fwd code:
            ! CONCR( 13 ) = X                  ! (NH4)3H(SO4)2 is MIN (X,Y)
            ! CONCR( 4 )  = Y - X              ! (NH4)2SO4
            ! adj code:
            Y_ADJ = CONCR_ADJ ( 4 )
            X_ADJ = - CONCR_ADJ ( 4 )
            CONCR_ADJ ( 4 ) = 0.0
            X_ADJ = X_ADJ + CONCR_ADJ ( 13 )
            CONCR_ADJ ( 13 ) = 0.0
         ELSE
            !------
            ! fwd code:
            ! CONCR( 13 ) = Y                  ! (NH4)3H(SO4)2 is MIN (X,Y)
            ! CONCR( 9 )  = X - Y              ! NH4HSO4
            ! adj code:
            X_ADJ = CONCR_ADJ ( 9 ) 
            Y_ADJ = - CONCR_ADJ ( 9 )
            CONCR_ADJ ( 9 ) = 0.0
            Y_ADJ = Y_ADJ + CONCR_ADJ ( 13 )
            CONCR_ADJ ( 13 ) = 0.0
         END IF
         !------
         ! fwd code:
         ! X = 2.0D0 * CONC( 2 ) - CONC( 3 )   ! 2SO4 - NH4
         ! Y = CONC( 3 ) - CONC( 2 )           ! NH4 - SO4
         ! adj code:
         CONC_ADJ ( 3 ) = Y_ADJ
         CONC_ADJ ( 2 ) = - Y_ADJ
         Y_ADJ = 0.0
         CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) - X_ADJ
         CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + 2.0D0 * X_ADJ
         X_ADJ = 0.0
         
      ELSE IF ( SC .EQ. 'M2' .OR. SC .EQ. 'P2' ) THEN   ! sulfate rich (free acid)
         ! Checkpointing
         CONC = CONC_P2_CHK

         !------
         ! fwd code:
         ! CONCR( 9 ) = CONC( 3 )                                    ! NH4HSO4
         ! CONCR( 7 ) = MAX( CONC( 2 )-CONC( 3 ), 0.0D0 )            ! H2SO4
         ! adj code:
         IF ( CONC( 2 ) - CONC( 3 ) > 0.0D0 ) THEN
            CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + CONCR_ADJ ( 7 )
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) - CONCR_ADJ ( 7 )
            CONCR_ADJ ( 7 ) = 0.0
         ELSE
            CONCR_ADJ ( 7 ) = 0.0
         END IF
         CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + CONCR_ADJ ( 9 )
         CONCR_ADJ ( 9 ) = 0.0
         
      ELSE IF ( SC .EQ. 'N3' ) THEN            ! sulfate poor (NH4-SO4-NO3 system)

         ! Checkpointing
         FNH4 = FNH4_CHK
         CONC = CONC_N3_CHK

         !------
         ! fwd code:
         ! CONCR( 4 ) = MIN( CONC( 2 ), 0.5D0 * CONC( 3 ) )          ! (NH4)2SO4
         ! FNH4       = MAX( CONC( 3 ) - 2.0D0 * CONCR( 4 ), 0.0D0 ) ! available NH4
         ! CONCR( 5 ) = MAX( MIN( FNH4, CONC( 4 ) ), 0.0D0 )         ! NH4NO3=MIN(NH4,NO3)
         ! adj code:
         IF ( MIN ( FNH4, CONC( 4 ) ) > 0.0D0 ) THEN
            IF ( FNH4 > CONC( 4 ) ) THEN
               CONC_ADJ ( 4 ) = CONC_ADJ ( 4 ) + CONCR_ADJ ( 5 )
               CONCR_ADJ ( 5 ) = 0.0
            ELSE
               FNH4_ADJ = FNH4_ADJ + CONCR_ADJ ( 5 )
               CONCR_ADJ ( 5 ) = 0.0
            END IF
         ELSE 
            CONCR_ADJ ( 5 ) = 0.0
         END IF
         ! Checkpointing
         CONCR(4) = CONCR_N3_CHK

         IF ( CONC( 3 ) - 2.0D0 * CONCR( 4 ) > 0.0D0 ) THEN
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + FNH4_ADJ
            CONCR_ADJ ( 4 ) = CONCR_ADJ ( 4 ) - 2.0D0 * FNH4_ADJ
            FNH4_ADJ = 0.0
         ELSE
            FNH4_ADJ = 0.0
         END IF
         IF ( CONC( 2 ) < 0.5D0 * CONC( 3 ) ) THEN
            CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + CONCR_ADJ ( 4 )
            CONCR_ADJ ( 4 ) = 0.0
         ELSE
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + 0.5D0 * CONCR_ADJ ( 4 )
            CONCR_ADJ ( 4 ) = 0.0
         END IF
        
      ELSE IF ( SC .EQ. 'Q5' ) THEN    ! sulfate poor, sodium poor (NH4-SO4-NO3-Cl-Na)
         !------
         ! fwd code:
         ! CONCR( 2 ) = 0.5D0 * CONC( 1 )                            ! Na2SO4
         ! FSO4       = MAX( CONC( 2 ) - CONCR( 2 ), 0.0D0 )         ! available SO4
         ! CONCR( 4 ) = MAX( MIN( FSO4, 0.5D0 * CONC( 3 ) ), SMALL ) ! NH42S4=MIN(NH4,S4)
         ! FNH4       = MAX( CONC( 3 ) - 2.0D0 * CONCR( 4 ), 0.0D0 ) ! available NH4
         ! CONCR( 5 ) = MIN( FNH4, CONC( 4 ) )                       ! NH4NO3=MIN(NH4,NO3)
         ! FNH4       = MAX( FNH4 - CONCR( 5 ), 0.0D0 )              ! avaialable NH4
         ! CONCR( 6 ) = MIN( FNH4, CONC( 5 ) )                       ! NH4Cl=MIN(NH4,Cl)
         ! adj code:

         ! Checkpointing
         FNH4 = FNH4_CHECK
         CONC = CONC_Q5_CHK

         IF ( FNH4 < CONC( 5 ) ) THEN
            FNH4_ADJ = FNH4_ADJ + CONCR_ADJ ( 6 )
            CONCR_ADJ ( 6 ) = 0.0
         ELSE
            CONC_ADJ ( 5 ) = CONC_ADJ ( 5 ) + CONCR_ADJ ( 6 ) 
            CONCR_ADJ ( 6 ) = 0.0
         END IF
         ! Checkpointing
         CONCR(5) = CONCR_Q5_5_CHK
         FNH4 = FNH4_CHK

         IF ( FNH4 - CONCR( 5 ) > 0.0D0 ) THEN
            CONCR_ADJ ( 5 ) = CONCR_ADJ ( 5 ) - FNH4_ADJ
         ELSE
            FNH4_ADJ = 0.0
         END IF
         IF ( FNH4 < CONC( 4 ) ) THEN
            FNH4_ADJ = FNH4_ADJ + CONCR_ADJ ( 5 )
            CONCR_ADJ ( 5 ) = 0.0
         ELSE
            CONC_ADJ ( 4 ) = CONC_ADJ ( 4 ) + CONCR_ADJ ( 5 )
            CONCR_ADJ ( 5 ) = 0.0
         END IF
         ! Checkpointing
         CONCR(4) = CONCR_Q5_4_CHK

         IF ( CONC ( 3 ) - 2.0D0 * CONCR( 4 ) > 0.0D0 ) THEN
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + FNH4_ADJ
            CONCR_ADJ ( 4 ) = CONCR_ADJ ( 4 ) - 2.0D0 * FNH4_ADJ
            FNH4_ADJ = 0.0
         ELSE
            FNH4_ADJ = 0.0
         END IF
         IF ( MIN ( FSO4, 0.5D0 * CONC( 3 ) ) > SMALL ) THEN
            IF ( FSO4 < 0.5D0 * CONC( 3 ) ) THEN
               FSO4_ADJ = FSO4_ADJ + CONCR_ADJ ( 4 )
               CONCR_ADJ ( 4 ) = 0.0
            ELSE
               CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + 0.5D0 * CONCR_ADJ ( 4 )
               CONCR_ADJ ( 4 ) = 0.0
            END IF
         ELSE
            CONCR_ADJ ( 4 ) = 0.0
         END IF
         ! Checkpointing
         CONCR(2) = CONCR_Q5_2_CHK

         IF ( CONC( 2 ) - CONCR( 2 ) > 0.0D0 ) THEN
            CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + FSO4_ADJ 
            CONCR_ADJ ( 2 ) = CONCR_ADJ ( 2 ) - FSO4_ADJ
            FSO4_ADJ = 0.0
         ELSE
            FSO4_ADJ = 0.0
         END IF
         CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) + 0.5D0 * CONCR_ADJ ( 2 ) 
         CONCR_ADJ ( 2 ) = 0.0

      ELSE IF ( SC .EQ. 'R6' ) THEN   ! sulfate poor, sodium rich (NH4-SO4-NO3-Cl-Na)

         !------
         ! fwd code: 
         ! CONCR( 6 ) = MIN( FCL, FNH4 )                             ! NH4Cl
         ! adj code:
         IF ( FCL < FNH4 ) THEN
            FCL_ADJ = FCL_ADJ + CONCR_ADJ ( 6 )
            CONCR_ADJ ( 6 ) = 0.0
         ELSE
            FNH4_ADJ = FNH4_ADJ + CONCR_ADJ ( 6 )
            CONCR_ADJ ( 6 ) = 0.0
         END IF

         !------
         ! fwd code:
         ! CONCR( 5 ) = MIN( FNO3, CONC( 3 ) )                       ! NH4NO3
         ! FNO3       = MAX( FNO3-CONCR( 5 ), 0.0D0 )
         ! FNH4       = MAX( CONC( 3 ) - CONCR( 5 ), 0.0D0 )
         ! adj code:
         IF ( CONC( 3 ) - CONCR( 5 ) > 0.0D0 ) THEN
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + FNH4_ADJ
            CONCR_ADJ ( 5 ) = CONCR_ADJ ( 5 ) - FNH4_ADJ
            FNH4_ADJ = 0.0
         ELSE
            FNH4_ADJ = 0.0
         END IF
         IF ( FNO3 - CONCR( 5 ) > 0.0D0 ) THEN
            CONCR_ADJ ( 5 ) = CONCR_ADJ ( 5 ) - FNO3_ADJ 
         ELSE
            FNO3_ADJ = 0.0
         END IF
         IF ( FNO3 < CONC( 3 ) ) THEN
            FNO3_ADJ = FNO3_ADJ + CONCR_ADJ ( 5 )
            CONCR_ADJ ( 5 ) = 0.0
         ELSE
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + CONCR_ADJ ( 5 )
            CONCR_ADJ ( 5 ) = 0.0
         END IF

         ! Checkpointing
         FNA = FNA_CHECK
 
         !------
         ! fwd code:
         ! CONCR( 1 ) = MIN( FNA, CONC( 5 ) )                        ! NaCl
         ! FCL        = MAX( CONC( 5 ) - CONCR( 1 ), 0.0D0 )
         ! FNA        = MAX( FNA-CONCR( 1 ), 0.0D0 )
         ! adj code:
         IF ( FNA - CONCR( 1 ) > 0.0D0 ) THEN
            CONCR_ADJ ( 1 ) = CONCR_ADJ ( 1 ) - FNA_ADJ
         ELSE 
            FNA_ADJ = 0.0
         END IF
         IF ( CONC( 5 ) - CONCR( 1 ) > 0.0 ) THEN
            CONC_ADJ ( 5 ) = CONC_ADJ ( 5 ) + FCL_ADJ 
            CONCR_ADJ ( 1 ) = CONCR_ADJ ( 1 ) - FCL_ADJ
            FCL_ADJ = 0.0
         ELSE
            FCL_ADJ = 0.0
         END IF
         IF ( FNA < CONC( 5 ) ) THEN
            FNA_ADJ = FNA_ADJ + CONCR_ADJ ( 1 )
            CONCR_ADJ ( 1 ) = 0.0
         ELSE
            CONC_ADJ ( 5 ) = CONC_ADJ ( 5 ) + CONCR_ADJ ( 1 )
            CONCR_ADJ ( 1 ) = 0.0
         END IF   

         ! Checkpointing
         FNA = FNA_CHECK1

         !------
         ! fwd code:
         ! CONCR( 3 ) = MIN( FNA, CONC( 4 ) )                        ! NaNO3
         ! FNO3       = MAX( CONC( 4 ) - CONCR( 3 ), 0.0D0 )
         ! FNA        = MAX( FNA-CONCR( 3 ), 0.0D0 )
         ! adj code:
         IF ( FNA - CONCR( 3 ) > 0.0D0 ) THEN
            CONCR_ADJ ( 3 ) = CONCR_ADJ ( 3 ) - FNA_ADJ
         ELSE
            FNA_ADJ = 0.0
         END IF
         IF ( CONC( 4 ) - CONCR( 3 ) > 0.0D0 ) THEN
            CONC_ADJ ( 4 ) = CONC_ADJ ( 4 ) + FNO3_ADJ
            CONCR_ADJ ( 3 ) = CONCR_ADJ ( 3 ) - FNO3_ADJ
            FNO3_ADJ = 0.0
         ELSE
            FNO3_ADJ = 0.0
         END IF
         IF ( FNA < CONC( 4 ) ) THEN
            FNA_ADJ = FNA_ADJ + CONCR_ADJ ( 3 )
            CONCR_ADJ ( 3 ) = 0.0
         ELSE
            CONC_ADJ ( 4 ) = CONC_ADJ ( 4 ) + CONCR_ADJ ( 3 )
            CONCR_ADJ ( 3 ) = 0.0
         END IF 

         !------
         ! fwd code:
         ! CONCR( 2 ) = CONC( 2 )                                 ! Na2SO4
         ! FNA        = MAX( CONC( 1 ) - 2.0D0 * CONCR( 2 ), 0.0D0 )
         ! adj code:
         IF ( CONC( 1 ) - 2.0D0 * CONCR( 2 ) > 0.0D0 ) THEN
            CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) + FNA_ADJ
            CONCR_ADJ ( 2 ) = CONCR_ADJ ( 2 ) - 2.0D0 * FNA_ADJ
            FNA_ADJ = 0.0
         ELSE
            FNA_ADJ = 0.0
         END IF
         CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + CONCR_ADJ ( 2 )
         CONCR_ADJ ( 2 ) = 0.0

      ELSE IF ( SC .EQ. 'S6' ) THEN       ! sulfate rich (no acid) (NH4-SO4-NO3-Cl-Na)

         IF ( FSO4 .LE. SMALL ) THEN          ! reduce (NH4)3H(SO4)2, add (NH4)2SO4
            !------
            ! fwd code:
            ! CONCR( 13 ) = MAX( CONCR( 13 ) - FNH4, 0.0D0 )         ! (NH4)3H(SO4)2
            ! CONCR( 4 )  = 2.0D0 * FNH4                             ! (NH4)2SO4
            ! adj code:
            FNH4_ADJ = 2.0D0 * CONCR_ADJ ( 4 )
            CONCR_ADJ ( 4 ) = 0.0
            IF ( CONCR( 13 ) - FNH4 > 0.0D0 ) THEN
               FNH4_ADJ = FNH4_ADJ - CONCR_ADJ ( 13 )
            ELSE
               CONCR_ADJ ( 13 ) = 0.0
            END IF
         ELSE IF ( FNH4 .LE. SMALL ) THEN     ! reduce (NH4)3H(SO4)2, add NH4HSO4
            IF ( CONCR( 2 ) .GT. SMALL ) THEN ! reduce Na2SO4, add NaHSO4
               !------
               ! fwd code:
               ! FSO4      = MAX( FSO4 - CONCR( 9 ) / 3.0D0, 0.0D0 )
               ! CONCR( 12 ) = 2.0D0 * FSO4                          ! NaHSO4
               ! CONCR( 2 )  = MAX( CONCR( 2 ) - FSO4, 0.0D0 )       ! Na2SO4
               ! adj code:
               IF ( CONCR( 2 ) - FSO4 > 0.0D0 ) THEN
                  FSO4_ADJ = FSO4_ADJ - CONCR_ADJ ( 2 ) 
               ELSE
                  CONCR_ADJ ( 2 ) = 0.0
               END IF
               FSO4_ADJ = FSO4_ADJ + 2.0D0 * CONCR_ADJ ( 12 )
               CONCR_ADJ ( 12 ) = 0.0 
               IF ( FSO4 - CONCR( 9 )  / 3.0D0 > 0.0D0 ) THEN
                  CONCR_ADJ ( 9 ) = CONCR_ADJ ( 9 ) - ( 1.0 / 3.0D0 ) *
     &                              FSO4_ADJ
               ELSE
                  FSO4_ADJ = 0.0
               END IF
            END IF
            !------
            ! fwd code:
            ! CONCR( 9 )  = 3.0D0 * MIN( FSO4, CONCR( 13 ) )         ! NH4HSO4
            ! CONCR( 13 ) = MAX( CONCR( 13 ) - FSO4, 0.0D0 )
            ! adj code:
            IF ( CONCR( 13 ) - FSO4 > 0.0D0 ) THEN
               FSO4_ADJ = FSO4_ADJ - CONCR_ADJ ( 13 )
            ELSE
               CONCR_ADJ ( 13 ) = 0.0
            END IF 
            IF ( FSO4 < CONCR( 13 ) ) THEN
               FSO4_ADJ = FSO4_ADJ + 3.0D0 * CONCR_ADJ ( 9 )
               CONCR_ADJ ( 9 ) = 0.0
            ELSE
               CONCR_ADJ ( 13 ) = CONCR_ADJ ( 13 ) + 3.0D0 * 
     &                            CONCR_ADJ ( 9 )
               CONCR_ADJ ( 9 ) = 0.0
            END IF
         END IF
         !------
         ! fwd code:
         ! CONCR( 2 )  = 0.5D0 * CONC( 1 )                           ! Na2SO4
         ! FSO4        = MAX( CONC( 2 ) - CONCR( 2 ), 0.0D0 )
         ! CONCR( 13 ) = MIN( CONC( 3 ) / 3.0D0, FSO4 / 2.0D0 )      ! (NH4)3H(SO4)2
         ! FSO4        = MAX( FSO4 - 2.0D0 * CONCR( 13 ), 0.0D0 )
         ! FNH4        = MAX( CONC( 3 ) - 3.0D0 * CONCR( 13 ), 0.0D0 )
         ! adj code:
         IF ( CONC( 3 ) - 3.0D0 * CONCR( 13 ) > 0.0D0 ) THEN
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + FNH4_ADJ
            CONCR_ADJ ( 13 ) = CONCR_ADJ ( 13 ) - 3.0D0 * FNH4_ADJ
            FNH4_ADJ = 0.0
         ELSE
            FNH4_ADJ = 0.0
         END IF
         IF ( FSO4 - 2.0D0 * CONCR( 13 ) > 0.0D0 ) THEN
            CONCR_ADJ ( 13 ) = CONCR_ADJ ( 13 ) - 2.0D0 * FSO4_ADJ
         ELSE
            FSO4_ADJ = 0.0
         END IF
         IF ( CONC( 3 ) / 3.0D0 < FSO4 / 2.0D0 ) THEN
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + ( 1.0 / 3.0D0 ) * 
     &                       CONCR_ADJ ( 13 )
            CONCR_ADJ ( 13 ) = 0.0
         ELSE
            FSO4_ADJ = FSO4_ADJ + ( 1.0 / 2.0D0 ) * CONCR_ADJ ( 13 )
            CONCR_ADJ ( 13 ) = 0.0
         END IF
         IF ( CONC( 2 ) - CONCR( 2 ) > 0.0D0 ) THEN
            CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + FSO4_ADJ 
            CONCR_ADJ ( 2 ) = CONCR_ADJ ( 2 ) - FSO4_ADJ
            FSO4_ADJ = 0.0
         ELSE
            FSO4_ADJ = 0.0
         END IF
         CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) + 0.5D0 * CONCR_ADJ ( 2 ) 
         CONCR_ADJ ( 2 ) = 0.0

      ELSE IF ( SC .EQ. 'T3' ) THEN      ! sulfate rich (free acid) (NH4-SO4-NO3-Cl-Na)
         !------
         ! fwd code:
         ! CONCR( 9 )  = CONC( 3 )                                   ! NH4HSO4
         ! CONCR( 12 ) = CONC( 1 )                                   ! NAHSO4
         ! CONCR( 7 )  = MAX( CONC( 2 ) - CONC( 3 ) - CONC( 1 ), 0.0D0 ) ! H2SO4
         ! adj code:
         IF ( CONC( 2 ) - CONC( 3 ) - CONC( 1 ) > 0.0D0 ) THEN
            CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + CONCR_ADJ ( 7 )
            CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) - CONCR_ADJ ( 7 )
            CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) - CONCR_ADJ ( 7 )
            CONCR_ADJ ( 7 ) = 0.0
         ELSE
            CONCR_ADJ ( 7 ) = 0.0D0
         END IF
         CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) + CONCR_ADJ ( 12 )
         CONCR_ADJ ( 12 ) = 0.0
         CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + CONCR_ADJ ( 9 )
         CONCR_ADJ ( 9 ) = 0.0

      END IF

      !------
      ! fwd code:
      ! CONCR = 0.0D0
      ! adj code:
      CONCR_ADJ = 0.0D0

      ! Checkpointing
      CONC = CONC_SC_CHK

      !------
      ! fwd code:
      ! CALL GETSC ( CONC, RH, T, SC )
      ! adj code:
      CALL GETSC_ADJ ( CONC, CONC_ADJ, RH, T, SC )

      !------
      ! fwd code:  
      ! CONC = WI
      ! adj code:
      WI_ADJ = WI_ADJ + CONC_ADJ

      IF ( WI( 1 ) + WI( 2 ) + WI( 3 ) + WI( 4 ) + WI( 5 ) .LE. SMALL )
     &                                                      THEN
         !------
         ! fwd code:
         ! H2O_NEW = SMALL
         ! adj code:
         H2O_NEW_ADJ = 0.0
      END IF
      
      END SUBROUTINE CALC_H2O_ADJ

