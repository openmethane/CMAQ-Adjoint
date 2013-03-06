      SUBROUTINE GETSC_ADJ   ( CONC, CONC_ADJ, RH, T, SC )

      IMPLICIT NONE

!     Parameters:
      INTEGER, PARAMETER :: NCMP  = 5    ! number of aerosol components
      REAL( 8 ), PARAMETER :: SMALL = 1.0D-20

!     Arguments:
      REAL( 8 ), INTENT( INOUT )    :: CONC( NCMP )
      REAL( 8 ), INTENT( IN )       :: RH, T
      CHARACTER( 2 )                :: SC

!     Local Variables:
      REAL( 8 ) :: T0, TCF                     ! DRH(T) factor
      REAL( 8 ) :: S4RAT, S4RATW, NaRAT, SRI   ! sulfate & sodium ratios
      REAL( 8 ) :: FSO4                        ! "free" sulfate
      REAL( 8 ) :: DNACL, DNH4CL, DNANO3, DNH4NO3, DNH42S4 ! DRH values

      REAL( 8 ) :: GETASR    ! ISORROPIA function for sulfate ratio

      LOGICAL :: NH4S4, NH4S4N3, ALLSP  ! concentration regime

! Adjoint variables
   
      REAL( 8 ), INTENT ( INOUT ) :: CONC_ADJ( NCMP )
      REAL( 8 ) :: S4RAT_ADJ
      REAL( 8 ) :: S4RATW_ADJ
      REAL( 8 ) :: NaRAT_ADJ
      REAL( 8 ) :: SRI_ADJ
      REAL( 8 ) :: FSO4_ADJ
      INTEGER   :: X

! initialize adjoint variables
      S4RAT_ADJ = 0.0
      S4RATW_ADJ = 0.0
      NaRAT_ADJ = 0.0
      SRI_ADJ = 0.0
      FSO4_ADJ = 0.0

!------------------------- Begin Execution -------------------------

! Forward Code:

      NH4S4   = .FALSE.  ! NH4-SO4       
      NH4S4N3 = .FALSE.  ! NH4-SO4-NO3
      ALLSP   = .FALSE.  ! all species

C     See if any components are negligible
      IF ( CONC( 1 ) + CONC( 4 ) + CONC( 5 ) .LE. SMALL ) THEN   ! Na,Cl,NO3=0
         NH4S4 = .TRUE.
      ELSE IF ( CONC( 1 ) + CONC( 5 ) .LE. SMALL ) THEN        ! Na,Cl=0
         NH4S4N3 = .TRUE.
      ELSE
         ALLSP = .TRUE.
      END IF

      CONC( : ) = MAX( CONC( : ), SMALL )

C     Deliquescence RH calculations
      DNH42S4 = 0.7997D0
      DNH4NO3 = 0.6183D0
      IF ( INT( T ) .NE. 298 ) THEN
         T0      = 298.15D0 
         TCF     = 1.0D0 / T - 1.0D0 / T0
         DNH4NO3 = DNH4NO3 * EXP( 852.0D0 * TCF ) 
         DNH42S4 = DNH42S4 * EXP(  80.0D0 * TCF ) 
         DNH4NO3 = MIN( DNH4NO3, DNH42S4 ) ! adjust for curves crossing T<271K
      END IF

C     Find sub-case "SC"
      IF ( NH4S4 ) THEN ! NH4-S04 system
           
         IF (RH >= DNH42S4) THEN
            S4RATW = GETASR(CONC(2), RH) ! aerosol sulfate ratio
         ELSE
            S4RATW = 2.0D0                ! dry aerosol sulfate ratio
         END IF
         S4RAT  = CONC(3) / CONC(2)     ! sulfate ratio (NH4/SO4)
        
         IF (S4RATW <= S4RAT) THEN      ! sulfate poor
            SC = 'K2'
         ELSE IF (1.0D0 <= S4RAT .AND. S4RAT < S4RATW) THEN ! sulfate rich (no acid)
            SC = 'L4'
         ELSE IF (S4RAT < 1.0D0) THEN   ! sulfate rich (free acid)
            SC = 'M2'
         END IF 

      ELSE IF ( NH4S4N3 ) THEN ! NH4-SO4-NO3 system
         
         IF (RH >= DNH4NO3) THEN
            S4RATW = GETASR(CONC(2), RH)
         ELSE
            S4RATW = 2.0D0               ! dry aerosol ratio
         END IF
         S4RAT = CONC(3) / CONC(2)


10       IF (S4RATW <= S4RAT) THEN     ! sulfate poor
            SC = 'N3'
         ELSE IF (1.0D0 <= S4RAT .AND. S4RAT < S4RATW) THEN  ! sulfate rich (no acid)
            SC = 'O4'
         ELSE IF (S4RAT < 1.0D0) THEN    ! sulfate rich (free acid)
            SC = 'P2'
         END IF
      
      ELSE IF ( ALLSP )  THEN ! all species

C     Adjust DRH of NH4NO3 for low temperature
         DNACL   = 0.7528D0
         DNANO3  = 0.7379D0
         DNH4CL  = 0.7710D0 
         IF ( INT( T ) .NE. 298 ) THEN
            DNACL   = DNACL  * EXP(  25.0D0 * TCF )
            DNANO3  = DNANO3 * EXP( 304.0D0 * TCF )
            DNH4CL  = DNH4Cl * EXP( 239.0D0 * TCF ) 
            DNH4NO3 = MIN( DNH4NO3, DNH4CL, DNANO3, DNACL )
         END IF
        
        IF ( RH .GE. DNH4NO3 ) THEN
           FSO4   = CONC( 2 ) - CONC( 1 ) / 2.0D0   ! sulfate unbound by Na+
           FSO4   = MAX( FSO4, SMALL )
           SRI    = GETASR( FSO4, RH ) 
           S4RATW = ( CONC( 1 ) + FSO4 * SRI ) / CONC( 2 )
           S4RATW = MIN( S4RATW, 2.0D0 )
        ELSE
           S4RATW = 2.0D0                       ! ratio for dry aerosol
        END IF
        S4RAT = ( CONC( 1 ) + CONC( 3 ) ) / CONC( 2 )
        NaRAT = CONC( 1 ) / CONC( 2 )
        
        IF ( S4RATW .LE. S4RAT .AND. NaRAT .LT. 2.0D0 ) THEN ! sulfate poor, sodium poor
           SC = 'Q5'
        ELSE IF ( S4RAT .GE. S4RATW .AND. NaRAT .GE. 2.0D0 ) THEN ! SO4 poor, Na rich
           SC = 'R6'
        ELSE IF ( 1.0D0 .LE. S4RAT .AND. S4RAT .LT. S4RATW ) THEN ! SO4 rich, no acid
           SC = 'S6'
        ELSE IF ( S4RAT .LT. 1.0D0 ) THEN                    ! sulfate rich, free acid

           SC = 'T3'
        END IF

      END IF

! Adjoint Code:

      IF ( NH4S4 ) THEN ! NH4-S04 system

         !------
         ! fwd code:
         ! S4RAT  = CONC(3) / CONC(2)     ! sulfate ratio (NH4/SO4)
         ! adj code:
         CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + ( 1.0 / CONC( 2 ) )
     &                  * S4RAT_ADJ
         CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) - ( CONC( 3 ) / ( CONC( 2 ) 
     &                  ** 2.0 ) ) * S4RAT_ADJ
!         WRITE(*,*) 'S4RAT_ADJ = ', S4RAT_ADJ    !DEBUG, MDT
         S4RAT_ADJ = 0.0
!         WRITE(*,*) '  adj:: conc( 2 ) = ', CONC( 2 )    !debug, mdt
!         WRITE(*,*) '  adj:: conc( 3 ) = ', CONC( 3 )    !debug, mdt

         IF (RH >= DNH42S4) THEN
            !------
            ! fwd code:
            ! S4RATW = GETASR(CONC(2), RH) ! aerosol sulfate ratio
            ! adj code:
!            CALL GETASR_ADJ ( CONC( 2 ) , CONC_ADJ ( 2 ), S4RATW_ADJ, RH )
!            WRITE(*,*) ' CONC_ADJ ( 2 ) = ', CONC_ADJ ( 2 ) !debug, mdt
         ELSE
            !------
            ! fwd code:
            ! S4RATW = 2.0D0                ! dry aerosol sulfate ratio
            ! adj code:
            S4RATW_ADJ = 0.0
         END IF

      ELSE IF ( NH4S4N3 ) THEN ! NH4-SO4-NO3 system

         !------
         ! fwd code:
         ! S4RAT = CONC(3) / CONC(2)
         ! adj code:
         CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + ( 1.0 / CONC( 2 ) ) * 
     &                    S4RAT_ADJ
         CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) - ( CONC( 3 ) / ( CONC( 2 ) **
     &                    2.0 ) ) * S4RAT_ADJ
!         WRITE(*,*) ' S4RAT_ADJ = ', S4RAT_ADJ    !debug, mdt
         S4RAT_ADJ = 0.0
!         WRITE(*,*) '  adj:: CONC( 2 ) = ', CONC( 2 )   !debug, mdt
!         WRITE(*,*) '  adj:: CONC( 3 ) = ', CONC( 3 )   !debug, mdt
         IF (RH >= DNH4NO3) THEN
            !------
            ! fwd code:
            ! S4RATW = GETASR(CONC(2), RH)
            ! adj code:
!            CALL GETASR_ADJ ( CONC( 2 ), CONC_ADJ ( 2 ), S4RATW_ADJ, RH )
            S4RATW_ADJ = 0.0
         ELSE
            !------
            ! fwd code:
            ! S4RATW = 2.0D0               ! dry aerosol ratio
            ! adj code:
            S4RATW_ADJ = 0.0
         END IF

      ELSE IF ( ALLSP )  THEN ! all species
      
        !------
        ! fwd code:
        ! S4RAT = ( CONC( 1 ) + CONC( 3 ) ) / CONC( 2 )
        ! NaRAT = CONC( 1 ) / CONC( 2 )
        ! adj code:
        CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) + ( 1.0 / CONC( 2 ) ) * 
     &                   NaRAT_ADJ
        CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) - ( CONC( 2 ) / ( CONC( 2 ) **
     &                   2.0 ) ) * NaRAT_ADJ
!         WRITE(*,*) ' NaRAT_ADJ  = ', NaRAT_ADJ   !debug, mdt
        NaRAT_ADJ = 0.0
        CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) + ( 1.0 / CONC( 2 ) ) * 
     &                   S4RAT_ADJ
        CONC_ADJ ( 3 ) = CONC_ADJ ( 3 ) + ( 1.0 / CONC( 2 ) ) * 
     &                   S4RAT_ADJ
        CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) - ( ( CONC( 1 ) + CONC ( 3 ) )
     &                   / ( CONC( 2 ) **  2.0 ) ) * S4RAT_ADJ
!         WRITE(*,*) ' S4RAT_ADJ  = ', S4RAT_ADJ   !debug, mdt
!         WRITE(*,*) '  adj:: CONC( 1 ) = ', CONC( 1 )   !debug, mdt
!         WRITE(*,*) '  adj:: CONC( 2 ) = ', CONC( 2 )   !debug, mdt
!         WRITE(*,*) '  adj:: CONC( 3 ) = ', CONC( 3 )   !debug, mdt
         S4RAT_ADJ = 0.0

        IF ( RH .GE. DNH4NO3 ) THEN
           !------
           ! fwd code:
           ! FSO4   = CONC( 2 ) - CONC( 1 ) / 2.0D0   ! sulfate unbound by Na+
           ! FSO4   = MAX( FSO4, SMALL )
           ! SRI    = GETASR( FSO4, RH )
           ! S4RATW = ( CONC( 1 ) + FSO4 * SRI ) / CONC( 2 )
           ! S4RATW = MIN( S4RATW, 2.0D0 )
           ! adj code:
           IF ( S4RATW > 2.0D0 ) THEN
              S4RATW_ADJ = 0.0
           END IF
           CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) + ( 1.0 / CONC( 2 ) ) * 
     &                      S4RATW_ADJ
           FSO4_ADJ = FSO4_ADJ + SRI / CONC( 2 ) * S4RATW_ADJ
           SRI_ADJ = SRI_ADJ + FSO4 / CONC( 2 ) * S4RATW_ADJ
           CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) - ( ( CONC( 1 ) + FSO4 * SRI 
     &                    ) / ( CONC( 2 ) ** 2.0 ) ) * S4RATW_ADJ
           S4RATW_ADJ = 0.0

!           WRITE(*,*) '  adj:: FSO4 = ', FSO4    !debug, mdt
!           WRITE(*,*) '  adj:: SRI  = ', SRI     !debug, mdt
!           CALL GETASR_ADJ ( FSO4, FSO4_ADJ, SRI_ADJ, RH )
           IF ( FSO4 < SMALL ) THEN
              FSO4_ADJ = 0.0
           END IF
           CONC_ADJ ( 2 ) = CONC_ADJ ( 2 ) + FSO4_ADJ
           CONC_ADJ ( 1 ) = CONC_ADJ ( 1 ) - ( 1.0 / 2.0D0 ) * FSO4_ADJ
!           WRITE(*,*) ' FSO4_ADJ = ', FSO4_ADJ    !debug, mdt
           FSO4_ADJ = 0.0
              
        ELSE
           !------
           ! fwd code:
           ! S4RATW = 2.0D0                       ! ratio for dry aerosol
           ! adj code:
           S4RATW_ADJ = 0.0
        END IF
            
      END IF

      !------
      ! fwd code:
      ! CONC( : ) = MAX( CONC( : ), SMALL )
      ! adj code:
      DO X = 1, NCMP
         IF ( CONC( X ) < SMALL ) THEN
            CONC_ADJ ( X ) = 0.0
         END IF 
      END DO

!      WRITE(*,*) '  adj:: CONC = ', CONC    !debug, mdt

!      WRITE(*,*) ' CONC_ADJ ( 1 ) = ', CONC_ADJ ( 1 )    !debug, mdt
!      WRITE(*,*) ' CONC_ADJ ( 2 ) = ', CONC_ADJ ( 2 )    !debug, mdt
!      WRITE(*,*) ' CONC_ADJ ( 3 ) = ', CONC_ADJ ( 3 )    !debug, mdt
!      WRITE(*,*) ' CONC_ADJ ( 4 ) = ', CONC_ADJ ( 4 )    !debug, mdt
!      WRITE(*,*) ' CONC_ADJ ( 5 ) = ', CONC_ADJ ( 5 )    !debug, mdt

      END SUBROUTINE GETSC_ADJ
