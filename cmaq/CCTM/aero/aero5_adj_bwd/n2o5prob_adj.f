
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
      SUBROUTINE N2O5PROB_ADJ ( TEMP, RH, GPARAM, GAMMA_ADJ )

!-----------------------------------------------------------------------
! Function:
!   Adjoint of function that claculates the N2O5 heterogeneous reaction
!   probability, which is the fraction of collisions between a gaseous
!   N2O5 molecule and a particle surface that leads to nitrate
!   production.

! Revision History:
!   Mar 2011 by Matthew Turner at UC-Boulder: created for adjoint/4dvar
!-----------------------------------------------------------------------

      USE aero_data_b
      USE met_data
                     
      IMPLICIT NONE

! *** ARGUMENTS

      REAL,    INTENT( IN ) :: TEMP     ! Air temperature [ K ]
      REAL,    INTENT( IN ) :: RH       ! Fractional relative humidity
      INTEGER, INTENT( IN ) :: GPARAM   ! switch to select among
                                        !  parameterizations
      REAL,    INTENT( INOUT ) :: GAMMA_ADJ 

! *** PARAMETERS

! *** switch for alternative parameterization of LAM1 & LAM2
!     when APNDX = .TRUE. (default), Eqs A1-A2 are used for reaction
!     probability on aqueous sulfate particles.  Alternatively, set
!     APNDX = .FALSE. to use Eqs 4-5.

      LOGICAL, PARAMETER :: APNDX = .TRUE.

! *** LOCAL VARIABLES
      
! *** chemical species concentrations [ug/m3]

      REAL      ANH4      ! i+j mode ammonium
      REAL      ANO3      ! i+j mode nitrate
      REAL      ASO4      ! i+j mode sulfate

! *** variables for computing N2O5PROB when GPARAM = 2 or 3

      REAL      FRACSO4   ! aerosol mass ratio of SO4/(SO4+NO3)
      REAL      GAMMA     ! upper limit of rxn prob
      REAL      GAMMA1    ! upper limit of rxn prob
      REAL      GAMMA2    ! lower limit of rxn prob
      REAL      ALPHA     ! RH-dependent parameter to compute GAMMA1
      REAL      BETA      ! TEMP-dependent parameter to compute GAMMA1

! *** variables for default parameterization of N2O5PROB 

      LOGICAL   CRHB      ! function to determine if RH is below CRH
      LOGICAL   CRYSTAL   ! true if ambient RH < CRH, false otherwise
      LOGICAL   IRHX      ! function to determine whether RH exceeds IRH
      LOGICAL   FROZEN    ! true if ambient RH > IRH, false otherwise
      REAL      NNO3      ! particle-phase nitrate [micromoles/m3]
      REAL      NSO4      ! particle-phase sulfate [micromoles/m3]
      REAL      NNH4      ! particle-phase ammonium [micromoles/m3]
      REAL      NANI      ! particle-phase anions [micromoles/m3]
      REAL      X1        ! mole fraction of ammonium bisulfate
      REAL      X2        ! mole fraction of ammonium sulfate
      REAL      X3        ! mole fraction of ammonium nitrate
      REAL      LAM1      ! logit transformation of N2O5PROB on 
      REAL      LAM2      !   aqueous NH4HSO4 [LAM1], aqueous (NH4)2SO4
      REAL      LAM3      !   [LAM2], aqueous NH4NO3 [LAM3], and dry
      REAL      LAMD      !   sulfate-containing particles [LAMD]
      REAL      GAM1      ! reaction probability on aqueous NH4HSO4
      REAL      GAM2      !    "          "      "     "    (NH4)2SO4
      REAL      GAM3      !    "          "      "     "    NH4NO3
      REAL      GAMD      !    "          "      " dry sulfate particles
      REAL      T293,T291 ! temperature threshold variables
      REAL      RH46      ! RH threshold variable

! *** statement function for inverting the logit transformation given
!     in Eq 7 by Davis et al (2008)

      REAL      LOGITINV  ! statement function
      REAL      XX        ! dummy argument for LOGITINV
      LOGITINV( XX ) = 1.0 / ( 1.0 + EXP( -XX ) )

! *** Adjoint variables

      REAL      X1_ADJ
      REAL      X2_ADJ
      REAL      X3_ADJ
      REAL      NNH4_ADJ
      REAL      NANI_ADJ
      REAL      NNO3_ADJ
      REAL      NSO4_ADJ
      REAL      ANH4_ADJ
      REAL      ASO4_ADJ
      REAL      ANO3_ADJ
      REAL      FRACSO4_ADJ

! *** Checkpointing variables

      REAL      ANO3_CHK
      REAL      ASO4_CHK
      REAL      GAMMA1_CHK
      REAL      GAMMA2_CHK
      REAL      X3_CHK
      REAL      NNH4_CHK
      REAL      NANI_CHK
      REAL      GAMD_CHK
      REAL      GAM1_CHK
      REAL      GAM2_CHK
      REAL      GAM3_CHK
      REAL      GAM3_CHK1
      
!.......................................................................

! initialize adj variables

      X1_ADJ = 0.0
      X2_ADJ = 0.0
      X3_ADJ = 0.0
      NNH4_ADJ = 0.0
      NANI_ADJ = 0.0
      NNO3_ADJ = 0.0
      NSO4_ADJ = 0.0
      ANH4_ADJ = 0.0
      ASO4_ADJ = 0.0
      ANO3_ADJ = 0.0
      FRACSO4_ADJ = 0.0

!--------------------------------------------------------------- 

! Forward code

c *** retrieve particle-phase ammonium, nitrate, and sulfate [ug/m3]

      ANH4 = aerospc_conc( ANH4_IDX,1 ) + aerospc_conc( ANH4_IDX,2 )
      ANO3 = aerospc_conc( ANO3_IDX,1 ) + aerospc_conc( ANO3_IDX,2 )
      ASO4 = aerospc_conc( ASO4_IDX,1 ) + aerospc_conc( ASO4_IDX,2 )
      
c *** User Option: GPARAM = 1
c     Dentener and Crutzen (1993) recommended a constant value of
c     N2O5PROB = 0.1, which was used in CMAQ prior to ver4.3.  In more
c     recent literature, this value has been recognized as an upper
c     estimate of N2O5PROB so it should not be used for routine
c     simulations.  It is included here only to facilitate sensitivity
c     studies by CMAQ model users.

      IF ( GPARAM .EQ. 1 ) THEN
         GAMMA    = 0.1
         RETURN
      END IF
      
c *** User Options: GPARAM = 2 and 3
c     These options both employ Eqs 2 and 3 by Riemer et al (2003), in
c     which N2O5PROB varies according to the particle-phase sulfate and 
c     nitrate concentrations.  In both options, the NO3 effect (i.e., 
c     GAMMA1/GAMMA2) is assumed to be a factor of 10 based on Mentel et
c     al (1999) and Riemer et al (2003).
c      - When GPARAM = 2, upper limit of N2O5PROB is fixed at 0.02.
c        This was the default setting in CMAQ ver4.3 through ver4.5.1.
c      - When GPARAM = 3, upper limit of N2O5PROB is a function of 
c        ambient TEMP & RH based on the "Sulfate" equation in Table 1
c        by Evans & Jacob (2005).  This was the default setting in CMAQ
c        ver4.6.  After that release, a typographical error was found
c        in the published equation of Evans & Jacob (2005) so this code
c        has been corrected accordingly.

      IF ( GPARAM .EQ. 2 ) THEN

         GAMMA1 = 0.02

      ELSE IF ( GPARAM .EQ. 3 ) THEN

c        In this function, RH is in fractional units whereas the
c        published equation by Evans&Jacob refers to RH as a percentage.

         ALPHA = 2.79E-4
     &         + RH * ( 1.3E-2 + RH * ( -3.43E-2 + 7.52E-2 * RH ) )

c        To fix the typographical error by Evans & Jacob (2005), the
c        sign of BETA has been switched in this code.

         IF ( TEMP .LT. 282.0 ) THEN
            GAMMA1 = 3.0199517 * ALPHA   ! (10.0 ** 0.48) * ALPHA
         ELSE
            BETA  = 0.04 * ( 294.0 - TEMP )
            GAMMA1 = ALPHA * ( 10.0 ** BETA )
         END IF

      END IF

      IF ( ( GPARAM .EQ. 2 ) .OR. ( GPARAM .EQ. 3 ) ) THEN

         ! Checkpointing
         ANO3_CHK = ANO3
         ASO4_CHK = ASO4

         IF ( ANO3 .GT. 0.0 ) THEN
            FRACSO4 = ASO4 / ( ASO4 + ANO3 )
         ELSE
            FRACSO4 = 1.0
         END IF

         GAMMA2 = 0.1 * GAMMA1

         ! Checkpointing
         GAMMA1_CHK = GAMMA1
         GAMMA2_CHK = GAMMA2

         GAMMA    = GAMMA2 + FRACSO4 * ( GAMMA1 - GAMMA2 )
         RETURN

      END IF

c *** Default setting in current version of CMAQ:
c     This code implements the paramaterization given in Eq 15 by Davis 
c     et al (2008), in which N2O5PROB is a function of RH, TEMP, 
c     particle composition, and phase state.  Note: In this function, RH
c     is in fractional units whereas the published equations refer to RH
c     as a percentage.

c *** Check whether the ambient RH is below the crystallization RH for 
c     the given inorganic particle composition.

      CRYSTAL = CRHB( RH, .TRUE. )

c *** Check whether the ambient RH exceeds the RH of ice formation.

      FROZEN = IRHX( TEMP, RH )

c *** Set N2O5PROB to constant value if particles contain ice, based on
c     Eq 14 by Davis et al (2008).

      IF ( FROZEN ) THEN
         GAMMA    = 0.02                               ! Eq 14

c *** Compute mole-fractional-composition of particles based on Eq 11 by
c     Davis et al (2008).

      ELSE
         NNO3 = ANO3 / aerospc_mw( ANO3_IDX )
         NSO4 = ASO4 / aerospc_mw( ASO4_IDX )
         NNH4 = ANH4 / aerospc_mw( ANH4_IDX )
         NANI = NNO3 + NSO4

         X3 = NNO3 / NANI

         ! Checkpointing
         X3_CHK = X3
         NNH4_CHK = NNH4
         NANI_CHK = NANI

         X2 = MAX( 0.0, MIN( 1.0 - X3, NNH4/NANI - 1.0 ) )
         X1 = 1.0 - ( X2 + X3 )

c *** Compute N2O5PROB on pure NH4NO3 particles using Eqs 6 and 8 by
c     Davis et al (2008).

         LAM3 = -8.10774 + 4.902 * RH                  ! Eq 6
         GAM3 = MIN( LOGITINV( LAM3 ), 0.0154 )        ! Eq 8

c *** Compute N2O5PROB on dry particles using Eqs 9, 10, and 13 by 
c     Davis et al (2008).
      
         IF ( CRYSTAL ) THEN
            T293     = MAX( 0.0, TEMP - 293.0 )
            LAMD     = -6.13376 + 3.592 * RH           ! Eq 9
     &                - 0.19688 * T293
            GAMD     = MIN( LOGITINV( LAMD ), 0.0124 ) ! Eq 10

            ! Checkpointing
            GAMD_CHK = GAMD
            GAM3_CHK = GAM3

            GAMMA    = ( X1 + X2 ) * GAMD              ! Eq 13
     &                + X3 * MIN( GAMD, GAM3 )
         
c *** Compute N2O5PROB on aqeuous particles using Eqs A1, A2, 8, and 12
c     by Davis et al (2008).  When APNDX = .TRUE. (default), Eqs A1-A2
c     are used for reaction probability on aqueous sulfate particles.
c     Switch to .FALSE. if Eqs 4-5 are desired.  See Appendix A by
c     Davis et al. (2008) for a discussion of these options.

         ELSE
            T291 = MAX( 0.0, TEMP - 291.0 )
            IF ( APNDX ) THEN
               RH46 = MIN( 0.0, RH - 0.46 )
               LAM2  = -3.64849 + 9.553 * RH46         ! Eq A2
               LAM1  = LAM2 + 0.97579                  ! Eqs A1 & A2
     &                - 0.20427 * T291
            ELSE
               LAM1  = -4.10612 + 2.386 * RH           ! Eq 4
     &                - 0.23771 * T291
               LAM2  = LAM1 - 0.80570                  ! Eqs 4 & 5
     &                + 0.10225 * T291
            END IF
            GAM1     = MIN( LOGITINV( LAM1 ), 0.08585 )! Eq 8
            GAM2     = MIN( LOGITINV( LAM2 ), 0.053 )  ! Eq 8

            ! Checkpointing
            GAM1_CHK = GAM1
            GAM2_CHK = GAM2
            GAM3_CHK1 = GAM3

            GAMMA    = ( X1 * GAM1 )                   ! Eq 12
     &               + ( X2 * GAM2 )
     &               + ( X3 * GAM3 )
      
         END IF

      END IF
     
! Adjoint code

      IF ( GPARAM .EQ. 1 ) THEN

         GAMMA_ADJ = 0.0

      ELSE IF ( ( GPARAM .EQ. 2 ) .OR. ( GPARAM .EQ. 3 ) ) THEN

         ! checkpointing
         GAMMA1 = GAMMA1_CHK
         GAMMA2 = GAMMA2_CHK

         !------
         ! fwd code:
         ! GAMMA2 = 0.1 * GAMMA1
         ! N2O5PROB = GAMMA2 + FRACSO4 * ( GAMMA1 - GAMMA2 )
         ! adj code:
         FRACSO4_ADJ = FRACSO4_ADJ + ( GAMMA1 - GAMMA2 ) * GAMMA_ADJ
         GAMMA_ADJ = 0.0

         ! Checkpointing
         ANO3 = ANO3_CHK
         ASO4 = ASO4_CHK

         IF ( ANO3 > 0.0 ) THEN
            ! fwd code:
            ! FRACSO4 = ASO4 / ( ASO4 + ANO3 )
            ! adj code:
            ASO4_ADJ = ASO4_ADJ + ( 1 / ( ASO4 + ANO3 ) - ASO4
     &               / ( ASO4 + ANO3 ) ** 2 ) * FRACSO4_ADJ
            ANO3_ADJ = ANO3_ADJ - ( ASO4 / ( ASO4 + ANO3 ) ** 2 )
     &               * FRACSO4_ADJ
            FRACSO4_ADJ = 0.0
         ELSE
            ! fwd code:
            ! FRACSO4 = 1.0
            ! adj code:
            FRACSO4_ADJ = 0.0
         END IF

      ELSE

         IF ( FROZEN ) THEN

            GAMMA_ADJ = 0.0

         ELSE

            IF ( CRYSTAL ) THEN

               ! Checkpointing
               GAMD = GAMD_CHK
               GAM3 = GAM3_CHK

               !------
               ! fwd code:
               ! N2O5PROB = ( X1 + X2 ) * GAMD              ! Eq 13
               !      + X3 * MIN( GAMD, GAM3 )
               ! adj code:
               X1_ADJ = X1_ADJ + GAMD * GAMMA_ADJ
               X2_ADJ = X2_ADJ + GAMD * GAMMA_ADJ
               X3_ADJ = X3_ADJ + MIN( GAMD, GAM3 ) * GAMMA_ADJ
               GAMMA_ADJ = 0.0

            ELSE

               ! Checkpointing
               GAM1 = GAM1_CHK
               GAM2 = GAM2_CHK
               GAM3 = GAM3_CHK1

               !------
               ! fwd cdoe:
               ! N2O5PROB = ( X1 * GAM1 )                   ! Eq 12
               !     + ( X2 * GAM2 )
               !     + ( X3 * GAM3 )
               ! adj code:
               X1_ADJ = X1_ADJ + GAM1 * GAMMA_ADJ
               X2_ADJ = X2_ADJ + GAM2 * GAMMA_ADJ
               X3_ADJ = X3_ADJ + GAM3 * GAMMA_ADJ
               GAMMA_ADJ = 0.0

            END IF

            ! Checkpointing
            X3 = X3_CHK
            NNH4 = NNH4_CHK
            NANI = NANI_CHK

            !------
            ! fwd code:
            ! X3 = NNO3 / NANI
            ! X2 = MAX( 0.0, MIN( 1.0 - X3, NNH4/NANI - 1.0 ) )
            ! X1 = 1.0 - ( X2 + X3 )
            ! adj code:
            X2_ADJ = X2_ADJ - X1_ADJ
            X3_ADJ = X3_ADJ - X1_ADJ
            X1_ADJ = 0.0
            IF ( MIN( 1.0 - X3, NNH4 / NANI - 1.0 ) > 0 ) THEN
               IF ( 1.0 - X3 > NNH4 / NANI - 1.0 ) THEN
                  ! fwd code:
                  ! x2 = NNH4 / NANI - 1.0
                  ! adj code:
                  NNH4_ADJ = NNH4_ADJ + 1 / NANI * X2_ADJ
                  NANI_ADJ = NANI_ADJ - NNH4 / ( NANI ** 2 ) * X2_ADJ
                  X2_ADJ = 0.0
               ELSE
                  ! fwd code:
                  ! x2 = 1 - x3
                  ! adj code:
                  X3_ADJ = X3_ADJ - X2_ADJ
                  X2_ADJ = 0.0
               END IF
            ELSE
               X2_ADJ = 0.0
            END IF
            NNO3_ADJ = NNO3_ADJ + ( 1 / NANI ) * X3_ADJ
            NANI_ADJ = NANI_ADJ - ( NNO3 / ( NANI ** 2 ) ) * X3_ADJ
            X3_ADJ = 0.0

            !------
            ! fwd code:
            ! NNO3 = ANO3 / aerospc_mw( ANO3_IDX )
            ! NSO4 = ASO4 / aerospc_mw( ASO4_IDX )
            ! NNH4 = ANH4 / aerospc_mw( ANH4_IDX )
            ! NANI = NNO3 + NSO4
            ! adj code:
            NNO3_ADJ = NNO3_ADJ + NANI_ADJ
            NSO4_ADJ = NSO4_ADJ + NANI_ADJ
            NANI_ADJ = 0.0
            ANH4_ADJ = ANH4_ADJ + NNH4_ADJ / aerospc_mw( ANH4_IDX )
            NNH4_ADJ = 0.0
            ASO4_ADJ = ASO4_ADJ + NSO4_ADJ / aerospc_mw( ASO4_IDX )
            NSO4_ADJ = 0.0
            ANO3_ADJ = ANO3_ADJ + NNO3_ADJ / aerospc_mw( ANO3_IDX )
            NSO4_ADJ = 0.0

         END IF ! end if frozen block

      END IF ! end if gparam == 1 block

      !------
      ! fwd code:
      ! ANH4 = aerospc_conc( ANH4_IDX,1 ) + aerospc_conc( ANH4_IDX,2 )
      ! ANO3 = aerospc_conc( ANO3_IDX,1 ) + aerospc_conc( ANO3_IDX,2 )
      ! ASO4 = aerospc_conc( ASO4_IDX,1 ) + aerospc_conc( ASO4_IDX,2 )
      ! adj code:
      aerospc_concb( ASO4_IDX, 1 ) = aerospc_concb( ASO4_IDX, 1 ) + ASO4_ADJ
      aerospc_concb( ASO4_IDX, 2 ) = aerospc_concb( ASO4_IDX, 2 ) + ASO4_ADJ
      ASO4_ADJ = 0.0 
      aerospc_concb( ANO3_IDX, 1 ) = aerospc_concb( ANO3_IDX, 1 ) + ANO3_ADJ
      aerospc_concb( ANO3_IDX, 2 ) = aerospc_concb( ANO3_IDX, 2 ) + ANO3_ADJ
      ANO3_ADJ = 0.0 
      aerospc_concb( ANH4_IDX, 1 ) = aerospc_concb( ANH4_IDX, 1 ) + ANH4_ADJ
      aerospc_concb( ANH4_IDX, 2 ) = aerospc_concb( ANH4_IDX, 2 ) + ANH4_ADJ
      ANH4_ADJ = 0.0

      END SUBROUTINE

