#include "cppdefs.h"
      MODULE ice_bulk_flux_mod
#ifdef ICE_BULK_FLUXES
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Bulk flux formula for ice thermal dynamics.                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: ice_bulk_flux
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_bulk_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ice
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17, __LINE__, MyFile)
# endif
      CALL ice_bulk_flux_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         liold(ng), linew(ng),                    &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
# endif
# ifdef WET_DRY
     &                         GRID(ng) % rmask_wet,                    &
     &                         GRID(ng) % umask_wet,                    &
     &                         GRID(ng) % vmask_wet,                    &
# endif
# ifdef ICESHELF
     &                         GRID(ng) % zice,                         &
# endif
# ifdef WIND_MINUS_CURRENT
     &                         ICE(ng) % ui,                            &
     &                         ICE(ng) % vi,                            &
# endif
     &                         FORCES(ng) % Hair,                       &
     &                         FORCES(ng) % Pair,                       &
     &                         FORCES(ng) % Tair,                       &
     &                         FORCES(ng) % Uwind,                      &
     &                         FORCES(ng) % Vwind,                      &
# ifdef CLOUDS
     &                         FORCES(ng) % cloud,                      &
# endif
     &                         ICE(ng) % ai,                            &
     &                         ICE(ng) % tis,                           &
# ifdef ICE_ALBEDO
     &                         FORCES(ng) % albo,                       &
     &                         FORCES(ng) % albi,                       &
     &                         FORCES(ng) % frswp,                      &
# endif
     &                         FORCES(ng) % taux_ai,                    &
     &                         FORCES(ng) % tauy_ai,                    &
     &                         FORCES(ng) % qswi,                       &
     &                         FORCES(ng) % qlwi,                       &
     &                         FORCES(ng) % qswo,                       &
     &                         FORCES(ng) % qai,                        &
# ifdef ICE_DIAGS
#  ifdef ICESHELF
     &                         FORCES(ng) % qswl,                       &
     &                         FORCES(ng) % qlwl,                       &
#  endif
     &                         FORCES(ng) % qswa,                       &
     &                         FORCES(ng) % qlwa,                       &
     &                         FORCES(ng) % qlai,                       &
     &                         FORCES(ng) % qsei,                       &
# endif
     &                         ICE(ng) % coef_ice_heat,                 &
     &                         ICE(ng) % rhs_ice_heat)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE ice_bulk_flux
!
!***********************************************************************
      SUBROUTINE ice_bulk_flux_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               liold, linew,                      &
# ifdef MASKING
     &                               rmask, umask, vmask,               &
# endif
# ifdef WET_DRY
     &                               rmask_wet, umask_wet, vmask_wet,   &
# endif
# ifdef ICESHELF
     &                               zice,                              &
# endif
# ifdef WIND_MINUS_CURRENT
     &                               ui, vi,                            &
# endif
     &                               Hair, Pair, Tair, Uwind, Vwind,    &
# ifdef CLOUDS
     &                               cloud,                             &
# endif
     &                               ai, tis,                           &
# ifdef ICE_ALBEDO
     &                               albo, albi, frswp,                 &
# endif
     &                               taux_ai, tauy_ai,                  &
     &                               qswi, qlwi, qswo, qai,             &
# ifdef ICE_DIAGS
#  ifdef ICESHELF
     &                               qswl, qlwl,                        &
#  endif
     &                               qswa, qlwa, qlai, qsei,            &
# endif
     &                               coef_ice_heat, rhs_ice_heat)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
#  endif
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
#  ifdef ICE_ALBEDO
      real(r8), intent(in) :: albo(LBi:,LBj:)
      real(r8), intent(in) :: albi(LBi:,LBj:)
      real(r8), intent(in) :: frswp(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: taux_ai(LBi:,LBj:)
      real(r8), intent(out) :: tauy_ai(LBi:,LBj:)
      real(r8), intent(inout) :: qswi(LBi:,LBj:)
      real(r8), intent(inout) :: qlwi(LBi:,LBj:)
      real(r8), intent(out) :: qswo(LBi:,LBj:)
      real(r8), intent(out) :: qai(LBi:,LBj:)
#  ifdef ICE_DIAGS
#   ifdef ICESHELF
      real(r8), intent(inout) :: qswl(LBi:,LBj:)
      real(r8), intent(inout) :: qlwl(LBi:,LBj:)
#   endif
      real(r8), intent(out) :: qswa(LBi:,LBj:)
      real(r8), intent(out) :: qlwa(LBi:,LBj:)
      real(r8), intent(out) :: qlai(LBi:,LBj:)
      real(r8), intent(out) :: qsei(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: coef_ice_heat(LBi:,LBj:)
      real(r8), intent(out) :: rhs_ice_heat(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vi(LBi:UBi,LBj:UBj,2)
#  endif
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
#  ifdef ICE_ALBEDO
      real(r8), intent(in) :: albo(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: albi(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: frswp(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: taux_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tauy_ai(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: qswi(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: qlwi(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qswo(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qai(LBi:UBi,LBj:UBj)
#  ifdef ICE_DIAGS
#   ifdef ICESHELF
      real(r8), intent(inout) :: qswl(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: qlwl(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(out) :: qswa(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qlwa(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qlai(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qsei(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: coef_ice_heat(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhs_ice_heat(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, listp

!     real(r8), parameter :: io_frac = 0.17_r8   ! from MU
      real(r8) :: PairM, RH, delW, Q, rhoAir, TairC, TairK, Wgus, Wmag

      real(r8) :: qlwi_i, qlai_i, qsei_i
      real(r8) :: dq_i, fqlat1, Qsati, sfc_temp, slp, vap_p_i

      real(r8) :: cff, cff1, cff2

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uair
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vair

# include "set_bounds.h"
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
      IF (PerfectRST(ng) .and. iic(ng).eq.ntstart(ng)) THEN
        listp=liold
      ELSE
        listp=linew
      END IF

# ifdef WIND_MINUS_CURRENT
!
!  Modify near-surface (2m or 10m) effective winds by subtracting the
!  ocean surface current (J Wilkin). See:
!
!  Bye, J.A.T. and J.-O. Wolff, 1999: Atmosphere-ocean momentum exchange
!     in general circulation models. J. Phys. Oceanogr. 29, 671-691.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,MIN(Iend+1,Lm(ng))
          Uair(i,j)=Uwind(i,j)-                                         &
     &              0.5_r8*(ui(i,j,listp)+ui(i+1,j,listp))
        END DO
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          Uair(Iend+1,j)=Uwind(Iend+1,j)-ui(Iend,j,listp)
        END IF
      END DO
      DO i=Istr-1,Iend+1
        DO j=Jstr-1,MIN(Jend+1,Mm(ng))
          Vair(i,j)=Vwind(i,j)-                                         &
     &              0.5_r8*(vi(i,j,listp)+vi(i,j+1,listp))
        END DO
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          Vair(i,Jend+1)=Vwind(i,Jend+1)-vi(i,Jend,listp)
        END IF
      END DO
# else
!
!  Load wind components to local arrays.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Uair(i,j)=Uwind(i,j)
          Vair(i,j)=Vwind(i,j)
        END DO
      END DO
# endif

      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
# ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) THEN
!
!  Zero out fluxes under iceshelf
!
            taux_ai(i,j)       = 0.0_r8
            tauy_ai(i,j)       = 0.0_r8
            qswi(i,j)          = 0.0_r8
            qlwi(i,j)          = 0.0_r8
            qswo(i,j)          = 0.0_r8
            qai(i,j)           = 0.0_r8
#  ifdef ICE_DIAGS
            qswl(i,j)          = 0.0_r8
            qlwl(i,j)          = 0.0_r8
            qswa(i,j)          = 0.0_r8
            qlwa(i,j)          = 0.0_r8
            qlai(i,j)          = 0.0_r8
            qsei(i,j)          = 0.0_r8
#  endif
            coef_ice_heat(i,j) = 0.0_r8
            rhs_ice_heat(i,j)  = 0.0_r8
          ELSE
# endif
!
!  Preparations
!
          Wmag  = SQRT(Uair(i,j)*Uair(i,j)+Vair(i,j)*Vair(i,j))
          Wgus  = 0.5_r8
          delW  = SQRT(Wmag*Wmag+Wgus*Wgus)
          PairM = Pair(i,j)
          TairC = Tair(i,j)
          TairK = TairC+273.16_r8
          RH    = Hair(i,j)
          cff = (1.0007_r8+3.46E-6_r8*PairM)*6.1121_r8*                 &
     &          EXP(17.502_r8*TairC/(240.97_r8+TairC))
          IF (RH.lt.2.0_r8) THEN
!
! RH is relative humidity. First calculate vapor pressure (mb), then
!   calculate specific humidity (kg/kg)
!
            cff = cff*RH
            Q   = 0.62197_r8*(cff/(PairM-0.378_r8*cff))
          ELSE
!
! RH is specific humidity in g/kg
!
            Q=RH*0.001_r8
          END IF
          rhoAir=PairM*100.0_r8/(blk_Rgas*TairK*(1.0_r8+0.61_r8*Q))
!
!  Ice model bulk flux formula.
!
!-----------------------------------------------------------------------
!  Momentum flux over ice
!-----------------------------------------------------------------------
!
          taux_ai(i,j) = rhoAir*cdai(ng)*Uwind(i,j)*Wmag
          tauy_ai(i,j) = rhoAir*cdai(ng)*Vwind(i,j)*Wmag
!
!-----------------------------------------------------------------------
!  Initialize ice heat budgets calculation.
!-----------------------------------------------------------------------
!
!    Note that ocean net short-wave radiation for ice cover are
!    converted to kinematic units.
!
          coef_ice_heat(i,j) = 0._r8
          rhs_ice_heat(i,j)  = 0._r8
!
          IF (ai(i,j,listp).gt.min_a(ng)) THEN
!
!-----------------------------------------------------------------------
!  Sensible heat flux over ice
!-----------------------------------------------------------------------
!
            sfc_temp = tis(i,j) + 273.16_r8
            qsei_i   = rhoAir*blk_Cpa*Ch_i*delW*                        &
     &                 (TairK+0.0098*blk_ZT(ng)-sfc_temp)
            coef_ice_heat(i,j) = coef_ice_heat(i,j) +                   &
     &                           rhoAir*blk_Cpa*Ch_i*delW
            rhs_ice_heat(i,j)  = rhs_ice_heat(i,j)  +                   &
     &                           rhoAir*blk_Cpa*Ch_i*delW*              &
     &                           (TairK+0.0098*blk_ZT(ng))
!
!-----------------------------------------------------------------------
!  Latent heat flux over ice
!-----------------------------------------------------------------------
!
!    Qsati is saturation specific humidity over the ice from
!    Parkinson and Washington (1979)
!
!  Convert back to Pascal from millibars
!
            slp = Pair(i,j)*100._r8
            cff = 611._r8*10._r8**(9.5_r8*(sfc_temp-273.16_r8) /        &
     &                                    (sfc_temp-  7.66_r8))
            Qsati  = 0.62197_r8*cff/(slp-(1._r8-.62197_r8)*cff)
            dq_i   = Q-Qsati
            fqlat1 = rhoAir*Ce_i*delW
            qlai_i = sublim_lat_h*fqlat1*dq_i
!
            rhs_ice_heat(i,j) = rhs_ice_heat(i,j) + qlai_i
!
!-----------------------------------------------------------------------
!  Net long-wave radiation from ice
!-----------------------------------------------------------------------
!
            cff     = Q/(1.0_r8+Q)
            vap_p_i = slp*cff/(0.62197 + cff)
# if defined LONGWAVE
            qlwi_i = emmiss*StefBo*(TairK**4)*                          &
     &               (0.39_r8-0.005_r8*SQRT(vap_p_i))*                  &
     &               (1.0_r8-0.65_r8*cloud(i,j)) -                      &
     &               4._r8*emmiss*StefBo*TairK**3*(sfc_temp-TairK)

            coef_ice_heat(i,j) = coef_ice_heat(i,j) +                   &
     &                           4._r8*emmiss*StefBo*TairK**3
            rhs_ice_heat(i,j)  = rhs_ice_heat(i,j)  -                   &
     &                           emmiss*StefBo*(TairK**4)*              &
     &                           ((0.39_r8-0.005_r8*SQRT(vap_p_i))*     &
     &                            (1.0_r8-0.65_r8*cloud(i,j)) - 4.0_r8)
# elif defined LONGWAVE_OUT
            qlwi_i = qlwi(i,j) - emmiss*StefBo*sfc_temp**4
            coef_ice_heat(i,j) = coef_ice_heat(i,j) +                   &
     &                           4._r8*StefBo*emmiss*sfc_temp**3
            rhs_ice_heat(i,j)  = rhs_ice_heat(i,j) + qlwi(i,j) +        &
     &                           3._r8*emmiss*StefBo*sfc_temp**4
# else
            qlwi_i = 0.0_r8
# endif
# ifdef ICE_DIAGS
            qlwa(i,j) = emmiss*StefBo*sfc_temp**4
            qlwi(i,j) = qlwi_i
# endif
!
!  Correct for Kelvin to Celsius in RHS
!
            rhs_ice_heat(i,j) = rhs_ice_heat(i,j) -                     &
     &                          coef_ice_heat(i,j)*273.16_r8
!
!-----------------------------------------------------------------------
!  Shortwave radiation flux absorbed by ice
!-----------------------------------------------------------------------
!
!  Initialize a few variables
!
# ifdef ICE_ALBEDO
#  ifdef ICE_DIAGS
!
!  Mean reflected shortwave flux into air (W/m2)
!
            qswa(i,j) = qswi(i,j)*(        ai(i,j,listp) *albi(i,j) +   &
     &                             (1.0_r8-ai(i,j,listp))*albo(i,j))
#  endif
!
!  Shortwave flux absorbed by ice and ocean (W/m2)
!
            qswi(i,j) = qswi(i,j)*(1.0_r8-albi(i,j))
!
!  Process shortwave radiation in ice and through ice
!    Penetrating SW = SW*(1-albedo_ice)*frswp (W/m2 ice)
!
#  ifdef ICE_I_O
!
!  Separate penetrating shorwave flux from qswi.
!
!  Shortwave flux at ice base
!
            qswo(i,j) = qswi(i,j) * frswp(i,j)
!
!  Old formula to compute shortwave penetrating ice. Not recommended.
!
!           IF (hsn(i,j,listp) .lt. 1.e-3_r8) THEN
!             qswo(i,j) = qswi(i,j) * io_frac
!           ELSE
!             qswo(i,j) = 0._r8
!           END IF
!
!  Update shortwave radiation absorbed by ice
!
            qswi(i,j) = qswi(i,j) - qswo(i,j)
#  endif
# endif
!
!  Update ice heat budgets.
!
            rhs_ice_heat(i,j) = rhs_ice_heat(i,j) + qswi(i,j)
# ifdef ICE_DIAGS
!
!-----------------------------------------------------------------------
!  Save ice model diagnostics
!-----------------------------------------------------------------------
!
            qlai(i,j) = qlai_i
            qsei(i,j) = qsei_i
# endif
!
!-----------------------------------------------------------------------
!  Net heat flux, positive from ice to atmosphere
!-----------------------------------------------------------------------
!
!  Compute net heat flux from ice to atmosphere - Mellor and Kantha (7)
!    Note that in Mellor and Kantha qai is positive **upwards**,
!    while in ocean model heat flux is positive **downwards**.
!
            qai(i,j)  = -(qsei_i + qlai_i + qswi(i,j) + qlwi_i)
          ELSE
            qai(i,j)  = 0._r8
            qswi(i,j) = 0._r8
            qswo(i,j) = 0._r8
# ifdef ICE_DIAGS
            qswa(i,j) = 0._r8
            qlwa(i,j) = 0._r8
            qlwi(i,j) = 0._r8
            qlai(i,j) = 0._r8
            qsei(i,j) = 0._r8
# endif
          END IF
# ifdef ICESHELF
          END IF
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          taux_ai)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tauy_ai)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qai)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qswo)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          coef_ice_heat)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhs_ice_heat)
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    taux_ai, tauy_ai)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    qai, qswo)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    coef_ice_heat, rhs_ice_heat)
# endif
!
      RETURN
      END SUBROUTINE ice_bulk_flux_tile
#endif
      END MODULE ice_bulk_flux_mod
