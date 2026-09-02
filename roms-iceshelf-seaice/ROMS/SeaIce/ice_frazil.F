#include "cppdefs.h"
      MODULE ice_frazil_mod
#if defined ICE_THERMO && !defined TS_FIXED
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the frazil ice growth in the water when the   !
!  water temperature gets below freezing. It adjusts the water         !
!  temperature and salinity accordingly.                               !
!                                                                      !
!  Reference: Steele et al. (1989). JPO, 19, 139-147.                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ice_frazil

      CONTAINS

      SUBROUTINE ice_frazil (ng, tile)

      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_ice
      USE mod_stepping

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 96, __LINE__, __FILE__)
# endif
!
      CALL ice_frazil_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nnew(ng),                                   &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
# endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % z_r,                             &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % t,                              &
     &                      ICE(ng) % wfr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 96, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE ice_frazil
!
!***********************************************************************
      subroutine ice_frazil_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nnew,                                 &
# ifdef MASKING
     &                            rmask,                                &
# endif
# ifdef WET_DRY
     &                            rmask_wet,                            &
# endif
     &                            Hz, z_r, rho,                         &
     &                            t, wfr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_iounits
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nnew

# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: wfr(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(out) :: wfr(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable definitions
!
      integer :: i, j, k, itrc

      real(r8), parameter :: Lhat = 79.2_r8
      real(r8), parameter :: r = 0.5_r8

      real(r8) :: t_freeze, s1, z1, sice, t_fr, gamma_k, delta_wfr
      real(r8) :: cff
!
!  A few choices for freezing temperature formula
!
!  Freezing temperature (Gill, 1982)
!     t_freeze(s,z) = -0.0575*s + 1.710523d-3*sqrt(s)**3 -
!                     2.154996d-4*s**2 + 0.000753*z
!
!  Freezing temperature (Steele et al. 1989) (default)
!     t_freeze(s,z) = frln*s + 0.000759*z
!
!  Another simple formula
!
!     t_freeze(s,z) = frln*s
!
# include "set_bounds.h"
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          wfr(i,j) = 0.0_r8
        END DO
      END DO
!
!  Old ice model
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          DO k=1,N(ng)
# ifdef WET_DRY
            IF (rmask_wet(i,j) .ne. 0.0_r8) THEN
# elif defined MASKING
            IF (rmask(i,j) .ne. 0.0_r8) THEN
# endif
              t_fr = frln*t(i,j,k,nnew,isalt) + 0.000759_r8*z_r(i,j,k)
              delta_wfr = 0._r8
              IF (t(i,j,k,nnew,itemp) .lt. t_fr) THEN
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                    (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r) -    &
     &                     t_fr)
                delta_wfr = gamma_k * Hz(i,j,k) *                       &
     &                      (rho0 + rho(i,j,k)) / rhoice(ng)
                wfr(i,j) = wfr(i,j) + delta_wfr
                t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) +             &
     &            gamma_k * (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r))
                sice = MIN(sice_ref,t(i,j,k,nnew,isalt))
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) *             &
     &                                (1.0_r8 + gamma_k) - sice*gamma_k
              ELSE IF (wfr(i,j)            .gt. 0._r8 .and.             &
     &                 t(i,j,k,nnew,itemp) .gt. t_fr) THEN
!
!  Use heat at this level to melt some ice from below.
!    gamma_k becomes negative here.
!
                gamma_k = (t_fr - t(i,j,k,nnew,itemp)) /                &
     &                    (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r) -    &
     &                     t_fr)
                delta_wfr = gamma_k * Hz(i,j,k) *                       &
     &                      (rho0 + rho(i,j,k)) / rhoice(ng)
                IF ((wfr(i,j) + delta_wfr) .gt. 0._r8) THEN
                  wfr(i,j) = wfr(i,j) + delta_wfr
                ELSE
                  gamma_k = -wfr(i,j) * rhoice(ng) /                    &
     &                      (Hz(i,j,k)*(rho0+rho(i,j,k)))
                  wfr(i,j) = 0.0_r8
                  delta_wfr = -wfr(i,j)
                END IF
                t(i,j,k,nnew,itemp) = t(i,j,k,nnew,itemp) +             &
     &            gamma_k * (Lhat + t(i,j,k,nnew,itemp)*(1.0_r8 - r))
                t(i,j,k,nnew,isalt) = t(i,j,k,nnew,isalt) *             &
     &            (1.0_r8 + gamma_k)
              END IF
!
! Sanity check.
!
              IF (wfr(i,j) .lt. 0.0_r8) THEN
                WRITE(stdout,*) 'Trouble in ice_frazil', i, j, k,       &
     &            t(i,j,k,nnew,itemp), t(i,j,k,nnew,isalt),             &
     &            t_fr, wfr(i,j), gamma_k, Hz(i,j,k)
                exit_flag = 10
              END IF
!
# ifdef T_PASSIVE
              cff = delta_wfr * rhoice(ng) / (rho0+rho(i,j,k))
              DO itrc=NAT+1,NT(ng)
                t(i,j,k,nnew,itrc) = t(i,j,k,nnew,itrc)*                &
     &                               (Hz(i,j,k)-cff)/Hz(i,j,k)
              END DO
# endif
# if defined WET_DRY || defined MASKING
            END IF
# endif
          END DO
!
          wfr(i,j) = wfr(i,j) / dt(ng)
        END DO
      END DO
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,itemp))
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          t(:,:,:,nnew,isalt))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wfr)
      END IF
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    wfr)
# endif
!
      RETURN
      END SUBROUTINE ice_frazil_tile
#endif
      END MODULE ice_frazil_mod
