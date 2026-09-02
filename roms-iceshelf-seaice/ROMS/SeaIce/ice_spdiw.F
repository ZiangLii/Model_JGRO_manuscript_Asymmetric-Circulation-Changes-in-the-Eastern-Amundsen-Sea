#include "cppdefs.h"
      MODULE ice_spdiw_mod
#ifdef ICE_MODEL
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module computes the magnitude of the shear between the ice     !
!  and the surface water. In this case, the surface water is defined   !
!  as the water in a surface mixed layer, so that velocity must be     !
!  computed first.                                                     !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ice_spdiw

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_spdiw (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
#  ifdef LMD_SKPP
      USE mod_mixing
#  endif
      USE mod_ice
      USE mod_stepping

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6, __LINE__, __FILE__)
# endif
      CALL ice_spdiw_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng), liuol(ng),                         &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
# ifdef LMD_SKPP
     &                     MIXING(ng) % hsbl,                           &
# endif
     &                     ICE(ng) % ui,                                &
     &                     ICE(ng) % vi,                                &
     &                     ICE(ng) % uwater,                            &
     &                     ICE(ng) % vwater,                            &
     &                     ICE(ng) % spd_iw)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE ice_spdiw
!
!***********************************************************************
      SUBROUTINE ice_spdiw_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs, liuol,                           &
     &                           z_r, z_w,                              &
     &                           u, v,                                  &
# ifdef LMD_SKPP
     &                           hsbl,                                  &
# endif
     &                           ui, vi,                                &
     &                           uwater, vwater, spd_iw)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, liuol

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(out) :: uwater(LBi:,LBj:)
      real(r8), intent(out) :: vwater(LBi:,LBj:)
      real(r8), intent(out) :: spd_iw(LBi:,LBj:)
# else
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  ifdef LMD_SKPP
      real(r8), intent(in) :: hsbl(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: uwater(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: vwater(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: spd_iw(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
      integer :: nlio, nbotu, nbotv, k

      real(r8) :: mlio, dml, totml

      integer, dimension(IminS:ImaxS,JminS:JmaxS) :: nbot
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uw
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vw
!
# include "set_bounds.h"
!
      DO j=MAX(Jstr-2,0),MIN(Jend+2,Mm(ng)+1)
        DO i=MAX(Istr-2,0),MIN(Iend+2,Lm(ng)+1)
!         sl_dpth = lmd_epsilon*(z_w(i,j,N(ng)) - hsbl(i,j))
# ifdef LMD_SKPP
!
! hsbl is now a positive quantity
!
          mlio = MIN(-hsbl(i,j),-10._r8)
# else
          mlio = -10._r8
# endif
          nbot(i,j) = 1
!
          DO k=N(ng),1,-1
            IF (z_r(i,j,k).lt.mlio) THEN
              nbot(i,j) = MIN(k,N(ng))
              nbot(i,j) = MAX(nbot(i,j),1)
              EXIT
            END IF
          END DO
!
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=MAX(Istr-1,1),Iend+1
          nlio = 0
          nbotu = NINT(0.5_r8*(nbot(i-1,j)+nbot(i,j)))
          nbotu = MAX(MIN(nbotu,N(ng)),1)
          uw(i,j) = 0._r8
          totml = 0._r8
          DO k=N(ng),nbotu,-1
            nlio = nlio + 1
            dml = 0.5_r8*(z_w(i-1,j,k)-z_w(i-1,j,k-1) +                 &
     &                    z_w(i  ,j,k)-z_w(i  ,j,k-1))
            uw(i,j) = uw(i,j) + u(i,j,k,nrhs)*dml
            totml = totml + dml
          END DO
          uw(i,j) = uw(i,j) / totml
!         uw(i,j) = u(i,j,N,nrhs)
        END DO
      END DO

      DO j=MAX(Jstr-1,1),Jend+1
        DO i=Istr,Iend
          nlio = 0
          nbotv = NINT(0.5_r8*(nbot(i,j-1)+nbot(i,j)))
          nbotv = MAX(MIN(nbotv,N(ng)),1)
          vw(i,j) = 0._r8
          totml = 0._r8
          DO k=N(ng),nbotv,-1
            nlio = nlio + 1
            dml = 0.5_r8*(z_w(i,j-1,k)-z_w(i,j-1,k-1) +                 &
     &                    z_w(i,j  ,k)-z_w(i,j  ,k-1))
            vw(i,j) = vw(i,j) + v(i,j,k,nrhs)*dml
            totml = totml + dml
          END DO
          vw(i,j) = vw(i,j) / totml
!         vw(i,j) = v(i,j,N,nrhs)
        END DO
      END DO

      DO j=Jstr,Jend
        DO i=Istr,Iend
          spd_iw(i,j) = 0.5*SQRT((uw(i  ,j)-ui(i  ,j,liuol) +           &
     &                            uw(i+1,j)-ui(i+1,j,liuol))**2 +       &
     &                           (vw(i,j  )-vi(i,j  ,liuol) +           &
     &                            vw(i,j+1)-vi(i,j+1,liuol))**2)
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrP,Iend
          uwater(i,j) = uw(i,j)
        END DO
      END DO
      DO j=JstrP,Jend
        DO i=Istr,Iend
          vwater(i,j) = vw(i,j)
        END DO
      END DO
!
!  Apply boundary conditions.
!
        CALL bc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    spd_iw)
        CALL bc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    uwater)
        CALL bc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    vwater)
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    spd_iw, uwater, vwater)
# endif
!
      RETURN
      END SUBROUTINE ice_spdiw_tile
#endif
      END MODULE ice_spdiw_mod
