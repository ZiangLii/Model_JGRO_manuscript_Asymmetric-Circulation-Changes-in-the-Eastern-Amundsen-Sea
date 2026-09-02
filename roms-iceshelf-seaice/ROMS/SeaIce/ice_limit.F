#include "cppdefs.h"
      MODULE ice_limit_mod
#ifdef ICE_MODEL
!
!=================================================== W. Paul Budgell ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Limit ice tracer variables, for instance to not allow negative ice  !
!  concentration or thickness.                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC ice_limit
!
      CONTAINS
!
      SUBROUTINE ice_limit (ng, tile)
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
!
! ---------------------------------------------------------------------
!  Limit ice model variables.
! ---------------------------------------------------------------------
!
      CALL ice_limit_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     liold(ng), linew(ng),                        &
     &                     ICE(ng) % ui,                                &
     &                     ICE(ng) % vi,                                &
     &                     ICE(ng) % ai,                                &
     &                     ICE(ng) % hi,                                &
     &                     ICE(ng) % hsn,                               &
     &                     ICE(ng) % ageice,                            &
     &                     ICE(ng) % ti,                                &
     &                     ICE(ng) % enthalpi)
!
      RETURN
      END SUBROUTINE ice_limit
!
!***********************************************************************
      SUBROUTINE ice_limit_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           liold, linew,                          &
     &                           ui, vi, ai, hi, hsn, ageice,           &
     &                           ti, enthalpi)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE i2dbc_mod, ONLY : i2dbc_ice, i2dbc_tice
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liold, linew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ageice(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpi(LBi:UBi,LBj:UBj,2)
# endif
!
! Local variable definitions
!
      integer :: i, j
!
# include "set_bounds.h"
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          ai(i,j,linew) = MAX(MIN(ai(i,j,linew),max_a(ng)),0.0_r8)
          hi(i,j,linew) = MAX(hi(i,j,linew),0.0_r8)
          hsn(i,j,linew) = MAX(hsn(i,j,linew),0.0_r8)
          ageice(i,j,linew) =                                           &
     &      MAX(0.0_r8,MIN(ageice(i,j,linew),                           &
     &                     ageice(i,j,liold)+dtice(ng)/86400._r8))
          ti(i,j,linew) = MAX(ti(i,j,linew),-50.0_r8)
          IF (ai(i,j,linew) .le. min_a(ng) .or.                         &
     &        hi(i,j,linew) .le. 0.0_r8) THEN
            ai(i,j,linew) = 0.0_r8
            hi(i,j,linew) = 0.0_r8
            hsn(i,j,linew) = 0.0_r8
            ageice(i,j,linew) = 0.0_r8
          END IF
        END DO
      END DO
!
!     CALL i2dbc_ice  (ng, tile,                                        &
!    &                 LBi, UBi, LBj, UBj,                              &
!    &                 IminS, ImaxS, JminS, JmaxS,                      &
!    &                 liold, linew,                                    &
!    &                 ui, vi, ai, hi, hsn)
!
!     CALL i2dbc_tice (ng, tile,                                        &
!    &                 LBi, UBi, LBj, UBj,                              &
!    &                 IminS, ImaxS, JminS, JmaxS,                      &
!    &                 liold, linew,                                    &
!    &                 ui, vi, hi, ti, enthalpi)
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ai(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hi(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hsn(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ti(:,:,linew))
!
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ai(:,:,linew), hi(:,:,linew),                 &
     &                    hsn(:,:,linew),ti(:,:,linew))
# endif
      RETURN
      END SUBROUTINE ice_limit_tile
#endif
      END MODULE ice_limit_mod
