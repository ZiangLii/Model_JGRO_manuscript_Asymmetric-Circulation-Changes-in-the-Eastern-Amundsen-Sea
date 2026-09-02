#include "cppdefs.h"
      MODULE ice_enthalpi_mod
#ifdef ICE_MODEL
!
!=================================================== W. Paul Budgell ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine computes ice enthalpi (ti*hi) prior to advection.   !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ice_enthalpi

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ice_enthalpi (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 95, __LINE__, __FILE__)
#endif
!
! ---------------------------------------------------------------------
!  Compute ice enthalpi.
! ---------------------------------------------------------------------
!
      CALL ice_enthalpi_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        ICE(ng)%hi,                               &
     &                        ICE(ng)%ti,                               &
     &                        ICE(ng)%enthalpi)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 95, __LINE__, __FILE__)
#endif
!
      RETURN
      END SUBROUTINE ice_enthalpi
!
!***********************************************************************
      SUBROUTINE ice_enthalpi_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              hi, ti,                             &
     &                              enthalpi)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpi(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpi(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
!
# include "set_bounds.h"
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          enthalpi(i,j,liold(ng)) = hi(i,j,liold(ng))*ti(i,j,liold(ng))
          enthalpi(i,j,linew(ng)) = hi(i,j,linew(ng))*ti(i,j,linew(ng))
          IF(hi(i,j,liold(ng)).le.min_h(ng)) THEN
            enthalpi(i,j,liold(ng)) = 0.0_r8
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          enthalpi(:,:,liold(ng)))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          enthalpi(:,:,linew(ng)))
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      enthalpi(:,:,liold(ng)),                    &
     &                      enthalpi(:,:,linew(ng)))
# endif
!
      RETURN
      END SUBROUTINE ice_enthalpi_tile
#endif
      END MODULE ice_enthalpi_mod
