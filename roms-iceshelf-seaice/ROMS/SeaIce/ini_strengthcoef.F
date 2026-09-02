#include "cppdefs.h"
      MODULE ini_strengthcoef_mod
#ifdef ICE_STRENGTH_QUAD
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes pstar_grid based on the grid spacing.      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ini_strengthcoef

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_strengthcoef (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ice
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL ini_strengthcoef_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % grdscl,                    &
     &                            ICE(ng) % pstar_grid)

      RETURN
      END SUBROUTINE ini_strengthcoef
!
!***********************************************************************
      SUBROUTINE ini_strengthcoef_tile (ng, tile, model,                &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  grdscl,                         &
     &                                  pstar_grid)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: grdscl(LBi:,LBj:)
      real(r8), intent(inout) :: pstar_grid(LBi:,LBj:)
# else
      real(r8), intent(in) :: grdscl(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: pstar_grid(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax

      integer :: i, j
      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Scale ice strength according to the grid size. Scale by 1 km.
!-----------------------------------------------------------------------
!
      cff=pstar_e(ng)/1000.0_r8
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          pstar_grid(i,j) = cff*grdscl(i,j)
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
     &                          pstar_grid)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pstar_grid)
# endif
!
      RETURN
      END SUBROUTINE ini_strengthcoef_tile
#endif
      END MODULE ini_strengthcoef_mod
