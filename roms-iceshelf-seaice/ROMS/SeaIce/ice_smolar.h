!
!***********************************************************************
      SUBROUTINE ice_advect (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 95, __LINE__, __FILE__)
#endif
!
! ---------------------------------------------------------------------
!  Advect the ice concentration.
! ---------------------------------------------------------------------
!
      CALL ice_advect_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      linew(ng), liold(ng), liunw(ng),            &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#endif
#ifdef WET_DRY
     &                      GRID(ng) % rmask_wet,                       &
     &                      GRID(ng) % umask_wet,                       &
     &                      GRID(ng) % vmask_wet,                       &
#endif
#ifdef ICESHELF
     &                      GRID(ng) % zice,                            &
#endif
#ifndef ICE_UPWIND
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
#endif
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % omn,                             &
     &                      ICE(ng) % ui,                               &
     &                      ICE(ng) % vi,                               &
     &                      ICE(ng) % ai,                               &
     &                      ICE(ng) % hi,                               &
     &                      ICE(ng) % hsn,                              &
     &                      ICE(ng) % ti,                               &
     &                      ICE(ng) % ageice,                           &
     &                      ICE(ng) % enthalpi,                         &
     &                      ICE(ng) % hage,                             &
     &                      ICE(ng) % wdiv)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 95, __LINE__, __FILE__)
#endif
!
      RETURN
      END SUBROUTINE ice_advect
!
!***********************************************************************
      SUBROUTINE ice_advect_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            linew, liold, liunw,                  &
#ifdef MASKING
     &                            rmask, umask, vmask,                  &
#endif
#ifdef WET_DRY
     &                            rmask_wet, umask_wet, vmask_wet,      &
#endif
#ifdef ICESHELF
     &                            zice,                                 &
#endif
#ifndef ICE_UPWIND
     &                            pm, pn,                               &
#endif
     &                            on_u, om_v, omn, ui, vi,              &
     &                            ai, hi, hsn, ti, ageice,              &
     &                            enthalpi, hage, wdiv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
      USE i2dbc_mod, ONLY : i2dbc_ice, i2dbc_tice
      USE bc_2d_mod, ONLY : bc_r2d_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: linew, liold, liunw

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
# endif
# ifndef ICE_UPWIND
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
# endif
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
      real(r8), intent(inout) :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
      real(r8), intent(inout) :: ti(LBi:,LBj:,:)
      real(r8), intent(inout) :: ageice(LBi:,LBj:,:)
      real(r8), intent(inout) :: enthalpi(LBi:,LBj:,:)
      real(r8), intent(inout) :: hage(LBi:,LBj:,:)
      real(r8), intent(inout) :: wdiv(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifndef ICE_UPWIND
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ti(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ageice(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: enthalpi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hage(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: wdiv(LBi:UBi,LBj:UBj)
#endif
!
!  Local variables
!
      integer :: i, j
!
#include "set_bounds.h"
!
! ---------------------------------------------------------------------
!  Advect the ice coverage.
! ---------------------------------------------------------------------
!
      CALL ice_adv_var_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       linew, liold, liunw,                       &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
#endif
#ifdef ICESHELF
     &                       zice,                                      &
#endif
#ifndef ICE_UPWIND
     &                       pm, pn,                                    &
#endif
     &                       on_u, om_v, omn, ui, vi,                   &
     &                       ai)
!
! ---------------------------------------------------------------------
!  Advect the ice thickness.
! ---------------------------------------------------------------------
!
      CALL ice_adv_var_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       linew, liold, liunw,                       &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
#endif
#ifdef ICESHELF
     &                       zice,                                      &
#endif
#ifndef ICE_UPWIND
     &                       pm, pn,                                    &
#endif
     &                       on_u, om_v, omn, ui, vi,                   &
     &                       hi)
!
! ---------------------------------------------------------------------
!  Advect the snow thickness.
! ---------------------------------------------------------------------
!
#ifdef ICE_THERMO
      CALL ice_adv_var_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       linew, liold, liunw,                       &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
#endif
#ifdef ICESHELF
     &                       zice,                                      &
#endif
#ifndef ICE_UPWIND
     &                       pm, pn,                                    &
#endif
     &                       on_u, om_v, omn, ui, vi,                   &
     &                       hsn)
!
! ---------------------------------------------------------------------
!  Advect the interior ice temperature.
! ---------------------------------------------------------------------
!
      CALL ice_adv_var_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       linew, liold, liunw,                       &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
#endif
#ifdef ICESHELF
     &                       zice,                                      &
#endif
#ifndef ICE_UPWIND
     &                       pm, pn,                                    &
#endif
     &                       on_u, om_v, omn, ui, vi,                   &
     &                       enthalpi)
!
! ---------------------------------------------------------------------
!  Advect the ice age.
! ---------------------------------------------------------------------
!
      CALL ice_adv_var_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       linew, liold, liunw,                       &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef WET_DRY
     &                       rmask_wet, umask_wet, vmask_wet,           &
#endif
#ifdef ICESHELF
     &                       zice,                                      &
#endif
#ifndef ICE_UPWIND
     &                       pm, pn,                                    &
#endif
     &                       on_u, om_v, omn, ui, vi,                   &
     &                       hage)
!
! ---------------------------------------------------------------------
!  Compute other variables.
! ---------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
!
!  Ice divergence.
!
          wdiv(i,j) = (hi(i,j,linew) - hi(i,j,liold)) / dt(ng)
!
!  Ice interior temperature.
!
          ti(i,j,linew) = enthalpi(i,j,linew) /                         &
     &                    MAX(hi(i,j,linew), 1.0E-6_r8)
          IF (hi(i,j,linew).le.min_h(ng)) THEN
            enthalpi(i,j,linew) = 0.0_r8
            ti(i,j,linew) = 0.0_r8
          END IF
!
!  Ice age.
!
          ageice(i,j,linew) = hage(i,j,linew) /                         &
     &                        MAX(hi(i,j,linew),1.0E-6_r8)
          IF (hi(i,j,linew).le.min_h(ng)) THEN
            hage(i,j,linew) = 0.0_r8
            ageice(i,j,linew) = 0.0_r8
          END IF
        END DO
      END DO
!
! ---------------------------------------------------------------------
!  Clear ice velocities where ice is not presented.
! ---------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=IstrP,Iend
          IF ((ai(i,j,linew)+ai(i-1,j,linew)).le.min_a(ng)) THEN
            ui(i,j,linew) = 0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrP,Jend
        DO i=Istr,Iend
          IF ((ai(i,j,linew)+ai(i,j-1,linew)).le.min_a(ng)) THEN
            vi(i,j,linew) = 0.0_r8
          END IF
        END DO
      END DO
!
! ---------------------------------------------------------------------
!  Boundary data exchange.
! ---------------------------------------------------------------------
!
!  Set lateral boundary conditions.
!
      CALL i2dbc_ice  (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 liold, linew,                                    &
     &                 ui, vi, ai, hi, hsn)
!
      CALL i2dbc_tice (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 liold, linew,                                    &
     &                 ui, vi, hi, ti, enthalpi)
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  ageice(:,:,linew))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  hage(:,:,linew))
!
!  Apply periodic boundary conditions.
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
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          enthalpi(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ageice(:,:,linew))
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          hage(:,:,linew))
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ui(:,:,linew))
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vi(:,:,linew))
      END IF
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ai(:,:,linew), hi(:,:,linew), hsn(:,:,linew))
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ti(:,:,linew), enthalpi(:,:,linew),           &
     &                    ageice(:,:,linew), hage(:,:,linew))
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ui(:,:,liunw), vi(:,:,liunw))
# endif
#endif
      RETURN
      END SUBROUTINE ice_advect_tile
!
!=======================================================================!
      SUBROUTINE ice_adv_var_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             linew, liold, liunw,                 &
#ifdef MASKING
     &                             rmask, umask, vmask,                 &
#endif
#ifdef WET_DRY
     &                             rmask_wet, umask_wet, vmask_wet,     &
#endif
#ifdef ICESHELF
     &                             zice,                                &
#endif
#ifndef ICE_UPWIND
     &                             pm, pn,                              &
#endif
     &                             on_u, om_v, omn, ui, vi,             &
     &                             scr)
!=======================================================================!

      USE mod_param
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: linew, liold, liunw
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
# endif
# ifndef ICE_UPWIND
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
# endif
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: ui(LBi:,LBj:,:)
      real(r8), intent(in) :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: scr(LBi:,LBj:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
# endif
# ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
# endif
# ifndef ICE_UPWIND
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: scr(LBi:UBi,LBj:UBj,2)
#endif
!
! Local variable definitions
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j

      real(r8) :: Cu_crss, Cu
      real(r8) :: rateu, ratev, rateyiu, ratexiv, uspeed, vspeed
      real(r8) :: cff

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ar
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aflxv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
!
#include "set_bounds.h"
!
#ifndef ICE_UPWIND
      IF (EWperiodic(ng)) THEN
        Imin=Istr-1
        Imax=Iend+1
      ELSE
        Imin=MAX(Istr-1,1     )
        Imax=MIN(Iend+1,Lm(ng))
      END IF
      IF (NSperiodic(ng)) THEN
        Jmin=Jstr-1
        Jmax=Jend+1
      ELSE
        Jmin=MAX(Jstr-1,1     )
        Jmax=MIN(Jend+1,Mm(ng))
      END IF
#else
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
#endif
!
! Upstream advection scheme
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax+1
          aflxu(i,j)=on_u(i,j)*                                         &
     &               (MAX(0.0_r8,ui(i,j,liunw))*scr(i-1,j,liold) +      &
     &                MIN(0.0_r8,ui(i,j,liunw))*scr(i  ,j,liold))
        END DO
      END DO
      DO j=Jmin,Jmax+1
        DO i=Imin,Imax
          aflxv(i,j)=om_v(i,j)*                                         &
     &               (MAX(0.0_r8,vi(i,j,liunw))*scr(i,j-1,liold) +      &
     &                MIN(0.0_r8,vi(i,j,liunw))*scr(i,j  ,liold))
        END DO
      END DO
!
! Compute new tracer concentration
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
!
          ar(i,j)=1.0_r8/omn(i,j)
          aif(i,j)=scr(i,j,liold) -                                     &
     &      dtice(ng)*(aflxu(i+1,j  ) - aflxu(i,j) +                    &
     &                 aflxv(i  ,j+1) - aflxv(i,j))*ar(i,j)
#ifdef MASKING
          aif(i,j) = aif(i,j)*rmask(i,j)
#endif
#ifdef WET_DRY
          aif(i,j) = aif(i,j)*rmask_wet(i,j)
#endif
#ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) aif(i,j) = 0.0_r8
#endif
        END DO
      END DO
!
! Set values at the open boundaries
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Istr-1,j)=aif(Istr,j)
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jmin,Jmax
            aif(Iend+1,j)=aif(Iend,j)
          END DO
        END IF
      END IF
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Imin,Imax
            aif(i,Jstr-1)=aif(i,Jstr)
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Imin,Imax
            aif(i,Jend+1)=aif(i,Jend)
          END DO
        END IF
      END IF
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          aif(Istr-1,Jstr-1)=0.5_r8*(aif(Istr-1,Jstr)+aif(Istr,Jstr-1))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          aif(Istr-1,Jend+1)=0.5_r8*(aif(Istr-1,Jend)+aif(Istr,Jend+1))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          aif(Iend+1,Jstr-1)=0.5_r8*(aif(Iend+1,Jstr)+aif(Iend,Jstr-1))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          aif(Iend+1,Jend+1)=0.5_r8*(aif(Iend+1,Jend)+aif(Iend,Jend+1))
        END IF
      END IF
#ifndef ICE_UPWIND
!
! Antidiffusive corrector step:
! ------------- --------- -----
!
! This is needed to avoid touching "aif" under land mask. Note that only
! aif(i,j) and aif(i-1,j) are allowed to appear explicitly in the code
! segment below. This is OK because if either of them masked, then "ui"
! is 0 at that point, and therefore no additional masking is required. 
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend+1
          FE(i,j)=0.5*                                                  &
# ifdef MASKING
     &            vmask(i,j)*                                           &
# endif
# ifdef WET_DRY
     &            vmask_wet(i,j)*                                       &
# endif
     &            (aif(i,j)-aif(i,j-1))
          FX(i,j)=0.5*                                                  &
# ifdef MASKING
     &            umask(i,j)*                                           &
# endif
# ifdef WET_DRY
     &            umask_wet(i,j)*                                       &
# endif
     &            (aif(i,j)-aif(i-1,j))
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend+1
          cff = aif(i,j)+aif(i-1,j)
          rateu=(aif(i,j)-aif(i-1,j)) /                                 &
     &          (SIGN(1._r8, cff)*MAX(ABS(cff), tol))
!
          cff = aif(i  ,j) + FE(i  ,j+1) - FE(i  ,j) +                  &
     &          aif(i-1,j) + FE(i-1,j+1) - FE(i-1,j)
          rateyiu=(FE(i,j+1) + FE(i,j) + FE(i-1,j+1) + FE(i-1,j)) /     &
     &            (SIGN(1._r8, cff)*MAX(ABS(cff), tol))
!
          Cu=0.5*dtice(ng)*(pm(i,j)+pm(i-1,j))*ui(i,j,liunw)
!
          Cu_crss=0.5*dtice(ng)*0.0625 *                                &
     &      (pn(i-1,j+1)+pn(i,j+1)+pn(i-1,j-1)+pn(i,j-1)) *             &
     &      (vi(i-1,j+1,liunw)+vi(i,j+1,liunw)+                         &
     &       vi(i-1,j  ,liunw)+vi(i,j  ,liunw))
!
          uspeed=rateu*(ABS(ui(i,j,liunw))-Cu*ui(i,j,liunw)) -          &
     &           rateyiu*Cu_crss*ui(i,j,liunw)
!
          aflxu(i,j)=on_u(i,j)*(MAX(0.,uspeed)*aif(i-1,j) +             &
     &                          MIN(0.,uspeed)*aif(i  ,j))
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
          cff = aif(i,j)+aif(i,j-1)
          ratev=(aif(i,j)-aif(i,j-1)) /                                 &
     &          (SIGN(1._r8, cff)*MAX(ABS(cff), tol))
!
          cff = aif(i,j  ) + FX(i+1,j  ) - FX(i,j  ) +                  &
     &          aif(i,j-1) + FX(i+1,j-1) - FX(i,j-1)
          ratexiv=(FX(i+1,j) + FX(i,j) + FX(i+1,j-1) + FX(i,j-1)) /     &
     &            (SIGN(1._r8, cff)*MAX(ABS(cff), tol))
!
          Cu=0.5*dtice(ng)*(pn(i,j)+pn(i,j-1))*vi(i,j,liunw)
!
          Cu_crss=0.5*dtice(ng)*0.0625 *                                &
     &      (pm(i+1,j)+pm(i+1,j-1)+pm(i-1,j)+pm(i-1,j-1)) *             &
     &      (ui(i,j  ,liunw)+ui(i+1,j  ,liunw)+                         &
     &       ui(i,j-1,liunw)+ui(i+1,j-1,liunw))
!
          vspeed=ratev*(ABS(vi(i,j,liunw))-Cu*vi(i,j,liunw)) -          &
     &           ratexiv*Cu_crss*vi(i,j,liunw)
!
          aflxv(i,j)=om_v(i,j)*(MAX(0.,vspeed)*aif(i,j-1) +             &
     &                          MIN(0.,vspeed)*aif(i,j  ))
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          aif(i,j)=aif(i,j) -                                           &
     &             dtice(ng)*pm(i,j)*pn(i,j)*                           &
     &             (aflxu(i+1,j)-aflxu(i,j)+aflxv(i,j+1)-aflxv(i,j))
# ifdef MASKING
          aif(i,j)=aif(i,j)*rmask(i,j)
# endif
# ifdef WET_DRY
          aif(i,j)=aif(i,j)*rmask_wet(i,j)
# endif
# ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) aif(i,j)=0.
# endif
        END DO
      END DO
#endif /* !defined ICE_UPWIND */

      DO j=Jstr,Jend
        DO i=Istr,Iend
            scr(i,j,linew)=aif(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE ice_adv_var_tile
