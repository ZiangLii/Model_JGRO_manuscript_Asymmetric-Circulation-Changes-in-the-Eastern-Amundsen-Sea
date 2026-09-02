#include "cppdefs.h"
      MODULE albedo_mod
#if defined ICE_ALBEDO && !defined ICE_ALB_FILE
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes albedo of ocean and ice.                      !
!                                                                      !
!  Additionally, it computes the fraction of shortwave radtion         !
!    penetrating ice base, so that shortwave reaching ice base is      !
!                                                                      !
!      SW_ice_base = SW*(1-albedo_ice)*frswp                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: albedo

      CONTAINS
!
!***********************************************************************
      SUBROUTINE albedo (ng, tile)
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17, __LINE__, __FILE__)
# endif
      CALL albedo_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  liold(ng), linew(ng),                           &
# ifdef ICESHELF
     &                  GRID(ng) % zice,                                &
# endif
     &                  ICE(ng) % ai,                                   &
     &                  ICE(ng) % hi,                                   &
     &                  ICE(ng) % hsn,                                  &
     &                  ICE(ng) % tis,                                  &
# ifdef CLOUDS
     &                  GRID(ng) % latr,                                &
     &                  GRID(ng) % lonr,                                &
     &                  FORCES(ng) % cloud,                             &
# endif
# ifdef ICE_LMD_SHORTWAVE
     &                  ICE(ng) % rfaci,                                &
# endif
     &                  FORCES(ng) % albi,                              &
     &                  FORCES(ng) % albo,                              &
     &                  FORCES(ng) % frswp)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE albedo
!
!***********************************************************************
      SUBROUTINE albedo_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        liold, linew,                             &
# ifdef ICESHELF
     &                        zice,                                     &
# endif
     &                        ai, hi, hsn, tis,                         &
# ifdef CLOUDS
     &                        latr, lonr, cloud,                        &
# endif
# ifdef ICE_LMD_SHORTWAVE
     &                        rfaci,                                    &
# endif
     &                        albi, albo, frswp)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE dateclock_mod,   ONLY : caldate
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
      integer, intent(in) :: liold
      integer, intent(in) :: linew
!
# ifdef ASSUMED_SHAPE
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
#  ifdef CLOUDS
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#  endif
#  ifdef ICE_LMD_SHORTWAVE
      real(r8), intent(out) :: rfaci(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: albi(LBi:,LBj:)
      real(r8), intent(out) :: albo(LBi:,LBj:)
      real(r8), intent(out) :: frswp(LBi:,LBj:)

# else
#  ifdef ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
#  ifdef CLOUDS
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: lonr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICE_LMD_SHORTWAVE
      real(r8), intent(out) :: rfaci(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: albi(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: albo(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: frswp(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, listp
!
# ifdef ICE_ALB_EC92
      real(r8) :: alb_s = 0.83_r8
# elif defined ICE_ALB_CSIM
!
! DuVivier et al., 2018
!
!   Visible ice albedo for h > ahmax; CESM2/CICE5.1 default
!     (cf. 0.73 in Briegleb et al., 2004)
!
      real(r8), parameter :: albicev=0.78_r8
!
!   Near-ir ice albedo for h > ahmax; CESM2/CICE5.1 default
!     (cf. 0.33 in Briegleb et al., 2004)
!
      real(r8), parameter :: albicei=0.36_r8
!
!   Cold snow albedo, visible; CESM2/CICE5.1 default
!     (cf. 0.96 in Briegleb et al., 2004)
!
      real(r8), parameter :: albsnowv=0.98_r8
!
!   Cold snow albedo, near IR; CESM2/CICE5.1 default
!     (cf. 0.68 in Briegleb et al., 2004)
!
      real(r8), parameter :: albsnowi=0.70_r8
!
!   Thickness above which ice albedo is const (m); CESM2/CICE5.1 default
!     (cf. 0.5 m in Briegleb et al., 2004)
!
      real(r8), parameter :: ahmax=0.3_r8
!
!   Parameter for fractional snow area (m); CESM2/CICE5.1 default
!     (cf. 0.02 m in Briegleb et al., 2004)
!
      real(r8), parameter :: snowpatch=0.02_r8
!
!   Change in temp to give dalb_mlt albedo (C); CESM2/CICE5.1 default
!     (cf. 1.0 C in Briegleb et al., 2004)
!
      real(r8), parameter :: dT_melt=1.0_r8
!
!   Albedo change per dT_melt in temp for ice; CESM2/CICE5.1 default
!     (cf. -0.075 in Briegleb et al., 2004)
!
      real(r8), parameter :: dalb_mlt=-0.075_r8
!
!   Albedo vis change per dT_melt in temp for snow; CESM2/CICE5.1 default
!     (cf. -0.1 in Briegleb et al., 2004)
!
      real(r8), parameter :: dalb_mltv=-0.1_r8
!
!   Albedo ir change per dT_melt in temp for snow; CESM2/CICE5.1 default
!     (cf. -0.15 in Briegleb et al., 2004)
!
      real(r8), parameter :: dalb_mlti=-0.15_r8
!
!   Melting temperature, ice top surface (C)
!
      real(r8), parameter :: Timelt=0.0_r8
!
!   Melting temperature, snow top surface (C)
!
      real(r8), parameter :: Tsmelt=0.0_r8
      real(r8), parameter :: puny=1.0e-11_r8
!
!   Visual fraction of total shortwave (Aumont et al., 2015)
!
      real(r8), parameter :: fv=0.43_r8
!
!   Fraction of penetrating solar rad (visible); CESM2/CICE5.1 default
!     (cf. 0.7 in Briegleb et al., 2004)
!
      real(r8), parameter :: i0vis=0.70_r8

      real(r8) :: fh, alb, alvni, alini, dTs, fT
      real(r8) :: alvns, alins, asnow, alvn, alin
# else
      real(r8), parameter :: alb_i_dry=0.65_r8
      real(r8), parameter :: alb_i_wet=0.60_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.72_r8
!
!  Was these values
!
!     real(r8), parameter :: alb_i_thick=0.54_r8
!     real(r8), parameter :: alb_s_dry=0.83_r8
!     real(r8), parameter :: alb_s_wet=0.70_r8
!
!  Values from Ungermann et al, 2017
!
!     real(r8), parameter :: alb_i_thick=0.71_r8
!     real(r8), parameter :: alb_s_dry=0.86_r8
!     real(r8), parameter :: alb_s_wet=0.79_r8
!
!  Try something intermediate
!
!     real(r8), parameter :: alb_i_thick=0.65_r8
!     real(r8), parameter :: alb_s_dry=0.86_r8
!     real(r8), parameter :: alb_s_wet=0.72_r8
# endif
!
!  Default water albedo
!
      real(r8), parameter :: alb_o=0.06_r8
      real(r8), parameter :: alb_oz=0.065_r8

      real(r8) :: i_thick, s_thick
      real(r8) :: hours, days, dangle, hangle, langle
      real(r8) :: sin2, cos2, cosz, albdir
# ifdef ICE_ALB_DMEAN
      integer, parameter :: ifrac = 24
      real(r8), parameter :: absh2o = 0.09_r8
      real(r8), parameter :: fraci = 1._r8/REAL(ifrac,r8)
      integer :: npar
      real(r8) :: stot, sads, ssurf, sdir, sdif, altdeg, cfac
# endif
      real(r8) :: cff1, cff2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
! PURPOSE:
!   computes albedo over snow/ice/water
!-----------------------------------------------------------------------
!

      IF (PerfectRST(ng) .and. iic(ng).eq.ntstart(ng)) THEN
        listp = liold
      ELSE
        listp = linew
      END IF
!
# if defined ICE_ALB_DIR || defined ICE_ALB_DMEAN
!
! add a correction term to correct the time difference between Gregorian
!   time and local time.
!
!   CW 2022-01-07
!
      CALL caldate (tdays(ng), yd_dp=days, h_dp=hours)
      days = MOD(days-1._r8,365._r8)
      dangle = 2.0_r8*pi*days/365._r8
      dangle = .006918_r8 +                                             &
     &  .070257_r8*SIN(dangle)       - .399912_r8*COS(dangle)       +   &
     &  .000907_r8*SIN(2._r8*dangle) - .006758_r8*COS(2._r8*dangle) +   &
     &  .001480_r8*SIN(3._r8*dangle) - .002697_r8*COS(3._r8*dangle)

      hangle = pi*(hours-12.0_r8)/12.0_r8
# endif

      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
# ifdef ICESHELF
          IF (zice(i,j).ne.0.0_r8) THEN
!
!  Zero out fluxes under iceshelf
!
            albi(i,j) = 0.0_r8
            albo(i,j) = 0.0_r8
            frswp(i,j) = 0.0_r8
#  ifdef ICE_LMD_SHORTWAVE
            rfaci(i,j) = 0.0_r8
#  endif
          ELSE
# endif
!
! Calculate the ice/snow albedo
!
          i_thick = hi(i,j,listp)/(ai(i,j,listp)+0.001)
          s_thick = hsn(i,j,listp)/(ai(i,j,listp)+0.001)
# if defined ICE_ALB_CSIM
!
! We follow ice_shortwave.F90/compute_albedos subroutine (CESM2/CICE5.1)
!
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
!
! bare ice, thickness dependence
!
            fh = MIN(ATAN(i_thick*4.0_r8)/ATAN(ahmax*4.0_r8), 1.0_r8)
            alb = alb_o*(1.0_r8-fh)
            alvni = albicev*fh + alb
            alini = albicei*fh + alb
!
! bare ice, temperature dependence
!
            dTs = Timelt - tis(i,j)
!
! fT=0 for dTs>=dT_melt, =-1 for dTs=0
!
            fT = MIN(dTs/dT_melt-1.0_r8, 0.0_r8)
            alvni = MAX(alb_o, alvni - dalb_mlt*fT)
!
! Ice albedo decreases with warming, capped by ocean albedo
!
            alini = MAX(alb_o, alini - dalb_mlt*fT)

            IF (s_thick > puny) THEN
!
! Snow on ice, temperature dependence
!
              alvns = MAX(alb_o, albsnowv - dalb_mltv*fT)
              alins = MAX(alb_o, albsnowi - dalb_mlti*fT)
!
! Note: the CESM2/CICE5.1 code (ice_shortwave.F90/compute_albedos) does
!       not impose lower cap here, but it seems safer to do so for the
!       case of tis>>Timelt (=> fT<<0)
!
              asnow = s_thick / (s_thick + snowpatch)
!
! Combined ice/snow visual albedo
!
              alvn = alvni*(1.0_r8-asnow) + alvns*asnow
!
! Combined ice/snow near IR albedo
!
              alin = alini*(1.0_r8-asnow) + alins*asnow
            ELSE
              asnow = 0.0_r8
!
! Combined ice/snow visual albedo
!
              alvn = alvni
!
! Combined ice/snow near IR albedo
!
              alin = alini
            END IF
!
! Combined ice/snow albedo wrt total incident SW
! (absorbed SW = (1-albedo_i)*SW)
!
            albi(i,j) = fv*alvn + (1.0_r8-fv)*alin
          ELSE
            albi(i,j) = alb_o
          END IF
# elif defined ICE_ALB_EC92
!
! Ice and snow albedo is calculated from Ebert and Curry,1992b
!
          IF (s_thick .gt. 0._r8) THEN
            albi(i,j) = alb_s
          ELSE
!
! Approximated values for ice alb (depends on COSZ, but small variation)
!
            IF (ai(i,j,listp) .gt. min_a(ng)) THEN
              i_thick = MAX(i_thick, 1.0E-6_r8)
              IF (i_thick.ge.2._r8) THEN
                albi(i,j) = 0.561632_r8
              ELSE IF (i_thick.ge.1._r8) THEN
                albi(i,j) = 0.07616_r8*i_thick+0.414492_r8
              ELSE
                albi(i,j) = 0.082409_r8*LOG(i_thick) + 0.485472_r8
              END IF
            ELSE
              albi(i,j) = alb_o
            END IF
          END IF
# else
          cff1 = alb_s_wet - alb_s_dry
          cff2 = alb_i_wet - alb_i_dry
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            IF (hsn(i,j,listp).gt.0._r8) THEN
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albi(i,j) = cff1*(tis(i,j)+1.0_r8) + alb_s_dry
              ELSE
                albi(i,j) = alb_s_dry
              END IF
            ELSE
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albi(i,j) = cff2*(tis(i,j)+1.0_r8) + alb_i_dry
              ELSE
                albi(i,j) = alb_i_dry
              END IF
            END IF
          ELSE
            albi(i,j) = alb_o
          END IF
# endif
!
! Compute ocean albedo
!
# if defined ICE_ALB_DIR
          sin2 = SIN(latr(i,j)*deg2rad)*SIN(dangle)
          cos2 = COS(latr(i,j)*deg2rad)*COS(dangle)
#  ifdef ICE_ALB_GREGORIAN
!
!  Correction for local time if Gregorian reference time is used
!
          langle = MOD(lonr(i,j), 360._r8)*rad2deg
#  else
          langle = 0.0_r8
#  endif
!
!  Cosine of the zenith angle
!
          cosz = MIN(1._r8,MAX(0._r8,sin2+cos2*COS(hangle + langle)))
!
!  Albedo for direct light
!
          albdir = 0.03_r8*EXP(0.742_r8*ACOS(cosz)**2.866_r8)
!
!  Albedo weighted by cloud fraction
!
          albo(i,j) = 1._r8 - ((1.0_r8-cloud(i,j))*(1.0_r8-albdir) +    &
     &                                 cloud(i,j) *(1.0_r8-alb_oz))
# elif defined ICE_ALB_DMEAN
!
! ----------------------------------------------------------------------
!  compute 24 hrs mean albedo at the marine surface layer
! ----------------------------------------------------------------------
!
          sin2 = SIN(latr(i,j)*deg2rad)*SIN(dangle)
          cos2 = COS(latr(i,j)*deg2rad)*COS(dangle)
          stot = 0._r8
          sads = 0._r8
          DO npart=1,ifrac
!
!  Hour angle, in radians
!
            hangle = pi*(REAL(npart,r8)/12._r8-1._r8)
!
!  Cosine of the zenith angle
!
            cosz = MIN(1._r8,MAX(0._r8,sin2+cos2*COS(hangle)))
!
!  Direct radiation component
!
            sdir = 0.7_r8**MIN(100._r8,1._r8/(cosz+tol))
!
!  Diffusive radiation component
!
            sdif = ((1._r8-absh2o)-sdir)*.5_r8
!
!  Solar noon altitude in degrees
!
            altdeg = MAX(0._r8,ASIN(sin2+cos2))*rad2deg
!
!  Cloudiness correction
!
            cfac = (1._r8-0.62_r8*cloud(i,j)+0.0019_r8*altdeg)
!
            ssurf = (sdir+sdif)*cfac+tol
!
            stot = stot+ssurf
!
!  Albedo for direct light
!
            albdir = 0.03_r8*EXP(0.742_r8*ACOS(cosz)**2.866_r8)
!
!  Radiation weighted sum of direct albedo
!
            sads = sads+albdir*ssurf
          END DO
!
!  Daily mean co-albedo over water
!
          albdir = sads/stot
!
!  Daily mean albedo weighted by cloud fraction
!
          albo(i,j) = 1._r8 - ((1.0_r8-cloud(i,j))*(1.0_r8-albdir) +    &
     &                                 cloud(i,j) *(1.0_r8-alb_oz))
# else
          albo(i,j) = alb_o
# endif
!
!  Calculate shortwave penetration factor
!
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
# ifdef ICE_ALB_CSIM
!
! Note: We assume that only the visual part incident on snow-uncovered
!       ice can penetrate through the surface layer (if not reflected)
!
            frswp(i,j) = (1.0_r8-alvni)*(1.0_r8-asnow)*i0vis*fv /       &
     &                   (1.0_r8-albi(i,j))
# else
            IF (i_thick .le. 0.1_r8) THEN
              cff1 = i_thick*5._r8
            ELSE
              cff1 = 0.4_r8 + i_thick
            END IF
            cff1 = cff1 + s_thick*20._r8
            frswp(i,j) = fv*EXP(-cff1)
# endif
          ELSE
            frswp(i,j) = 1.0_r8
          END IF
# ifdef ICE_LMD_SHORTWAVE
!
!  Calculate correction factor for lmd_r1
!
          cff1 =  1.0_r8 - albo(i,j)
          cff2 = (1.0_r8 - albi(i,j))*frswp(i,j)
          rfaci(i,j) = (1.0_r8-ai(i,j,listp))*cff1 /                    &
     &      MAX(tol,((1.0_r8-ai(i,j,listp))*cff1 + ai(i,j,listp)*cff2))
# endif
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
# ifdef ICE_LMD_SHORTWAVE
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rfaci)
# endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albo)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albi)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          frswp)
      END IF
# ifdef DISTRIBUTE
#  ifdef ICE_LMD_SHORTWAVE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rfaci)
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albo, albi, frswp)
# endif
      RETURN
!
      END SUBROUTINE albedo_tile
#endif
      END module albedo_mod
