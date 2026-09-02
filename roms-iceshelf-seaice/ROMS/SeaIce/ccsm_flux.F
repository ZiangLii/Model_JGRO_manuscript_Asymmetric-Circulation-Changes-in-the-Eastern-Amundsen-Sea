#include "cppdefs.h"
      MODULE ccsm_flux_mod
#if defined BULK_FLUXES && defined CCSM_FLUXES
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!       CCSM's flux_mod.F90 by B. Kauffman                             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ccsm_flux

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ccsm_flux (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
# ifdef EVAP_SHALLOW
      USE mod_coupling
# endif
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
      CALL ccsm_flux_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nrhs(ng),                                    &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
# endif
# ifdef WET_DRY
     &                     GRID(ng) % rmask_wet,                        &
     &                     GRID(ng) % umask_wet,                        &
     &                     GRID(ng) % vmask_wet,                        &
# endif
     &                     OCEAN(ng) % t,                               &
# ifdef WIND_MINUS_CURRENT
     &                     OCEAN(ng) % u,                               &
     &                     OCEAN(ng) % v,                               &
# endif
     &                     FORCES(ng) % Hair,                           &
     &                     FORCES(ng) % Pair,                           &
     &                     FORCES(ng) % Tair,                           &
     &                     FORCES(ng) % Uwind,                          &
     &                     FORCES(ng) % Vwind,                          &
# ifdef CLOUDS
     &                     FORCES(ng) % cloud,                          &
# endif
# ifdef COARE_TAYLOR_YELLAND
     &                     FORCES(ng) % Hwave,                          &
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
     &                     FORCES(ng) % Pwave_top,                      &
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
     &                     FORCES(ng) % Lwave,                          &
# endif
# ifdef EVAP_SHALLOW
     &                     GRID(ng) % h,                                &
     &                     COUPLING(ng) % Zt_avg1,                      &
# endif
     &                     FORCES(ng) % rain,                           &
     &                     FORCES(ng) % lhflx,                          &
     &                     FORCES(ng) % lrflx,                          &
     &                     FORCES(ng) % shflx,                          &
     &                     FORCES(ng) % srflx,                          &
# ifdef ICE_MODEL
#  ifdef ICE_SNOWFALL
     &                     FORCES(ng) % snow,                           &
#  endif
#  ifdef ICE_ALBEDO
     &                     FORCES(ng) % albo,                           &
#  endif
     &                     FORCES(ng) % wsf,                            &
     &                     FORCES(ng) % sustr_aw,                       &
     &                     FORCES(ng) % svstr_aw,                       &
     &                     FORCES(ng) % qswi,                           &
#  ifdef LONGWAVE_OUT
     &                     FORCES(ng) % qlwi,                           &
#  endif
#  ifdef ICE_DIAGS
     &                     FORCES(ng) % qswl,                           &
     &                     FORCES(ng) % qlwl,                           &
#  endif
# endif
# ifdef EMINUSP
     &                     FORCES(ng) % evap,                           &
# endif
     &                     FORCES(ng) % stflux,                         &
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17, __LINE__, __FILE__)
# endif
!
      RETURN
      END SUBROUTINE ccsm_flux
!
!***********************************************************************
      SUBROUTINE ccsm_flux_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nrhs,                                  &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef WET_DRY
     &                           rmask_wet, umask_wet, vmask_wet,       &
# endif
     &                           t,                                     &
# ifdef WIND_MINUS_CURRENT
     &                           u, v,                                  &
# endif
     &                           Hair, Pair, Tair, Uwind, Vwind,        &
# ifdef CLOUDS
     &                           cloud,                                 &
# endif
# ifdef COARE_TAYLOR_YELLAND
     &                           Hwave,                                 &
# endif
# if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
     &                           Pwave_top,                             &
# endif
# if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
     &                           Lwave,                                 &
# endif
# ifdef EVAP_SHALLOW
     &                           h, Zt_avg1,                            &
# endif
     &                           rain, lhflx, lrflx, shflx, srflx,      &
# ifdef ICE_MODEL
#  ifdef ICE_SNOWFALL
     &                           snow,                                  &
#  endif
#  ifdef ICE_ALBEDO
     &                           albo,                                  &
#  endif
     &                           wsf,                                   &
     &                           sustr_aw, svstr_aw,                    &
     &                           qswi,                                  &
#  ifdef LONGWAVE_OUT
     &                           qlwi,                                  &
#  endif
#  ifdef ICE_DIAGS
     &                           qswl, qlwl,                            &
#  endif
# endif
# ifdef EMINUSP
     &                           evap,                                  &
# endif
     &                           stflux, sustr, svstr)
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
      integer, intent(in) :: nrhs
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
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#  endif
      real(r8), intent(in) :: Hair(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
      real(r8), intent(in) :: Tair(LBi:,LBj:)
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(inout) :: Uwind(LBi:,LBj:)
      real(r8), intent(inout) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  endif
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:,LBj:)
#  endif
#  ifdef COARE_TAYLOR_YELLAND
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
#  endif
#  if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
#  endif
#  if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: rain(LBi:,LBj:)
#  ifdef EVAP_SHALLOW
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: Zt_avg1(LBi:,LBj:)
#  endif

      real(r8), intent(inout) :: lhflx(LBi:,LBj:)
      real(r8), intent(inout) :: lrflx(LBi:,LBj:)
      real(r8), intent(inout) :: shflx(LBi:,LBj:)
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflux(LBi:,LBj:,:)

#  ifdef EMINUSP
      real(r8), intent(out) :: evap(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)

#  ifdef ICE_MODEL
#   ifdef ICE_SNOWFALL
      real(r8), intent(in) :: snow(LBi:,LBj:)
#   endif
#   ifdef ICE_ALBEDO
      real(r8), intent(in) :: albo(LBi:,LBj:)
#   endif
      real(r8), intent(out) :: wsf(LBi:,LBj:)
      real(r8), intent(out) :: sustr_aw(LBi:,LBj:)
      real(r8), intent(out) :: svstr_aw(LBi:,LBj:)
      real(r8), intent(out) :: qswi(LBi:,LBj:)
#   ifdef LONGWAVE_OUT
      real(r8), intent(out) :: qlwi(LBi:,LBj:)
#   endif
#   ifdef ICE_DIAGS
      real(r8), intent(out) :: qswl(LBi:,LBj:)
      real(r8), intent(out) :: qlwl(LBi:,LBj:)
#   endif
#  endif

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
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),3)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),3)
#  endif
      real(r8), intent(in) :: Hair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Tair(LBi:UBi,LBj:UBj)
#  ifdef WIND_MINUS_CURRENT
      real(r8), intent(inout) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  endif
#  ifdef CLOUDS
      real(r8), intent(in) :: cloud(LBi:UBi,LBj:UBj)
#  endif
#  ifdef COARE_TAYLOR_YELLAND
      real(r8), intent(in) :: Hwave(LBi:UBi,LBj:UBj)
#  endif
#  if defined COARE_TAYLOR_YELLAND || defined COARE_OOST
      real(r8), intent(in) :: Pwave_top(LBi:UBi,LBj:UBj)
#  endif
#  if !defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)
      real(r8), intent(in) :: Lwave(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: rain(LBi:UBi,LBj:UBj)
#  ifdef EVAP_SHALLOW
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Zt_avg1(LBi:UBi,LBj:UBj)
#  endif

      real(r8), intent(inout) :: lhflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: lrflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: shflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: stflux(LBi:UBi,LBj:UBj,NT(ng))

#  ifdef EMINUSP
      real(r8), intent(out) :: evap(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)

#  ifdef ICE_MODEL
#   ifdef ICE_SNOWFALL
      real(r8), intent(in) :: snow(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ICE_ALBEDO
      real(r8), intent(in) :: albo(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(out) :: wsf(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: sustr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr_aw(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qswi(LBi:UBi,LBj:UBj)
#   ifdef LONGWAVE_OUT
      real(r8), intent(out) :: qlwi(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ICE_DIAGS
      real(r8), intent(out) :: qswl(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: qlwl(LBi:UBi,LBj:UBj)
#   endif
#  endif

# endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8), parameter :: r3 = 1.0_r8/3.0_r8
!
!  Ratio of molecular weight of water to dry air
!
      real(r8), parameter :: epsilon = 0.622_r8

      real(r8) :: Hscale, Hscale_inv
      real(r8) :: RH, cff, cff1, cff2, cff3
# ifdef EVAP_SHALLOW
      real(r8), parameter :: coef_shallow = 5.28_r8
!     real(r8), parameter :: coef_shallow = 4.69_r8
      real(r8) :: cffe
# endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uair
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vair
      real(r8) :: fsnow

! ------------------------------------------------------------------
! !DESCRIPTION:
!     wrapper to atm/ocn flux calculation
!
! !REMARKS:
!     All data must be on the ocean domain (note: a domain includes a
!     particular decomposition).
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - first version
!
! ------------------------------------------------------------------
!
! !IROUTINE: srfflx_ao -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - brought in from cpl5.
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as in Bill Large
!
!INTERFACE: ------------------------------------------------------------
!
!SUBROUTINE srfflx_ao_tile(blk_ZX, Uwind, Vwind, Tair, Hair, rhoAir,    &
!    &                     u, v, t, rmask, shflx, lhflx, lwup,          &
!    &                     evap, taux, tauy, tref, qref, duu10n)
!
! INPUT/OUTPUT PARAMETERS:
!
!--- input arguments --------------------------------
! integer(IN) ,intent(in) :: rmask(imax) ! ocn domain mask 0
! real(R8), intent(in) :: blk_ZX (imax)  ! atm level height  (m)
! real(R8), intent(in) :: Uwind (imax)   ! atm u wind        (m/s)
! real(R8), intent(in) :: Vwind (imax)   ! atm v wind        (m/s)
! real(R8), intent(in) :: Tair (imax)    ! atm potential T   (K)
! real(R8), intent(in) :: Hair (imax)    ! atm rel humidity  (frac)
! real(R8), intent(in) :: rhoAir (imax)  ! atm air density   (kg/m^3)
! real(R8), intent(in) :: u (imax)       ! ocn u-velocity    (m/s)
! real(R8), intent(in) :: v (imax)       ! ocn v-velocity    (m/s)
! real(R8), intent(in) :: t (imax)       ! ocn temperature   (K)
!
!--- output arguments -------------------------------
! real(R8), intent(out) :: shflx (imax)  ! heat flux: sensible    (W/m^2)
! real(R8), intent(out) :: lhflx (imax)  ! heat flux: latent      (W/m^2)
! real(R8), intent(out) :: lwup (imax)   ! heat flux: lw upward   (W/m^2)
! real(R8), intent(out) :: evap (imax)   ! water flux: evap  ((kg/s)/m^2)
! real(R8), intent(out) :: taux (imax)   ! surface stress, zonal      (N)
! real(R8), intent(out) :: tauy (imax)   ! surface stress, maridional (N)
! real(R8), intent(out) :: tref (imax)   ! diag:  2m ref height T     (K)
! real(R8), intent(out) :: qref (imax)   ! diag:  2m ref humidity (kg/kg)
! real(R8), intent(out) :: duu10n (imax) ! diag: 10m wspeed sqrd (m/s)^2
!
!EOP
!
!--- local constants --------------------------------
!
      real(R8),parameter :: umin  =  0.5  ! minimum wind speed   (m/s)
      real(R8),parameter :: zref  = 10.0  ! ref height           (m)
      real(R8),parameter :: ztref =  2.0  ! ref height for air T (m)
!
!--- local variables --------------------------------
!
      real(r8) :: lwup    ! upward longwave radiation
      real(r8) :: vmag    ! surface wind magnitude   (m/s)
      real(r8) :: thvbot  ! virtual temperature      (K)
      real(r8) :: ssq     ! sea surface humidity     (kg/kg)
      real(r8) :: Tseak   ! SST in Kelvins           (K)
      real(r8) :: delt    ! potential T difference   (K)
      real(r8) :: delq    ! humidity difference      (kg/kg)
      real(r8) :: stable  ! stability factor
      real(r8) :: rdn     ! sqrt of neutral exchange coeff (momentum)
      real(r8) :: rhn     ! sqrt of neutral exchange coeff (heat)
      real(r8) :: ren     ! sqrt of neutral exchange coeff (water)
      real(r8) :: rd      ! sqrt of exchange coefficient (momentum)
      real(r8) :: re      ! sqrt of exchange coefficient (water)
      real(r8) :: ustar   ! ustar
      real(r8) :: qstar   ! qstar
      real(r8) :: tstar   ! tstar
      real(r8) :: hol     ! H (at blk_ZX) over L
      real(r8) :: xsq     ! ?
      real(r8) :: xqq     ! ?
      real(r8) :: psimh   ! stability function at blk_ZX (momentum)
      real(r8) :: psixh   ! stability function at blk_ZX (heat/water)
      real(r8) :: psix2   ! stability function at ztref ref height
      real(r8) :: alz     ! ln(blk_ZX/zref)
      real(r8) :: u10n    ! 10m neutral wind
      real(r8) :: tau     ! stress at blk_ZX
      real(r8) :: cpair   ! specific heat of moist air
      real(r8) :: fac     ! vertical interpolation factor
      real(r8) :: Hlv     ! Latent heat of vaporization at surface
      real(r8) :: Qair    ! Specific humidity [kg/kg]
!
      real(r8), dimension(IminS:ImaxS) :: rhoAir
      real(r8), dimension(IminS:ImaxS) :: TairK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: LHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: SHeat
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Taux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Tauy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u10n_tmp
!
!--- local functions --------------------------------
!
      real(r8)    :: qsat    ! function: satur humididty of air (kg/m^3)
      real(r8)    :: cdn     ! function: neutral drag coeff at 10m
      real(r8)    :: psimhu  ! function: unstable part of psimh
      real(r8)    :: psixhu  ! function: unstable part of psimx
      real(r8)    :: tsk     ! water temperature (K)
      real(r8)    :: tak     ! air temperature (K)
      real(r8)    :: Tk      ! dummy arg ~ air temperature (K)
      real(r8)    :: Umps    ! dummy arg ~ wind velocity (m/s)
      real(r8)    :: xd      ! dummy arg ~ ?

!
! Stefan-Boltzmann constant ~ W/m^2/K^4
!
      real(r8), parameter :: shr_const_stebol  = 5.67e-8_r8
!
! Boltzmann's constant ~ J/K/molecule
!
      real(r8), parameter :: shr_const_boltz   = 1.38065e-23_r8
!
! Avogadro's number ~ molecules/kmole
!
      real(r8), parameter :: shr_const_avogad  = 6.02214e26_r8
!
! Universal gas constant ~ J/K/kmole
!
      real(r8), parameter :: shr_const_rgas =                           &
     &         shr_const_avogad*shr_const_boltz
!
! molecular weight dry air ~ kg/kmole
!
      real(r8), parameter :: shr_const_mwdair  = 28.966_r8
!
! molecular weight water vapor
!
      real(r8), parameter :: shr_const_mwwv = 18.016_r8
!
! Water vapor gas constant ~ J/K/kg
!
      real(r8), parameter :: shr_const_rwv =                            &
     &         shr_const_rgas/shr_const_mwwv
!
! Dry air gas constant     ~ J/K/kg
!
      real(r8), parameter :: shr_const_rdair   =                        &
     &         shr_const_rgas/shr_const_mwdair
      real(r8), parameter :: shr_const_zvir    =                        &
     &         (shr_const_rwv/shr_const_rdair)-1.0_r8
!
! Von Karman constant
!
      real(r8), parameter :: shr_const_karman  = 0.4_r8
!
! specific heat of dry air   ~ J/kg/K
!
      real(r8), parameter :: shr_const_cpdair  = 1.00464e3_r8
!
! specific heat of water vap ~ J/kg/K
!
      real(r8), parameter :: shr_const_cpwv    = 1.810e3_R8
!
! CPWV/CPDAIR - 1.0
!
      real(r8), parameter :: shr_const_cpvir   =                        &
     &         (shr_const_cpwv/shr_const_cpdair)-1.0_r8
!
! latent heat of evaporation ~ J/kg
!
      real(r8), parameter :: shr_const_latvap  = 2.501e6_r8


      qsat(Tk)   = 640380.0 / EXP(5107.4/Tk)
      cdn(Umps)  = 0.0027 / Umps + 0.000142 + 0.0000764 * Umps
      psimhu(xd) = LOG((1.0+xd*(2.0+xd))*(1.0+xd*xd)/8.0) -             &
     &             2.0*ATAN(xd) + 1.571
      psixhu(xd) = 2.0 * LOG((1.0 + xd*xd)/2.0)

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
! PURPOSE:
!   computes atm/ocn surface fluxes
!
! NOTES:
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + shflx + lhflx
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
!
! ASSUMPTIONS:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!-----------------------------------------------------------------------
!
# ifdef WIND_MINUS_CURRENT
!
!  Modify near-surface (2m or 10m) effective winds by subtracting the
!  ocean surface current (J Wilkin). See:
!
!  Bye, J.A.T. and J.-O. Wolff, 1999: Atmosphere-ocean momentum exchange
!  in general circulation models. J. Phys. Oceanogr. 29, 671-691.
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          Uair(i,j)=Uwind(i,j)-                                         &
     &              0.5_r8*(u(i,j,N(ng),nrhs)+u(i+1,j,N(ng),nrhs))
          Vair(i,j)=Vwind(i,j)-                                         &
     &              0.5_r8*(v(i,j,N(ng),nrhs)+v(i,j+1,N(ng),nrhs))
        END DO
      END DO
# else
!
!  Load wind components to local arrays.
!
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
          Uair(i,j)=Uwind(i,j)
          Vair(i,j)=Vwind(i,j)
        END DO
      END DO
# endif

      Hscale=rho0*Cp
      Hscale_inv=1.0_r8/Hscale
!
!  Compute Atmosphere-ocean fluxes using a bulk flux parameterization.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
!
!  Compute specific humidities (kg/kg).
!
          IF (Hair(i,j).lt.2.0_r8) THEN  !RH fraction
            cff  = (1.0007_r8+3.46E-6_r8*Pair(i,j))*6.1121_r8*          &
     &             EXP(17.502_r8*Tair(i,j)/(240.97_r8+Tair(i,j)))
            cff  = cff*Hair(i,j)
            Qair = 0.62197_r8*(cff/(Pair(i,j)-0.378_r8*cff))
          ELSE  !RH input was actually specific humidity in g/kg
            Qair=Hair(i,j)/1000.0_r8  !Spec Hum (kg/kg)
          END IF
!
!--- compute some needed quantities ---
!
          vmag     = MAX(umin, SQRT( (Uair(i,j)**2 + Vair(i,j)**2)) )
          TairK(i) = Tair(i,j) + 273.16_r8
          TseaK    = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          thvbot   = TairK(i) * (1.0 + shr_const_zvir * Qair)
!
!  Virtual temp (K)
!
          rhoAir(i) = Pair(i,j)*100.0_r8 /                              &
     &                (blk_Rgas*TairK(i)*(1.0_r8+0.61_r8*Qair))
!
!  Moist air density (kg/m3).
!
          ssq    = 0.98_r8 * qsat(TseaK) / rhoAir(i)
!
!  Sea surf hum (kg/kg)
!
          delt   = Tair(i,j) - t(i,j,N(ng),nrhs,itemp)
!
!  Pot temp diff (K)
!
          delq   = Qair - ssq
!
!  Spec hum dif (kg/kg)
!
          alz    = LOG(blk_ZW(ng)/zref)
          cpair  = shr_const_cpdair*(1.0 + shr_const_cpvir*ssq)
!
!  first estimate of Z/L and ustar, tstar and qstar
!
!--- neutral coefficients, z/L = 0.0 ---
!
          stable = 0.5_r8 + SIGN(0.5_r8 , delt)
          rdn    = SQRT(cdn(vmag))
          rhn    = (1.0_r8-stable) * 0.0327_r8 + stable * 0.018_r8
          ren    = 0.0346_r8
!
!--- ustar, tstar, qstar ---
!
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq
!
!--- compute stability & evaluate all stability functions ---
!
          hol    = shr_const_karman*g*blk_ZT(ng)*                       &
     &             (tstar/thvbot+qstar/(1.0_r8/shr_const_zvir+Qair)) /  &
     &             ustar**2
          hol    = SIGN( MIN(ABS(hol),10.0_r8), hol )
          stable = 0.5_r8 + SIGN(0.5_r8 , hol)
          xsq    = MAX(SQRT(ABS(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
          xqq    = SQRT(xsq)
          psimh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psimhu(xqq)
          psixh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psixhu(xqq)
!
!--- shift wind speed using old coefficient ---
!
          rd   = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          u10n = vmag * rd / rdn
!
!--- update transfer coeffs at 10m and neutral stability ---
!
          rdn = SQRT(cdn(u10n))
          ren = 0.0346_r8
          rhn = (1.0_r8-stable)*0.0327_r8 + stable * 0.018_r8
!
!--- shift all coeffs to measurement height and stability ---
!
          rd = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          rh = rhn / (1.0_r8 + rhn/shr_const_karman*(alz-psixh))
          re = ren / (1.0_r8 + ren/shr_const_karman*(alz-psixh))
!
!--- update ustar, tstar, qstar using updated, shifted coeffs --
!
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq
!
!------------------------------------------------------------
! iterate to converge on Z/L, ustar, tstar and qstar
!------------------------------------------------------------
!
!--- compute stability & evaluate all stability functions ---
!
          hol    = shr_const_karman*g*blk_ZQ(ng)*                       &
     &             (tstar/thvbot+qstar/(1.0/shr_const_zvir+Qair)) /     &
     &             ustar**2
          hol    = SIGN( MIN(ABS(hol),10.0_r8), hol )
          stable = 0.5_r8 + SIGN(0.5_r8 , hol)
          xsq    = MAX(SQRT(ABS(1.0_r8 - 16.0_r8*hol)) , 1.0_r8)
          xqq    = SQRT(xsq)
          psimh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psimhu(xqq)
          psixh  = -5.0_r8*hol*stable + (1.0_r8-stable)*psixhu(xqq)
!
!--- shift wind speed using old coeffs ---
!
          rd   = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          u10n = vmag * rd/rdn
!
!--- update transfer coeffs at 10m and neutral stability ---
!
          rdn = SQRT(cdn(u10n))
          ren = 0.0346_r8
          rhn = (1.0_r8 - stable)*0.0327_r8 + stable * 0.018_r8
!
!--- shift all coeffs to measurement height and stability ---
!
          rd = rdn / (1.0_r8 + rdn/shr_const_karman*(alz-psimh))
          rh = rhn / (1.0_r8 + rhn/shr_const_karman*(alz-psixh))
          re = ren / (1.0_r8 + ren/shr_const_karman*(alz-psixh))
!
!--- update ustar, tstar, qstar using updated, shifted coeffs ---
!
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq
!
!--- RD: save u10n to u10n_tmp of dim IminS:ImaxS,JminS:JmaxS
! otherwise there is an array overflow !!!
!
          u10n_tmp(i,j) = u10n
!
!-----------------------------------------------------------------------
!  Compute the fluxes
!-----------------------------------------------------------------------
!
          tau = rhoAir(i) * ustar * ustar
!
!-----------------------------------------------------------------------
!  Momentum flux
!-----------------------------------------------------------------------
!
          Taux(i,j) = tau * Uair(i,j) / vmag
          Tauy(i,j) = tau * Vair(i,j) / vmag
!
!  Sensible and latent heat flux
!
          SHeat(i,j) = cpair * tau * tstar / ustar
          LHeat(i,j) = shr_const_latvap * tau * qstar / ustar
!
# ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Save original SRad
!-----------------------------------------------------------------------
!
          qswi(i,j) = srflx(i,j)*Hscale
#  ifdef LONGWAVE_OUT
          qlwi(i,j) = lrflx(i,j)*Hscale
#  endif
#  ifdef ICE_ALBEDO
!
!  Correct shortwave radiation for ocean albedo
!
          srflx(i,j) = srflx(i,j)*(1.0_r8-albo(i,j))
#  endif
# endif
        END DO
      END DO

      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          shflx (i,j) = SHeat(i,j)
          lhflx (i,j) = LHeat(i,j)
          tsk         = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          lwup        = -shr_const_stebol * tsk**4

# if defined LONGWAVE
          tsk = t(i,j,N(ng),nrhs,itemp) + 273.16_r8
          tak = Tair(i,j)               + 273.16_r8

          cff3 = (0.7859_r8+0.03477_r8*Tair(i,j)) /                     &
     &           (1.0_r8   +0.00412_r8*Tair(i,j))
          cff3 = 10.0_r8**cff3   ! saturation vapor pre (hPa or mbar)
          cff3 = cff3*Hair(i,j)  ! water vapor pressure (hPa or mbar)

          cff2 = tak*tak*tak
          cff1 = cff2*tak
          lrflx(i,j) = -shr_const_stebol*                               &
     &                 (cff1*(0.39_r8-0.05_r8*SQRT(cff3))  *            &
     &                       (1.0_r8-0.6823_r8*cloud(i,j)*cloud(i,j))+  &
     &                  cff2*4.0_r8*(tsk-tak))
          lrflx(i,j) = lrflx(i,j)*Hscale_inv
# elif defined LONGWAVE_OUT
          lrflx(i,j) = lwup*Hscale_inv + lrflx(i,j)
# endif
          lhflx(i,j) = lhflx(i,j)*Hscale_inv
          shflx(i,j) = shflx(i,j)*Hscale_inv
          stflux(i,j,itemp) = (srflx(i,j)+lrflx(i,j)+                   &
     &                         lhflx(i,j)+shflx(i,j))
# ifdef MASKING
          stflux(i,j,itemp)=stflux(i,j,itemp)*rmask(i,j)
# endif
# ifdef WET_DRY
          stflux(i,j,itemp)=stflux(i,j,itemp)*rmask_wet(i,j)
# endif

!
!--- water flux ---
!     evap(i,j) = lhflx(i,j)/shr_const_latvap
!
!------------------------------------------------------------
! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
!------------------------------------------------------------
!     hol       = hol*ztref/blk_ZQ(ng)
!     xsq       = max( 1.0, sqrt(abs(1.0-16.0*hol)) )
!     xqq       = sqrt(xsq)
!     psix2     = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
!     fac       = (rh/shr_const_karman) * (alz + al2 - psixh + psix2 )
!     tref(i,j) = TairK - delt*fac
!     tref(i,j) = tref(i,j) - 0.01*ztref -- pot temp to temp correction
!     fac       = (re/shr_const_karman) * (alz + al2 - psixh + psix2 )
!     qref(i,j) = Qair(i,j) - delq*fac
!
!     duu10n(i,j) = u10n*u10n ! 10m wind speed squared
!
!     real(r8) :: al2     ! ln(zref/ztref)
!     al2 = LOG(zref/ztref)
!
# ifdef EMINUSP
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
          Hlv = (2.501_r8-0.00237_r8*t(i,j,N(ng),nrhs,itemp))*1.0e+6_r8
          evap(i,j) = -LHeat(i,j)/Hlv
#  ifdef EVAP_SHALLOW
          cffe = 1._r8/(h(i,j)+Zt_avg1(i,j))
          cffe = MIN(cffe,0.2_r8)
          IF (evap(i,j).gt.0._r8) THEN
            evap(i,j) = (1.0_r8+coef_shallow*cffe)*evap(i,j)
          END IF
#  endif
#  ifdef MASKING
          evap(i,j)=evap(i,j)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          evap(i,j)=evap(i,j)*rmask_wet(i,j)
#  endif
#  if defined ICE_MODEL && defined ICE_SNOWFALL
          stflux(i,j,isalt)=cff*(evap(i,j)-(rain(i,j) + snow(i,j)))
#  else
          stflux(i,j,isalt)=cff*(evap(i,j)-rain(i,j))
#  endif
#  ifdef MASKING
          stflux(i,j,isalt)=stflux(i,j,isalt)*rmask(i,j)
#  endif
#  ifdef WET_DRY
          stflux(i,j,isalt)=stflux(i,j,isalt)*rmask_wet(i,j)
#  endif
# endif
        END DO
      END DO

# ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Save some model states for later use by ice model
!-----------------------------------------------------------------------
!
      cff = 1.0_r8/rhosnow_dry(ng)
      DO j=JstrR,JendR
        DO i=Istr,IendR
!
!  Net heat flux from ocean to atmosphere
!    Note that in Mellor and Kantha qai is positive **upwards**,
!    while in ocean model heat flux is positive **downwards**.
!
#  ifdef ICE_DIAGS
          qswl(i,j) = srflx(i,j)*Hscale
          qlwl(i,j) = lrflx(i,j)*Hscale
#  endif
!
!  Snow fall rate in m/s
!
          IF (Tair(i,j).ge.0.0) THEN
            wsf(i,j) = 0.0_r8
          ELSE
#  ifdef ICE_SNOWFALL
            wsf(i,j) = snow(i,j)*cff
#  else
            wsf(i,j) = rain(i,j)*cff
#  endif
          END IF
        END DO
      END DO
# endif
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      cff=1._r8/rho0
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=Taux(i,j)*cff
# ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
# endif
# ifdef WET_DRY
          sustr(i,j)=sustr(i,j)*umask_wet(i,j)
# endif
# ifdef ICE_MODEL
          sustr_aw(i,j)=sustr(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=Tauy(i,j)*cff
# ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
# endif
# ifdef WET_DRY
          svstr(i,j)=svstr(i,j)*vmask_wet(i,j)
# endif
# ifdef ICE_MODEL
          svstr_aw(i,j)=svstr(i,j)
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
     &                          lrflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          lhflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          shflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflux(:,:,itemp))
# ifdef EMINUSP
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          evap)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflux(:,:,isalt))
# endif
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
# ifdef ICE_MODEL
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wsf)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr_aw)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr_aw)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qswi)
#  ifdef LONGWAVE_OUT
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          qlwi)
#  endif
# endif
      END IF

# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 4,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    lrflx, lhflx, shflx,                          &
     &                    stflux(:,:,itemp))
#  ifdef EMINUSP
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    evap,                                         &
     &                    stflux(:,:,isalt))
#  endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
#  ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    wsf, sustr_aw, svstr_aw)
#   ifdef LONGWAVE_OUT
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx, qswi, qlwi)
#   else
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx, qswi)
#   endif
#  endif
# endif
!
      RETURN
      END SUBROUTINE ccsm_flux_tile
#endif
      END module ccsm_flux_mod
