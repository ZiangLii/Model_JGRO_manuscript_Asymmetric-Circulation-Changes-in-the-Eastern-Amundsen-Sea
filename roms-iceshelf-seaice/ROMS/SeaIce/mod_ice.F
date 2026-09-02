#include "cppdefs.h"
      MODULE mod_ice
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Ice basic states.                                                   !
!                                                                      !
!  ui, vi     - [m/s] Ice velocities.                                  !
!  uie, vie   - [m/s] Ice velocities in EVP steps.                     !
!  ai         - [non] Sea ice concentration.                           !
!  hi         - [m]   sea ice thickness, averaged over whole grid cell.!
!  hsn        - [m]   snow thickness, averaged over whole grid cell.   !
!  ageice     - [day] Sea ice age.                                     !
!                                                                      !
!  Mass fluxes.                                                        !
!                                                                      !
!  iomflx - [m/s] net ice mass flux rate at ocean/ice boundary.        !
!  wfr    - [m/s] frazil ice growth rate.                              !
!  wai    - [m/s] ice surface melt rate.                               !
!  wao    - [m/s] ice growth rate at sea surface.                      !
!  wio    - [m/s] ice growth/melt rate at ice base.                    !
!  wro    - [m/s] ice surface runoff rate.                             !
!                                                                      !
!  Ice thermodynamic states.                                           !
!                                                                      !
!  tis      - [degC]   Snow top temperature (t3 im Mellor & Kantha)    !
!  ti       - [degC]   Ice interior temperature (t1 im Mellor & Kantha)!
!  t2       - [degC]   Ice/snow interface temperature.                 !
!  enthalpi - [degC m] Scaled perturbation ice heat content.           !
!  hage     - [yr m]   Ice age times thickness.                        !
!                                                                      !
!  Stress calculation.                                                 !
!                                                                      !
!  utau_iw - [m2/s2]   Ice/water friction velocity.                    !
!  chu_iw  - [m/s]     Ice/water stress velocity.                      !
!  spd_iw  - [m's]     Ice/water relative velocity.                    !
!  pice    - [Pa]      Ice pressure.                                   !
!  zetai   - [m2/s]    Ice bulk viscosity.                             !
!  eta     - [m2/s]    Ice shear viscosity.                            !
!  tauaiu  - [kg/m/s2] Ice/water stress in XI-direction.               !
!  tauaiv  - [kg/m/s2] Ice/water stress in ETA-direction.              !
!                                                                      !
!  Heat budgets.                                                       !
!                                                                      !
!  coef_ice_heat - Gradient coefficient for heat conductivity.         !
!  rhs_ice_heat  - RHS of ice heat budget equations.                   !
!                                                                      !
!  Other variables.                                                    !
!                                                                      !
!  Uwater, Vwater      - [m/s]       Ocean surface velocity.           !
!  sealev              - [m]         Sea surface height.               !
!  s0mk, t0mk          - [PSU, degC] Ice/ocean boundary S/T.           !
!  sig11, sig12, sig13 -             Ice stress terms.                 !
#ifdef ICE_STRENGTH_QUAD
!  pstar_grid          - Grid dependent pStar values.                  !
#endif
#ifdef ICE_LANDFAST
!  h_loadu, h_loadv    - term (h_u - h_cu)*exp(-C(1-A_u))              !
#endif
!                                                                      !
!=======================================================================
!
#ifdef ICE_MODEL
        USE mod_kinds
!
        implicit none
!
        PUBLIC :: allocate_ice
        PUBLIC :: deallocate_ice
        PUBLIC :: initialize_ice
!
!-----------------------------------------------------------------------
!  Define ICE structure.
!-----------------------------------------------------------------------
!
        TYPE T_ICE
!
!  Ice model state.
!
          real(r8), pointer :: ui(:,:,:)
          real(r8), pointer :: vi(:,:,:)
          real(r8), pointer :: uie(:,:,:)
          real(r8), pointer :: vie(:,:,:)
          real(r8), pointer :: ai(:,:,:)
          real(r8), pointer :: hi(:,:,:)
          real(r8), pointer :: hsn(:,:,:)
          real(r8), pointer :: ageice(:,:,:)
          real(r8), pointer :: wdiv(:,:)
!
          real(r8), pointer :: iomflx(:,:)
          real(r8), pointer :: wfr(:,:)
          real(r8), pointer :: wai(:,:)
          real(r8), pointer :: wao(:,:)
          real(r8), pointer :: wio(:,:)
          real(r8), pointer :: wro(:,:)
!
          real(r8), pointer :: tis(:,:)
          real(r8), pointer :: ti(:,:,:)
          real(r8), pointer :: t2(:,:)
          real(r8), pointer :: enthalpi(:,:,:)
          real(r8), pointer :: hage(:,:,:)
          real(r8), pointer :: utau_iw(:,:)
          real(r8), pointer :: chu_iw(:,:)
          real(r8), pointer :: spd_iw(:,:)
          real(r8), pointer :: pice(:,:)
          real(r8), pointer :: zetai(:,:)
          real(r8), pointer :: eta(:,:)
          real(r8), pointer :: tauaiu(:,:)
          real(r8), pointer :: tauaiv(:,:)
!
          real(r8), pointer :: coef_ice_heat(:,:)
          real(r8), pointer :: rhs_ice_heat(:,:)
          real(r8), pointer :: rfaci(:,:)
!
          real(r8), pointer :: uwater(:,:)
          real(r8), pointer :: vwater(:,:)
          real(r8), pointer :: sealev(:,:)
          real(r8), pointer :: s0mk(:,:)
          real(r8), pointer :: t0mk(:,:)
          real(r8), pointer :: sig11(:,:,:)
          real(r8), pointer :: sig22(:,:,:)
          real(r8), pointer :: sig12(:,:,:)
# ifdef ICE_STRENGTH_QUAD
          real(r8), pointer :: pstar_grid(:,:)
# endif
# ifdef ICE_LANDFAST
          real(r8), pointer :: h_loadu(:,:)
          real(r8), pointer :: h_loadv(:,:)
# endif

        END TYPE T_ICE

        TYPE (T_ICE), allocatable :: ICE(:)

      CONTAINS

      SUBROUTINE allocate_ice (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( ICE(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Ice model state
!
      allocate ( ICE(ng) % ui(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % vi(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % uie(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % vie(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % ai(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % hi(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % hsn(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % ageice(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % wdiv(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % iomflx(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % wfr(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % wai(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % wao(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % wio(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % wro(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % tis(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % ti(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % t2(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % enthalpi(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % hage(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d*2
      allocate ( ICE(ng) % utau_iw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % chu_iw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % spd_iw(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % pice(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % zetai(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % eta(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % tauaiu(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % tauaiv(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % coef_ice_heat(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % rhs_ice_heat(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % rfaci(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % uwater(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % vwater(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % sealev(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % s0mk(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % t0mk(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % sig11(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % sig22(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % sig12(LBi:UBi,LBj:UBj,2) )
      Dmem(ng)=Dmem(ng)+size2d
# ifdef ICE_STRENGTH_QUAD
      allocate ( ICE(ng) % pstar_grid(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif
# ifdef ICE_LANDFAST
      allocate ( ICE(ng) % h_loadu(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
      allocate ( ICE(ng) % h_loadv(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d
# endif

      RETURN
      END SUBROUTINE allocate_ice
!
      SUBROUTINE deallocate_ice (ng)
!
!=======================================================================
!                                                                      !
!  This routine deallocates all variables in the module for all nested !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param,   ONLY : Ngrids
# ifdef SUBOBJECT_DEALLOCATION
      USE destroy_mod, ONLY : destroy
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", deallocate_forces"
# ifdef SUBOBJECT_DEALLOCATION
!
!-----------------------------------------------------------------------
!  Deallocate each variable in the derived-type ICE structure
!  separately.
!-----------------------------------------------------------------------
!
!  Ice model state
!
      IF (.not.destroy(ng, ICE(ng)%ui, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%ui')) RETURN

      IF (.not.destroy(ng, ICE(ng)%vi, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%vi')) RETURN

      IF (.not.destroy(ng, ICE(ng)%uie, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%uie')) RETURN

      IF (.not.destroy(ng, ICE(ng)%vie, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%vie')) RETURN

      IF (.not.destroy(ng, ICE(ng)%ai, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%ai')) RETURN

      IF (.not.destroy(ng, ICE(ng)%hi, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%yi')) RETURN

      IF (.not.destroy(ng, ICE(ng)%hsn, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%hsn')) RETURN

      IF (.not.destroy(ng, ICE(ng)%ageice, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%ageice')) RETURN

      IF (.not.destroy(ng, ICE(ng)%wdiv, MyFile,                        &
     &                 __LINE__, 'ICE(ng)%wdiv')) RETURN

      IF (.not.destroy(ng, ICE(ng)%iomflx, MyFile,                    &
     &                 __LINE__, 'ICE(ng)%iomflx')) RETURN

      IF (.not.destroy(ng, ICE(ng)%wfr, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%wfr')) RETURN

      IF (.not.destroy(ng, ICE(ng)%wai, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%wai')) RETURN

      IF (.not.destroy(ng, ICE(ng)%wao, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%wao')) RETURN

      IF (.not.destroy(ng, ICE(ng)%wio, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%wio')) RETURN

      IF (.not.destroy(ng, ICE(ng)%wro, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%wro')) RETURN

      IF (.not.destroy(ng, ICE(ng)%tis, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%tis')) RETURN

      IF (.not.destroy(ng, ICE(ng)%ti, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%ti')) RETURN

      IF (.not.destroy(ng, ICE(ng)%t2, MyFile,                          &
     &                 __LINE__, 'ICE(ng)%t2')) RETURN

      IF (.not.destroy(ng, ICE(ng)%enthalpi, MyFile,                    &
     &                 __LINE__, 'ICE(ng)%enthalpi')) RETURN

      IF (.not.destroy(ng, ICE(ng)%hage, MyFile,                        &
     &                 __LINE__, 'ICE(ng)%hage')) RETURN

      IF (.not.destroy(ng, ICE(ng)%utau_iw, MyFile,                     &
     &                 __LINE__, 'ICE(ng)%utau_iw')) RETURN

      IF (.not.destroy(ng, ICE(ng)%chu_iw, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%chu_iw')) RETURN

      IF (.not.destroy(ng, ICE(ng)%spd_iw, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%spd_iw')) RETURN

      IF (.not.destroy(ng, ICE(ng)%pice, MyFile,                        &
     &                 __LINE__, 'ICE(ng)%pice')) RETURN

      IF (.not.destroy(ng, ICE(ng)%zetai, MyFile,                       &
     &                 __LINE__, 'ICE(ng)%zetai')) RETURN

      IF (.not.destroy(ng, ICE(ng)%eta, MyFile,                         &
     &                 __LINE__, 'ICE(ng)%eta')) RETURN

      IF (.not.destroy(ng, ICE(ng)%tauaiu, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%tauaiu')) RETURN

      IF (.not.destroy(ng, ICE(ng)%tauaiv, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%tauaiv')) RETURN

      IF (.not.destroy(ng, ICE(ng)%coef_ice_heat, MyFile,               &
     &                 __LINE__, 'ICE(ng)%coef_ice_heat')) RETURN

      IF (.not.destroy(ng, ICE(ng)%rhs_ice_heat, MyFile,                &
     &                 __LINE__, 'ICE(ng)%rhs_ice_heat')) RETURN

      IF (.not.destroy(ng, ICE(ng)%rfaci, MyFile,                       &
     &                 __LINE__, 'ICE(ng)%rfaci')) RETURN

      IF (.not.destroy(ng, ICE(ng)%uwater, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%uwater')) RETURN

      IF (.not.destroy(ng, ICE(ng)%vwater, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%vwater')) RETURN

      IF (.not.destroy(ng, ICE(ng)%sealev, MyFile,                      &
     &                 __LINE__, 'ICE(ng)%sealev')) RETURN

      IF (.not.destroy(ng, ICE(ng)%s0mk, MyFile,                        &
     &                 __LINE__, 'ICE(ng)%s0mk')) RETURN

      IF (.not.destroy(ng, ICE(ng)%t0mk, MyFile,                        &
     &                 __LINE__, 'ICE(ng)%t0mk')) RETURN

      IF (.not.destroy(ng, ICE(ng)%sig11, MyFile,                       &
     &                 __LINE__, 'ICE(ng)%sig11')) RETURN

      IF (.not.destroy(ng, ICE(ng)%sig22, MyFile,                       &
     &                 __LINE__, 'ICE(ng)%sig22')) RETURN

      IF (.not.destroy(ng, ICE(ng)%sig12, MyFile,                       &
     &                 __LINE__, 'ICE(ng)%sig12')) RETURN
#  ifdef ICE_STRENGTH_QUAD

      IF (.not.destroy(ng, ICE(ng)%pstar_grid, MyFile,                  &
     &                 __LINE__, 'ICE(ng)%pstar_grid')) RETURN
#  endif

      IF (.not.destroy(ng, ICE(ng)%h_loadu, MyFile,                     &
     &                 __LINE__, 'ICE(ng)%h_loadu')) RETURN

      IF (.not.destroy(ng, ICE(ng)%h_loadv, MyFile,                     &
     &                 __LINE__, 'ICE(ng)%h_loadv')) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type ICE structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(ICE)) deallocate ( ICE )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_ice
!
      SUBROUTINE initialize_ice (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
!
      real(r8), parameter :: IniVal = 0.0_r8
!
# include "tile.h"
!
!  Set array initialization range.
!
#ifdef DISTRIBUTE
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Ice model state.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ICE(ng) % ui(i,j,1) = IniVal
          ICE(ng) % ui(i,j,2) = IniVal
          ICE(ng) % vi(i,j,1) = IniVal
          ICE(ng) % vi(i,j,2) = IniVal
          ICE(ng) % uie(i,j,1) = IniVal
          ICE(ng) % uie(i,j,2) = IniVal
          ICE(ng) % vie(i,j,1) = IniVal
          ICE(ng) % vie(i,j,2) = IniVal
          ICE(ng) % ai(i,j,1) = IniVal
          ICE(ng) % ai(i,j,2) = IniVal
          ICE(ng) % hi(i,j,1) = IniVal
          ICE(ng) % hi(i,j,2) = IniVal
          ICE(ng) % hsn(i,j,1) = IniVal
          ICE(ng) % hsn(i,j,2) = IniVal
          ICE(ng) % ageice(i,j,1) = IniVal
          ICE(ng) % ageice(i,j,2) = IniVal
          ICE(ng) % wdiv(i,j) = IniVal
          ICE(ng) % iomflx(i,j) = IniVal
          ICE(ng) % wfr(i,j) = IniVal
          ICE(ng) % wai(i,j) = IniVal
          ICE(ng) % wao(i,j) = IniVal
          ICE(ng) % wio(i,j) = IniVal
          ICE(ng) % wro(i,j) = IniVal
          ICE(ng) % tis(i,j) = IniVal
          ICE(ng) % ti(i,j,1) = IniVal
          ICE(ng) % ti(i,j,2) = IniVal
          ICE(ng) % t2(i,j) = IniVal
          ICE(ng) % enthalpi(i,j,1) = IniVal
          ICE(ng) % enthalpi(i,j,2) = IniVal
          ICE(ng) % hage(i,j,1) = IniVal
          ICE(ng) % hage(i,j,2) = IniVal
          ICE(ng) % utau_iw(i,j) = IniVal
          ICE(ng) % chu_iw(i,j) = IniVal
          ICE(ng) % spd_iw(i,j) = IniVal
          ICE(ng) % pice(i,j) = IniVal
          ICE(ng) % eta(i,j) = IniVal
          ICE(ng) % tauaiu(i,j) = IniVal
          ICE(ng) % tauaiv(i,j) = IniVal
          ICE(ng) % coef_ice_heat(i,j) = IniVal
          ICE(ng) % rhs_ice_heat(i,j) = IniVal
          ICE(ng) % rfaci(i,j) = 1.0_r8
          ICE(ng) % zetai(i,j) = IniVal
          ICE(ng) % uwater(i,j) = IniVal
          ICE(ng) % vwater(i,j) = IniVal
          ICE(ng) % sealev(i,j) = IniVal
          ICE(ng) % s0mk(i,j) = IniVal
          ICE(ng) % t0mk(i,j) = IniVal
          ICE(ng) % sig11(i,j,1) = IniVal
          ICE(ng) % sig11(i,j,2) = IniVal
          ICE(ng) % sig22(i,j,1) = IniVal
          ICE(ng) % sig22(i,j,2) = IniVal
          ICE(ng) % sig12(i,j,1) = IniVal
          ICE(ng) % sig12(i,j,2) = IniVal
# ifdef ICE_STRENGTH_QUAD
          ICE(ng) % pstar_grid(i,j) = IniVal
# endif
# ifdef ICE_LANDFAST
          ICE(ng) % h_loadu(i,j) = IniVal
          ICE(ng) % h_loadv(i,j) = IniVal
# endif

        END DO
      END DO

      RETURN
      END SUBROUTINE initialize_ice
#endif
      END MODULE mod_ice
