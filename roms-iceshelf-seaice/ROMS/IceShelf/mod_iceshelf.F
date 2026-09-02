#include "cppdefs.h"
      MODULE mod_iceshelf
#ifdef ICESHELF 
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!========================================== Benjamin K. Galton-Fenzi ===
!                                                                      !
!  Ice Shelf Model Kernel Variables:                                   !
!  This module declares Ice Shelf Model Kernel constants and structure !
!  containing several variables describing ice draft and thickness     !
!  properties.                                                         !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
!
        implicit none

        PUBLIC :: allocate_iceshelfvar
        PUBLIC :: deallocate_iceshelfvar
        PUBLIC :: initialize_iceshelfvar
!
!======================================================================!
!                                                                      !
! Model parameters                                                     !
!                                                                      !
! a, b, c - freezing point EOS constants [degC PSU^-1, degC, degC m^-1]!
! L       - latent heat of melting [J kg^-1]                           !
! gamma   - 2-EQ formula turbulent exchange coefficient                !
! refSalt - reference salinity [PSU]                                   !
!                                                                      !
! Ti, Si  - ice temperature/salinity                                   !
!                                                                      !
!======================================================================!
!
        real(r8), parameter :: a = -0.057_r8
        real(r8), parameter :: b = 0.0939_r8
        real(r8), parameter :: c = -7.61e-4_r8
!
# if defined ICESHELF_2EQN_VBC
        real(r8), parameter :: L = 3.34e5_r8
        real(r8), parameter :: gamma = 0.0001_r8
        real(r8), parameter :: refSalt = 34.4_r8
# elif defined ICESHELF_3EQN_VBC
        real(r8), parameter :: Pr = 13.8_r8
        real(r8), parameter :: Sc = 2432.2_r8
        real(r8), parameter :: Cd = 5.0e-3_r8
        real(r8), parameter :: visc = 1.95e-6_r8
        real(r8), parameter :: L = 3.33e5_r8
!       real(r8), parameter :: small = 1.0e-3_r8 ! Diffusion + tidal
        real(r8), parameter :: small = 8.0e-6_r8 ! Diffusion only
        real(r8), parameter :: dt_i = 1.54e-6_r8
        real(r8), parameter :: cp_w = 3947.0_r8
        real(r8), parameter :: rho_i = 920.0_r8
!
        real(r8), parameter :: Ti = -20.0_r8
        real(r8), parameter :: Si = 0.0_r8
# endif
        real(r8), parameter :: eps = 1.0E-14_r8
        real(r8), parameter :: trcIce = 1.0_r8
!
! =====================================================================!
!                                                                      !
! Model parameters                                                     !
!                                                                      !
! m - iceshelf melt rate, wate equvilent [m s-1]                       !
! gammaT, gammaS - heat/salt turbulent exchange coefficient [m s-1]    !
! Tstar - temperature gradient from water to boundary layer            !
! Ustar - friction velocity [m s^-1]                                   !
! Tb, Sb - boundary layer temperature/salinity                         !
! iceshelf_draft - draft depth of iceshelf [m]                         !
!                                                                      !
! =====================================================================!
!
        TYPE T_ICESHELFVAR
!
!  Nonlinear model state.
!
          real(r8), pointer :: m(:,:)
# ifdef ICESHELF_3EQN_VBC
          real(r8), pointer :: gammaT(:,:)
          real(r8), pointer :: gammaS(:,:)
          real(r8), pointer :: Tstar(:,:)
          real(r8), pointer :: Ustar(:,:)
          real(r8), pointer :: Tb(:,:)
          real(r8), pointer :: Sb(:,:)
# endif
# ifdef ICESHELF_SOF
          real(r8), pointer :: SOFs(:,:)
          real(r8), pointer :: SOFh(:,:)
# endif
# ifdef ICESHELF_MORPH
          real(r8), pointer :: iceshelf_draft0(:,:)
          real(r8), pointer :: iceshelf_draft(:,:,:)
#  ifdef FISOC_DDDT
          real(r8), pointer :: iceshelf_dddt(:,:)
#  endif
# endif
# ifdef ICESHELF_TRACER
          real(r8), pointer :: istrc(:,:,:)
          real(r8), pointer :: isvol(:,:)
# endif
        END TYPE T_ICESHELFVAR

        TYPE (T_ICESHELFVAR), allocatable :: ICESHELFVAR(:)
!
      CONTAINS
!
      SUBROUTINE allocate_iceshelfvar (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      real(r8) :: size2d
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( ICESHELFVAR(Ngrids) )
!
!  Set horizontal array size.
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Nonlinear model state.
!
      allocate ( ICESHELFVAR(ng) % m(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# ifdef ICESHELF_3EQN_VBC
      allocate ( ICESHELFVAR(ng) % gammaT(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( ICESHELFVAR(ng) % gammaS(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( ICESHELFVAR(ng) % Tstar(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( ICESHELFVAR(ng) % Ustar(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( ICESHELFVAR(ng) % Tb(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( ICESHELFVAR(ng) % Sb(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# endif
# ifdef ICESHELF_SOF
      allocate ( ICESHELFVAR(ng) % SOFs(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( ICESHELFVAR(ng) % SOFh(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# endif
# ifdef ICESHELF_MORPH
      allocate ( ICESHELFVAR(ng) % iceshelf_draft0(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

      allocate ( ICESHELFVAR(ng) % iceshelf_draft(LBi:UBi,LBj:UBj,1:2) )
      Dmem(ng)=Dmem(ng)+2.0_r8*size2d

#  ifdef FISOC_DDDT
      allocate ( ICESHELFVAR(ng) % iceshelf_dddt(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

#  endif
# endif
# ifdef ICESHELF_TRACER
      allocate ( ICESHELFVAR(ng) % istrc(LBi:UBi,LBj:UBj,NAT+1:NT(ng)) )
      Dmem(ng)=Dmem(ng)+REAL(NT(ng)-2)*size2d

      allocate ( ICESHELFVAR(ng) % isvol(LBi:UBi,LBj:UBj) )
      Dmem(ng)=Dmem(ng)+size2d

# endif
!
      RETURN
      END SUBROUTINE allocate_iceshelfvar
!
      SUBROUTINE deallocate_iceshelfvar (ng)
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
!  Deallocate each variable in the derived-type ICESHELFVAR structure
!  separately.
!-----------------------------------------------------------------------
!
!  Iceshelf model state
!
      IF (.not.destroy(ng, ICESHELFVAR(ng)%m, MyFile,                   &
     &                 __LINE__, 'ICESHELFVAR(ng)%m')) RETURN
#  ifdef ICESHELF_3EQN_VBC

      IF (.not.destroy(ng, ICESHELFVAR(ng)%gammaT, MyFile,              &
     &                 __LINE__, 'ICESHELFVAR(ng)%gammaT')) RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%gammaS, MyFile,              &
     &                 __LINE__, 'ICESHELFVAR(ng)%gammaS')) RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%Tstar, MyFile,               &
     &                 __LINE__, 'ICESHELFVAR(ng)%Tstar')) RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%Ustar, MyFile,               &
     &                 __LINE__, 'ICESHELFVAR(ng)%Ustar')) RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%Tb, MyFile,                  &
     &                 __LINE__, 'ICESHELFVAR(ng)%Tb')) RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%Sb, MyFile,                  &
     &                 __LINE__, 'ICESHELFVAR(ng)%Sb')) RETURN
#  endif
#  ifdef ICESHELF_SOF

      IF (.not.destroy(ng, ICESHELFVAR(ng)%SOFs, MyFile,                &
     &                 __LINE__, 'ICESHELFVAR(ng)%SOFh')) RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%SOFs, MyFile,                &
     &                 __LINE__, 'ICESHELFVAR(ng)%SOFh')) RETURN
#  endif
#  ifdef ICESHELF_MORPH

      IF (.not.destroy(ng, ICESHELFVAR(ng)%iceshelf_draft0, MyFile,     &
     &                 __LINE__, 'ICESHELFVAR(ng)%iceshelf_draft0'))    &
     &  RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%iceshelf_draft, MyFile,      &
     &                 __LINE__, 'ICESHELFVAR(ng)%iceshelf_draft'))     &
     &  RETURN
#   ifdef FISOC_DDDT

      IF (.not.destroy(ng, ICESHELFVAR(ng)%iceshelf_dddt, MyFile,       &
     &                 __LINE__, 'ICESHELFVAR(ng)%iceshelf_dddt'))      &
     &  RETURN
#   endif
#  endif
#  ifdef ICESHELF_TRACER

      IF (.not.destroy(ng, ICESHELFVAR(ng)%istrc, MyFile,               &
     &                 __LINE__, 'ICESHELFVAR(ng)%istrc'))              &
     &  RETURN

      IF (.not.destroy(ng, ICESHELFVAR(ng)%isvol, MyFile,               &
     &                 __LINE__, 'ICESHELFVAR(ng)%isvol'))              &
     &  RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Deallocate derived-type ICESHELFVAR structure.
!-----------------------------------------------------------------------
!
      IF (ng.eq.Ngrids) THEN
        IF (allocated(ICESHELFVAR)) deallocate ( ICESHELFVAR )
      END IF
!
      RETURN
      END SUBROUTINE deallocate_iceshelfvar
!
!
      SUBROUTINE initialize_iceshelfvar (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the "shared     !
!  arrays" across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, itrc

      real(r8), parameter :: IniVal = 0.0_r8

# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# else
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
# endif
!
!-----------------------------------------------------------------------
!  Initialize iceshelf structure variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ICESHELFVAR(ng) % m(i,j) = IniVal
# ifdef ICESHELF_3EQN_VBC 
            ICESHELFVAR(ng) % gammaT(i,j) = IniVal
            ICESHELFVAR(ng) % gammaS(i,j) = IniVal
            ICESHELFVAR(ng) % Tstar(i,j) = IniVal
            ICESHELFVAR(ng) % Ustar(i,j) = IniVal
            ICESHELFVAR(ng) % Tb(i,j) = IniVal
            ICESHELFVAR(ng) % Sb(i,j) = IniVal
# endif
# ifdef ICESHELF_SOF
            ICESHELFVAR(ng) % SOFs(i,j) = IniVal
            ICESHELFVAR(ng) % SOFh(i,j) = IniVal
# endif
# ifdef ICESHELF_MORPH
            ICESHELFVAR(ng) % iceshelf_draft0(i,j) = IniVal
            ICESHELFVAR(ng) % iceshelf_draft(i,j,1) = IniVal
            ICESHELFVAR(ng) % iceshelf_draft(i,j,2) = IniVal
#  ifdef FISOC_DDDT
            ICESHELFVAR(ng) % iceshelf_dddt(i,j) = IniVal
#  endif
# endif
# ifdef ICESHELF_TRACER
            DO itrc=NAT+1,NT(ng)
              ICESHELFVAR(ng) % istrc(i,j,itrc) = IniVal
            END DO
            ICESHELFVAR(ng) % isvol(i,j) = IniVal
# endif
          END DO
        END DO
      END IF
!
      RETURN
      END SUBROUTINE initialize_iceshelfvar
#endif
      END MODULE mod_iceshelf
