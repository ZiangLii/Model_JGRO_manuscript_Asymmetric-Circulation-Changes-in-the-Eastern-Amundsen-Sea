#include "cppdefs.h"
      SUBROUTINE seaice
#ifdef ICE_MODEL
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This is the main driver routine for the sea ice portion of the      !
!  model.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping

      USE mod_forces
      USE mod_ice

      USE ice_spdiw_mod,    ONLY : ice_spdiw
      USE ice_vbc_mod,      ONLY : ice_vbc
# ifdef ICE_MOMENTUM
      USE ice_evp_mod,      ONLY : ice_evp
      USE ice_evp_sig_mod,  ONLY : ice_evp_sig
      USE ice_elastic_mod,  ONLY : ice_elastic
# endif
# ifdef ICE_ADVECT
      USE ice_advect_mod,   ONLY : ice_advect
      USE ice_enthalpi_mod, ONLY : ice_enthalpi
# endif
# ifdef ICE_THERMO
      USE ice_thermo_mod,   ONLY : ice_thermo
# endif
# if defined ICE_ADVECT || defined ICE_THERMO
      USE ice_limit_mod,    ONLY : ice_limit
# endif
!
      implicit none
!
!  Local variable declarations.
!
      integer :: tile
      integer :: ig, ng, nl, my_ievp
      real(r8), parameter :: dt_large = 1.0E+23_r8
!
      NEST_LAYER : DO nl=1,NestLayers
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          liold(ng) = linew(ng)
          linew(ng) = 3 - liold(ng)
!
! ----------------------------------------------------------------------
!  Compute the ice-ocean shear.
! ----------------------------------------------------------------------
!
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_spdiw (ng, tile)
          END DO
        END DO
!$OMP BARRIER
!
! ----------------------------------------------------------------------
!  Compute the stresses on the ice from the air and water.
! ----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_vbc (ng, tile)
          END DO
!$OMP BARRIER
        END DO
# ifdef ICE_MOMENTUM
!
! ----------------------------------------------------------------------
!  Compute the internal ice stresses according to the 
!  Elastic-Viscous-Plastic rheology (EVP).
! ----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)

          liuol(ng) = liunw(ng)
          liunw(ng) = 3-liuol(ng)

          dte(ng) = dtice(ng)/FLOAT(nevp(ng))

          DO my_ievp=1,nevp(ng)

            lieol(ng) = lienw(ng)
            lienw(ng) = 3-lieol(ng)

            ievp(ng)=my_ievp
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_evp (ng, tile)
            END DO
!$OMP BARRIER
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_evp_sig (ng, tile)
            END DO
!$OMP BARRIER
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ice_elastic (ng, tile)
            END DO
!$OMP BARRIER
          END DO
        END DO
# endif
# ifdef ICE_ADVECT
!
! ----------------------------------------------------------------------
!  Compute the ice enthalpi before advection.
! ----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_enthalpi (ng, tile)
          END DO
!$OMP BARRIER
        END DO
!
! ----------------------------------------------------------------------
!  Compute the advection of the ice tracer fields.
! ----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_advect (ng, tile)
            CALL ice_limit (ng, tile)
          END DO
!$OMP BARRIER
        END DO
# endif
# ifdef ICE_THERMO
!
! ----------------------------------------------------------------------
!  Compute the ice thermodynamics.
! ----------------------------------------------------------------------
!
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ice_thermo (ng, tile)
            CALL ice_limit (ng, tile)
          END DO
!$OMP BARRIER
        END DO
# endif
      END DO NEST_LAYER
#endif
      RETURN
      END SUBROUTINE seaice
