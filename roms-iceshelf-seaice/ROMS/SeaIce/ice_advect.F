#include "cppdefs.h"
      MODULE ice_advect_mod
#ifdef ICE_ADVECT
!
!=================================================== W. Paul Budgell ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine performs advection of ice using one of 2 options:   !
!  (1) Smolarkiewicz second-order upwind scheme. Reference:            !
!      Smolarkiewicz and Grabowski (1990)                              !
!  (2) MPDATA advection scheme.                                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ice_advect

      CONTAINS

# if defined ICE_SMOLAR
#  include "ice_smolar.h"
# endif

#endif
      END MODULE ice_advect_mod
