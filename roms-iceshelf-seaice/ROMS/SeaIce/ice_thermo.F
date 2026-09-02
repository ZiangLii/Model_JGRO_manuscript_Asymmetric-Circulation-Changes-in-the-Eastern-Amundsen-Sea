#include "cppdefs.h"
      MODULE ice_thermo_mod
#ifdef ICE_THERMO
!
!=================================================== W. Paul Budgell ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine evaluates the ice thermodynamic growth and decay    !
!  term based on :                                                     !
!    Mellor and Kantha (1989) and Parkinson and Washington (1979)      !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ice_thermo

      CONTAINS

# ifdef ICE_MK
#  include "ice_mk.h"
# endif

#endif
      END MODULE ice_thermo_mod
