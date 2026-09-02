#include "cppdefs.h"
      MODULE uibc_mod
#ifdef ICE_MODEL
! 
!=======================================================================
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license               Chuning Wang   !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine computes lateral boundary conditions for U-velocity.!
!                                                                      !
!=======================================================================
!

      implicit none

      PRIVATE
      PUBLIC  :: uibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE uibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile
!
#include "tile.h"
!
      CALL  uibc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 liuol(ng), liunw(ng),                            &
     &                 ICE(ng) % ui)
      RETURN
      END SUBROUTINE uibc
!
!***********************************************************************
      SUBROUTINE uibc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liunw,                               &
     &                      ui)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_boundary
      USE mod_grid
      USE mod_scalars

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liuol, liunw

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, Imax, Imin, j
      real(r8) :: Ce, Cx
      real(r8) :: cff, dUde, dUdt, dUdx, phi, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isUice,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=ui(Istr  ,j  ,liuol)-ui(Istr  ,j-1,liuol)
            grad(Istr+1,j)=ui(Istr+1,j  ,liuol)-ui(Istr+1,j-1,liuol)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dUdt=ui(Istr+1,j,liuol)-ui(Istr+1,j,liunw)
              dUdx=ui(Istr+1,j,liunw)-ui(Istr+2,j,liunw)

              IF (LBC(iwest,isUice,ng)%nudging) THEN
                IF ((dUdt*dUdx).lt.0.0_r8) THEN
                  tau=M2obc_in (ng,iwest)
                ELSE
                  tau=M2obc_out(ng,iwest)
                END IF
# ifdef IMPLICIT_NUDGING
                IF (tau.gt.0.0_r8) tau=1.0_r8/tau
# else
                tau=tau*dt(ng)
# endif
              END IF

              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Istr+1,j  )+                              &
     &                   grad(Istr+1,j+1))).gt.0.0_r8) THEN
                dUde=grad(Istr+1,j  )
              ELSE
                dUde=grad(Istr+1,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,tol)
              Cx=dUdt*dUdx
# ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dUdt*dUde,-cff))
# else
              Ce=0.0_r8
# endif
              ui(Istr,j,liunw)=(cff*ui(Istr  ,j,liuol)+                 &
     &                          Cx *ui(Istr+1,j,liunw)-                 &
     &                          MAX(Ce,0.0_r8)*grad(Istr,j  )-          &
     &                          MIN(Ce,0.0_r8)*grad(Istr,j+1)) /        &
     &                         (cff+Cx)

              IF (LBC(iwest,isUice,ng)%nudging) THEN
# ifdef IMPLICIT_NUDGING
                phi=dt(ng)/(tau+dt(ng))
                ui(Istr,j,liunw)=(1.0_r8-phi)*ui(Istr,j,liunw)+         &
     &                           phi*BOUNDARY(ng)%ui_west(j)
# else
                ui(Istr,j,liunw)=ui(Istr,j,liunw)+                      &
     &                           tau*(BOUNDARY(ng)%ui_west(j)-          &
     &                           ui(Istr,j,liuol))
# endif
              END IF
# ifdef MASKING
              ui(Istr,j,liunw)=ui(Istr,j,liunw)*                        &
     &                         GRID(ng)%umask(Istr,j)
# endif
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isUice,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ui(Istr,j,liunw)=BOUNDARY(ng)%ui_west(j)
# ifdef MASKING
              ui(Istr,j,liunw)=ui(Istr,j,liunw)*                        &
     &                         GRID(ng)%umask(Istr,j)
# endif
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isUice,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ui(Istr,j,liunw)=ui(Istr+1,j,liunw)
# ifdef MASKING
              ui(Istr,j,liunw)=ui(Istr,j,liunw)*                        &
     &                         GRID(ng)%umask(Istr,j)
# endif
            END IF
          END DO
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(iwest,isUice,ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              IF (ui(Istr+1,j,liunw).ge.0._r8) THEN
                ui(Istr,j,liunw)=ui(Istr+1,j,liunw)
              ELSE
                ui(Istr,j,liunw)=BOUNDARY(ng)%ui_west(j)
              END IF
# ifdef MASKING
              ui(Istr,j,liunw)=ui(Istr,j,liunw)*                        &
     &                         GRID(ng)%umask(Istr,j)
# endif
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isUice,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ui(Istr,j,liunw)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isUice,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=ui(Iend  ,j  ,liuol)-ui(Iend  ,j-1,liuol)
            grad(Iend+1,j)=ui(Iend+1,j  ,liuol)-ui(Iend+1,j-1,liuol)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dUdt=ui(Iend,j,liuol)-ui(Iend  ,j,liunw)
              dUdx=ui(Iend,j,liunw)-ui(Iend-1,j,liunw)

              IF (LBC(ieast,isUice,ng)%nudging) THEN
                IF ((dUdt*dUdx).lt.0.0_r8) THEN
                  tau=M2obc_in (ng,ieast)
                ELSE
                  tau=M2obc_out(ng,ieast)
                END IF
# ifdef IMPLICIT_NUDGING
                IF (tau.gt.0.0_r8) tau=1.0_r8/tau
# else
                tau=tau*dt(ng)
# endif
              END IF

              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Iend,j  )+                                &
     &                   grad(Iend,j+1))).gt.0.0_r8) THEN
                dUde=grad(Iend,j)
              ELSE
                dUde=grad(Iend,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,tol)
              Cx=dUdt*dUdx
# ifdef RADIATION_2D
              Ce=MIN(cff,MAX(dUdt*dUde,-cff))
# else
              Ce=0.0_r8
# endif
              ui(Iend+1,j,liunw)=(cff*ui(Iend+1,j,liuol)+               &
     &                            Cx *ui(Iend  ,j,liunw)-               &
     &                            MAX(Ce,0.0_r8)*grad(Iend+1,j  )-      &
     &                            MIN(Ce,0.0_r8)*grad(Iend+1,j+1)) /    &
     &                           (cff+Cx)

              IF (LBC(ieast,isUice,ng)%nudging) THEN
#ifdef IMPLICIT_NUDGING
                phi=dt(ng)/(tau+dt(ng))
                ui(Iend+1,j,liunw)=(1.0_r8-phi)*ui(Iend+1,j,liunw)+     &
     &                             phi*BOUNDARY(ng)%ui_east(j)
# else
                ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)+                  &
     &                             tau*(BOUNDARY(ng)%ui_east(j)-        &
     &                             ui(Iend+1,j,liuol))
# endif
              END IF
# ifdef MASKING
              ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                    &
     &                           GRID(ng)%umask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isUice,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ui(Iend+1,j,liunw)=BOUNDARY(ng)%ui_east(j)
# ifdef MASKING
              ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                    &
     &                           GRID(ng)%umask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isUice,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ui(Iend+1,j,liunw)=ui(Iend,j,liunw)
# ifdef MASKING
              ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                    &
     &                           GRID(ng)%umask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        ELSE IF (LBC(ieast,isUice,ng)%mixed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              IF (ui(Iend,j,liunw).le.0._r8) THEN
                ui(Iend+1,j,liunw)=BOUNDARY(ng)%ui_east(j)
              ELSE
                ui(Iend+1,j,liunw)=ui(Iend,j,liunw)
              END IF
# ifdef MASKING
              ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                    &
     &                           GRID(ng)%rmask(Iend+1,j)
# endif
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isUice,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ui(Iend+1,j,liunw)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isUice,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=ui(i+1,Jstr-1,liuol)-ui(i  ,Jstr-1,liuol)
            grad(i,Jstr  )=ui(i+1,Jstr  ,liuol)-ui(i  ,Jstr  ,liuol)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dUdt=ui(i,Jstr,liuol)-ui(i,Jstr  ,liunw)
              dUde=ui(i,Jstr,liunw)-ui(i,Jstr+1,liunw)

              IF (LBC(isouth,isUice,ng)%nudging) THEN
                IF ((dUdt*dUde).lt.0.0_r8) THEN
                  tau=M2obc_in (ng,isouth)
                ELSE
                  tau=M2obc_out(ng,isouth)
                END IF
# ifdef IMPLICIT_NUDGING
                IF (tau.gt.0.0_r8) tau=1.0_r8/tau
# else
                tau=tau*dt(ng)
# endif
              END IF

              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jstr)+                                &
     &                   grad(i  ,Jstr))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jstr)
              ELSE
                dUdx=grad(i  ,Jstr)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,tol)
# ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
# else
              Cx=0.0_r8
# endif
              Ce=dUdt*dUde
              ui(i,Jstr-1,liunw)=(cff*ui(i,Jstr-1,liuol)+               &
     &                            Ce *ui(i,Jstr  ,liunw)-               &
     &                            MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-      &
     &                            MIN(Cx,0.0_r8)*grad(i  ,Jstr-1)) /    &
     &                           (cff+Ce)

              IF (LBC(isouth,isUice,ng)%nudging) THEN
# ifdef IMPLICIT_NUDGING
                phi=dt(ng)/(tau+dt(ng))
                ui(i,Jstr-1,liunw)=(1.0_r8-phi)*ui(i,Jstr-1,liunw)+     &
     &                             phi*BOUNDARY(ng)%ui_south(i)
# else
                ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)+                  &
     &                             tau*(BOUNDARY(ng)%ui_south(i)-       &
     &                                  ui(i,Jstr-1,liuol))
# endif
              END IF
# ifdef MASKING
              ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isUice,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ui(i,Jstr-1,liunw)=BOUNDARY(ng)%ui_south(i)
# ifdef MASKING
              ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isUice,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ui(i,Jstr-1,liunw)=ui(i,Jstr,liunw)
# ifdef MASKING
              ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isUice,ng)%mixed) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ui(i,Jstr-1,liunw)=ui(i,Jstr,liunw)
# ifdef MASKING
              ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                    &
     &                           GRID(ng)%rmask(i,Jstr-1)
# endif
            END IF
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isUice,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%south(i)) THEN
              ui(i,Jstr-1,liunw)=gamma2(ng)*ui(i,Jstr,liunw)
# ifdef MASKING
              ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jstr-1)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isUice,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jend  )=ui(i+1,Jend  ,liuol)-ui(i  ,Jend  ,liuol)
            grad(i,Jend+1)=ui(i+1,Jend+1,liuol)-ui(i  ,Jend+1,liuol)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dUdt=ui(i,Jend,liuol)-ui(i,Jend  ,liunw)
              dUde=ui(i,Jend,liunw)-ui(i,Jend-1,liunw)

              IF (LBC(inorth,isUice,ng)%nudging) THEN
                IF ((dUdt*dUde).lt.0.0_r8) THEN
                  tau=M2obc_in (ng,inorth)
                ELSE
                  tau=M2obc_out(ng,inorth)
                END IF
# ifdef IMPLICIT_NUDGING
                IF (tau.gt.0.0_r8) tau=1.0_r8/tau
# else
                tau=tau*dt(ng)
# endif
              END IF

              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jend)+                                &
     &                   grad(i  ,Jend))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jend)
              ELSE
                dUdx=grad(i  ,Jend)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,tol)
# ifdef RADIATION_2D
              Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
# else
              Cx=0.0_r8
# endif
              Ce=dUdt*dUde
              ui(i,Jend+1,liunw)=(cff*ui(i,Jend+1,liuol)+               &
     &                            Ce *ui(i,Jend  ,liunw)-               &
     &                            MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-      &
     &                            MIN(Cx,0.0_r8)*grad(i  ,Jend+1)) /    &
     &                           (cff+Ce)
              IF (LBC(inorth,isUice,ng)%nudging) THEN
# ifdef IMPLICIT_NUDGING
                phi=dt(ng)/(tau+dt(ng))
                ui(i,Jend+1,liunw)=(1.0_r8-phi)*ui(i,Jend+1,liunw)+     &
     &                             phi*BOUNDARY(ng)%ui_north(i)
# else
                ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)+                  &
     &                             tau*(BOUNDARY(ng)%ui_north(i)-       &
     &                                  ui(i,Jend+1,liuol))
# endif
              END IF
# ifdef MASKING
              ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isUice,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ui(i,Jend+1,liunw)=BOUNDARY(ng)%ui_north(i)
# ifdef MASKING
              ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isUice,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ui(i,Jend+1,liunw)=ui(i,Jend,liunw)
# ifdef MASKING
              ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isUice,ng)%mixed) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ui(i,Jend+1,liunw)=ui(i,Jend,liunw)
# ifdef MASKING
              ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                    &
     &                           GRID(ng)%rmask(i,Jend+1)
# endif
            END IF
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isUice,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%north(i)) THEN
              ui(i,Jend+1,liunw)=gamma2(ng)*ui(i,Jend,liunw)
# ifdef MASKING
              ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                    &
     &                           GRID(ng)%umask(i,Jend+1)
# endif
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            ui(Istr  ,Jstr-1,liunw)=0.5_r8*(ui(Istr+1,Jstr-1,liunw)+    &
     &                                      ui(Istr  ,Jstr  ,liunw))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            ui(Iend+1,Jstr-1,liunw)=0.5_r8*(ui(Iend  ,Jstr-1,liunw)+    &
     &                                      ui(Iend+1,Jstr  ,liunw))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            ui(Istr  ,Jend+1,liunw)=0.5_r8*(ui(Istr  ,Jend  ,liunw)+    &
     &                                      ui(Istr+1,Jend+1,liunw))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            ui(Iend+1,Jend+1,liunw)=0.5_r8*(ui(Iend+1,Jend  ,liunw)+    &
     &                                      ui(Iend  ,Jend+1,liunw))
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE uibc_tile
#endif
      END MODULE uibc_mod
