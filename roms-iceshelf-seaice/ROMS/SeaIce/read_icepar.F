#include "cppdefs.h"
#if defined ICE_MODEL || defined ICESHELF
      SUBROUTINE read_IcePar (inp, out, Lwrite)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group          Chuning Wang   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports ice model input parameters.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, ng, status
      integer :: igrid, nline
      integer :: i, itrc

      logical, allocatable       :: Lswitch(:)
#if defined T_PASSIVE && ( defined ICESHELF_TRACER || \
                           defined ICE_FRZ_TRACER )
      logical, allocatable       :: LtracerI(:,:)
#endif

      real(r8)                   :: Rvalue(1)
      real(r8), dimension(nRval) :: Rval

      character (len=40)  :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Read in ice model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (.not.allocated(Lswitch)) THEN
            allocate ( Lswitch(Ngrids) )
          END IF
#if defined T_PASSIVE && ( defined ICESHELF_TRACER || \
                           defined ICE_FRZ_TRACER )
          IF (.not.allocated(LtracerI)) THEN
            allocate ( LtracerI(NPT,Ngrids) )
          END IF
#endif
          SELECT CASE (TRIM(KeyWord))
!
! Switches to turn ice model on/off of each nested grid.
!
# ifdef ICE_MODEL
            CASE ('Lice')
              Npts=load_l(Nval, Cval, Ngrids, Lice)
!
! Boundary options.
!
            CASE ('LBC(isUice)')
              Npts=load_lbc(Nval, Cval, line, nline, isUice, igrid,     &
     &                      0, 0, Vname(1,idIuice), LBC)
            CASE ('LBC(isVice)')
              Npts=load_lbc(Nval, Cval, line, nline, isVice, igrid,     &
     &                      0, 0, Vname(1,idIvice), LBC)
            CASE ('LBC(isAice)')
              Npts=load_lbc(Nval, Cval, line, nline, isAice, igrid,     &
     &                      0, 0, Vname(1,idIaice), LBC)
            CASE ('LBC(isHice)')
              Npts=load_lbc(Nval, Cval, line, nline, isHice, igrid,     &
     &                      0, 0, Vname(1,idIhice), LBC)
            CASE ('LBC(isTice)')
              Npts=load_lbc(Nval, Cval, line, nline, isTice, igrid,     &
     &                      0, 0, Vname(1,idItice), LBC)
            CASE ('LBC(isIsig)')
              Npts=load_lbc(Nval, Cval, line, nline, isIsig, igrid,     &
     &                      0, 0, Vname(1,idIsig), LBC)
!
! Model constants.
!
            CASE ('NSTRS')
              Npts=load_i(Nval, Rval, Ngrids, nstrs)
            CASE ('NEVP')
              Npts=load_i(Nval, Rval, Ngrids, nevp)
            CASE ('RHOICE')
              Npts=load_r(Nval, Rval, Ngrids, rhoice)
            CASE ('CDIW')
              Npts=load_r(Nval, Rval, Ngrids, cdiw)
            CASE ('CDAI')
              Npts=load_r(Nval, Rval, Ngrids, cdai)
            CASE ('RHO_AIR')
              Npts=load_r(Nval, Rval, Ngrids, rho_air)
            CASE ('RHOSNOW_DRY')
              Npts=load_r(Nval, Rval, Ngrids, rhosnow_dry)
            CASE ('RHOSNOW_WET')
              Npts=load_r(Nval, Rval, Ngrids, rhosnow_wet)
            CASE ('PSTAR')
              Npts=load_r(Nval, Rval, Ngrids, pstar)
            CASE ('PSTAR_E')
              Npts=load_r(Nval, Rval, Ngrids, pstar_e)
            CASE ('ASTREN')
              Npts=load_r(Nval, Rval, Ngrids, astren)
            CASE ('ZETAMAX')
              Npts=load_r(Nval, Rval, Ngrids, zetamax)
            CASE ('ZETAMIN')
              Npts=load_r(Nval, Rval, Ngrids, zetamin)
            CASE ('ELLIP_SQ')
              Npts=load_r(Nval, Rval, Ngrids, ellip_sq)
            CASE ('ALPHAI')
              Npts=load_r(Nval, Rval, Ngrids, alphai)
              do ng=1,Ngrids
                 alphai(ng) = alphai(ng)*deg2rad
              enddo
            CASE ('TOL')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              tol=Rvalue(1)
            CASE ('MIN_H')
              Npts=load_r(Nval, Rval, Ngrids, min_h)
            CASE ('MIN_A')
              Npts=load_r(Nval, Rval, Ngrids, min_a)
            CASE ('MAX_A')
              Npts=load_r(Nval, Rval, Ngrids, max_a)
            CASE ('STRESSANG')
              Npts=load_r(Nval, Rval, Ngrids, stressang)
              do ng=1,Ngrids
                 stressang(ng) = stressang(ng)*deg2rad
              enddo
#  ifdef ICE_LANDFAST
            CASE ('LF_K1')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              lf_k1=Rvalue(1)
            CASE ('LF_K2')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              lf_k2=Rvalue(1)
            CASE ('LF_U0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              lf_u0=Rvalue(1)
#  endif
!
! Nudging options.
!
            CASE ('LmiCLM')
              Npts=load_l(Nval, Cval, Ngrids, LmiCLM)
            CASE ('LaiCLM')
              Npts=load_l(Nval, Cval, Ngrids, LaiCLM)
            CASE ('LsiCLM')
              Npts=load_l(Nval, Cval, Ngrids, LsiCLM)
            CASE ('LnudgeMICLM')
              Npts=load_l(Nval, Cval, Ngrids, LnudgeMICLM)
            CASE ('LnudgeAICLM')
              Npts=load_l(Nval, Cval, Ngrids, LnudgeAICLM)
            CASE ('LnudgeSICLM')
              Npts=load_l(Nval, Cval, Ngrids, LnudgeSICLM)
# endif
# ifdef ICESHELF
            CASE ('GAMMA_U')
              Npts=load_r(Nval, Rval, Ngrids, gamU)
            CASE ('USTAR_MIN')
              Npts=load_r(Nval, Rval, Ngrids, ustar_min)
# endif
# if defined ICESHELF && defined ICESHELF_TRACER && defined T_PASSIVE
            CASE ('LtracerIce')
              Npts=load_l(Nval, Cval, NPT, Ngrids, LtracerI)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  LtracerIce(itrc,ng)=LtracerI(i,ng)
                END DO
              END DO
# endif
# if defined ICE_MODEL && defined ICE_FRZ_TRACER && defined T_PASSIVE
            CASE ('LtracerFrz')
              Npts=load_l(Nval, Cval, NPT, Ngrids, LtracerI)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=inert(i)
                  LtracerFrz(itrc,ng)=LtracerI(i,ng)
                END DO
              END DO
# endif
!
! Output options.
!
#include "read_icepar_check.h"
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
# ifdef ICE_MODEL
!
! Set ice time step to ocean time step
!
      DO ng = 1,Ngrids
        dtice(ng) = dt(ng)
      END DO
!
!-----------------------------------------------------------------------
!  Process input parameters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
#  ifndef ANA_MICLIMA
        IF (LmiCLM(ng)) CLM_FILE(ng)=.TRUE.
#  endif
#  ifndef ANA_AICLIMA
        IF (LaiCLM(ng)) CLM_FILE(ng)=.TRUE.
#  endif
#  ifndef ANA_SICLIMA
        IF (LsiCLM(ng)) CLM_FILE(ng)=.TRUE.
#  endif
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        DO ng=1,Ngrids
# ifdef ICE_MODEL
          IF (.not.Lice(ng)) RETURN
          WRITE (out,40) ng
          WRITE(out,*) 'Ice time step = ocean time step'
          WRITE (out,60) dtice(ng), 'DTICE',                            &
     &          'Ice model time step (s).'
          WRITE (out,50) nstrs(ng), 'NSTRS',                            &
     &          'Number of iterations for nonlinear ice dynamics.'
          WRITE (out,50) nevp(ng), 'NEVP',                              &
     &          'Number of elastic steps per plastic step in EVP.'
          WRITE (out,60) rhoice(ng), 'RHOICE',                          &
     &          'Density of sea ice (kg/m3).'
          WRITE (out,60) cdiw(ng), 'CDIW',                              &
     &          'Ice-water drag coefficient (nondimensional).'
          WRITE (out,60) cdai(ng), 'CDAI',                              &
     &          'Air-ice drag coefficient (nondimensional).'
          WRITE (out,60) rho_air(ng), 'RHO_AIR',                        &
     &          'Air density (kg/m3).'
          WRITE (out,60) rhosnow_dry(ng), 'RHOSNOW_DRY',                &
     &          'Dry snow density (kg/m3).'
          WRITE (out,60) rhosnow_wet(ng), 'RHOSNOW_WET',                &
     &          'Wet snow density (kg/m3).'
          WRITE (out,60) pstar(ng), 'PSTAR',                            &
     &          'Ice stength coefficient (N/m2).'
          WRITE (out,60) pstar_e(ng), 'PSTAR_E',                        &
     &          'Ice stength coefficient in quadratic form (N/m).'
          WRITE (out,60) alphai(ng)*rad2deg, 'ALPHAI',                  &
     &          'Mohr-Coulomb stress angle (degrees).'
          WRITE (out,60) min_h(ng), 'MIN_H',                            &
     &          'Minimum average ice thickness (m).'
          WRITE (out,60) min_a(ng), 'MIN_A',                            &
     &          'Minimum ice concentration (nondimensional).'
          WRITE (out,60) max_a(ng), 'MAX_A',                            &
     &          'Maximum ice concentration (nondimensional).'
          WRITE (out,60) stressang(ng)*rad2deg, 'STRESSANG',            &
     &          'Turning angle for ice-water drag (degrees).'
          WRITE (out,60) tol, 'TOL',                                    &
     &          'Numerical tolerance in rheology calculations .'
#  ifdef ICE_LANDFAST
          WRITE (out,60) lf_k1, 'LF_K1',                                &
     &          'Landfast ice tunable k1.'
          WRITE (out,60) lf_k2, 'LF_K2',                                &
     &          'Landfast ice tunable k2.'
          WRITE (out,60) lf_u0, 'LF_U0',                                &
     &          'Landfast ice minimal velocity.'
#  endif

          IF (LmiCLM(ng)) THEN
            WRITE (out,170) LmiCLM(ng), 'LmiCLM',                       &
     &          'Turning ON  processing of ice momentum climatology.'
          ELSE
            WRITE (out,170) LmiCLM(ng), 'LmiCLM',                       &
     &          'Turning OFF processing of ice momentum climatology.'
          END IF
          IF (LaiCLM(ng)) THEN
            WRITE (out,170) LaiCLM(ng), 'LaiCLM',                       &
     &          'Turning ON  processing of ice tracer climatology.'
          ELSE
            WRITE (out,170) LaiCLM(ng), 'LaiCLM',                       &
     &          'Turning OFF processing of ice tracer climatology.'
          END IF
          IF (LsiCLM(ng)) THEN
            WRITE (out,170) LsiCLM(ng), 'LsiCLM',                       &
     &          'Turning ON  processing of ice stress climatology.'
          ELSE
            WRITE (out,170) LsiCLM(ng), 'LsiCLM',                       &
     &          'Turning OFF processing of ice stress climatology.'
          END IF
          IF (LnudgeMICLM(ng)) THEN
            WRITE (out,170) LnudgeMICLM(ng), 'LnudgeMICLM',             &
     &          'Turning ON  nudging of ice momentum climatology.'
          ELSE
            WRITE (out,170) LnudgeMICLM(ng), 'LnudgeMICLM',             &
     &          'Turning OFF nudging of ice momentum climatology.'
          END IF
          IF (LnudgeAICLM(ng)) THEN
            WRITE (out,170) LnudgeAICLM(ng), 'LnudgeAICLM',             &
     &          'Turning ON  nudging of ice tracer climatology.'
          ELSE
            WRITE (out,170) LnudgeAICLM(ng), 'LnudgeAICLM',             &
     &          'Turning OFF nudging of ice tracer climatology.'
          END IF
          IF (LnudgeSICLM(ng)) THEN
            WRITE (out,170) LnudgeSICLM(ng), 'LnudgeSICLM',             &
     &          'Turning ON  nudging of ice stress climatology.'
          ELSE
            WRITE (out,170) LnudgeSICLM(ng), 'LnudgeSICLM',             &
     &          'Turning OFF nudging of ice stress climatology.'
          END IF
# endif
# ifdef ICESHELF
          WRITE (out,60) gamU(ng), 'GAMMA_U',                           &
     &          'Iceshelf Ustar scaling factor (nondimensional).'
          WRITE (out,60) ustar_min(ng), 'USTAR_MIN',                    &
     &          'Iceshelf minimum Ustar value.'
# endif
# if defined ICESHELF && defined ICESHELF_TRACER && defined T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (LtracerIce(itrc,ng)) THEN
              WRITE (out,185) LtracerIce(itrc,ng), 'LtracerIce', itrc,  &
     &            'Turning ON  point Sources/Sinks on iceshelf tracer ',&
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            ELSE
              WRITE (out,185) LtracerIce(itrc,ng), 'LtracerIce', itrc,  &
     &            'Turning OFF point Sources/Sinks on iceshelf tracer ',&
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END IF
          END DO
# endif
# if defined ICE_MODEL && defined ICE_FRZ_TRACER && defined T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (LtracerFrz(itrc,ng)) THEN
              WRITE (out,185) LtracerFrz(itrc,ng), 'LtracerFrz', itrc,  &
     &            'Turning ON  point Sources/Sinks on iceshelf tracer ',&
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            ELSE
              WRITE (out,185) LtracerFrz(itrc,ng), 'LtracerFrz', itrc,  &
     &            'Turning OFF point Sources/Sinks on iceshelf tracer ',&
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END IF
          END DO
# endif
!
#include "read_icepar_print.h"
!
        END DO
      END IF

  30  FORMAT (/,' READ_ICEPAR - Error while processing line: ',/,a)
  40  FORMAT (/,/,'      Ice Parameters, Grid: ',i2.2,                  &
     &        /,' ============================',/)
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (1p,e11.4,2x,a,t28,a)
  70  FORMAT (10x,l1,2x,a,t30,a)
  80  FORMAT (/,' READ_ICEPAR - Variable index not yet loaded, ', a)
 170  FORMAT (10x,l1,2x,a,t32,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
      RETURN
      END SUBROUTINE read_IcePar
#endif
