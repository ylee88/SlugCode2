subroutine sim_initBlock()

#include "definition.h"

  use sim_data
  use grid_data,  only : gr_V, gr_U,          &
                         gr_xbeg, gr_xend,    &
                         gr_i0, gr_imax,      &
                         gr_xCoord, gr_dx
  use primconsflux, only : prim2cons
  use bc

  implicit none

  integer :: i
  integer :: i0, imax
  real :: x, beta, gamm
  real :: rt3, xcntr, small

  ! for shock problem
  real, dimension(NSYS_VAR) :: Q1, Q2, Q3

  real :: x0
  real :: xx

  small = 0.01
  call RANDOM_SEED()

  i0   = gr_i0(XDIM)
  imax = gr_imax(XDIM)

  xcntr = (gr_xend(XDIM) + gr_xbeg(XDIM))/2.

  !some handy constants
  gamm = sim_gamma

  rt3 = sqrt(3.)

  do i = i0, imax
    xx = gr_xCoord(i)

    if (sim_icType == 'gauss') then
      beta = 100.
      x = xx - xcntr

      gr_V(DENS_VAR,i) = 1. + EXP(-beta*x**2)
      gr_V(VELX_VAR,i) = 1.
      gr_V(PRES_VAR,i) = 1./gamm

    elseif (sim_icType == 'sine') then
      x = xx - xcntr

      gr_V(DENS_VAR,i) = 1.5 - 0.5*sin(2.*PI*x)
      gr_V(VELX_VAR,i) = 1.
      gr_V(PRES_VAR,i) = 1./gamm

    elseif (sim_icType == 'sod') then
      x = xx

      !        DENS      VELX     PRES
      Q1 = (/   1.0,     0.0,     1.0 /)    ! left
      Q2 = (/ 0.125,     0.0,     0.1 /)    ! right

      if (x <= xcntr) then
        gr_V(DENS_VAR:PRES_VAR,i) = Q1(:)
      else
        gr_V(DENS_VAR:PRES_VAR,i) = Q2(:)
      end if

    elseif (sim_icType == 'blast') then
      x = xx

      !        DENS      VELX     PRES
      Q1 = (/   1.0,     0.0,   1000.0 /)    ! left
      Q2 = (/   1.0,     0.0,     0.01 /)    ! middle
      Q3 = (/   1.0,     0.0,    100.0 /)    ! right

      if (x <= 0.1) then
        gr_V(DENS_VAR:PRES_VAR,i) = Q1(:)
      elseif (x <= 0.9 .and. x > 0.1) then
        gr_V(DENS_VAR:PRES_VAR,i) = Q2(:)
      else
        gr_V(DENS_VAR:PRES_VAR,i) = Q3(:)
      end if

    elseif (sim_icType == 'shu') then
      x = xx - xcntr

      !             DENS            VELX        PRES
      Q1 = (/     3.857143,       2.629369,  10.33333 /)    ! left
      Q2 = (/ 1 + .2*SIN(5.*x),      0.0,       1.0   /)    ! right

      if (x < -4.0) then
        gr_V(DENS_VAR:PRES_VAR,i) = Q1(:)
      else
        gr_V(DENS_VAR:PRES_VAR,i) = Q2(:)
      end if

    elseif (sim_icType == 'tt') then
      x = xx

      !                  DENS                  VELX        PRES
      Q1 = (/          1.515695,             0.523346,    1.805 /)    ! left
      Q2 = (/ 1. + 0.1*SIN(20.*PI*(x-5.)),      0.0,       1.0   /)    ! right

      if (x < 0.5) then
        gr_V(DENS_VAR:PRES_VAR,i) = Q1(:)
      else
        gr_V(DENS_VAR:PRES_VAR,i) = Q2(:)
      end if

    end if    ! icType

    gr_V(GAMC_VAR,i) = sim_gamma
    gr_V(GAME_VAR,i) = sim_gamma
    gr_V(EINT_VAR,i) = gr_V(PRES_VAR,i)/(gr_V(GAME_VAR,i)-1.)/gr_V(DENS_VAR,i)

    call prim2cons(gr_V(:,i), gr_U(:,i))

  end do

  call bc_apply(gr_V)

end subroutine sim_initBlock
