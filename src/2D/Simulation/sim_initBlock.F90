subroutine sim_initBlock()

#include "definition.h"

  use sim_data
  use grid_data,  only : gr_V, gr_U,          &
                         gr_xbeg, gr_xend,    &
                         gr_i0, gr_imax,      &
                         gr_xCoord, gr_dx,    &
                         gr_yCoord, gr_dy
  use primconsflux, only : prim2cons
  use bc

  implicit none

  integer :: i, j
  integer :: i0, imax, j0, jmax
  real :: x, y, ranx, rany, small, r2, beta, T, dr2, E, gamm
  real :: rt3, xmin, ymin, xcntr, ycntr

  ! for riemann problem
  real, dimension(NSYS_VAR) :: Q1, Q2, Q3, Q4

  real :: x0, y0
  real :: xx, yy

  small = 0.01
  call RANDOM_SEED()

  i0   = gr_i0(XDIM)
  imax = gr_imax(XDIM)
  j0   = gr_i0(YDIM)
  jmax = gr_imax(YDIM)

  xcntr = (gr_xend(XDIM) + gr_xbeg(XDIM))/2.
  ycntr = (gr_xend(YDIM) + gr_xbeg(YDIM))/2.

  !some handy constants
  gamm = sim_gamma
  !these are used for sedov
  dr2 = (3.5*MIN(gr_dx, gr_dy))**2
  E = 1.

  rt3 = sqrt(3.)
  xmin = 1./6.
  x0 = xmin + 1./rt3

  do j = j0, jmax
    yy = gr_yCoord(j)
    do i = i0, imax
      xx = gr_xCoord(i)

      if (sim_icType == 'vortex') then
        beta = 5.
        ! DBG
        ! xcntr = 10.
        ! ycntr = 10.
        x = xx - xcntr
        y = yy - ycntr
        r2 = x**2+y**2
        T = 1. - (gamm-1.)*beta*beta*EXP(1.-r2)/(8.*gamm*PI*PI)

        gr_V(DENS_VAR,i,j) = (T)**(1./(gamm-1.))
        gr_V(VELX_VAR,i,j) = 1. - y*beta/(2.*PI)*EXP(0.5*(1.-r2))
        gr_V(VELY_VAR,i,j) = 1. + x*beta/(2.*PI)*EXP(0.5*(1.-r2))
        gr_V(PRES_VAR,i,j) = gr_V(DENS_VAR,i,j)*T

      elseif (sim_icType == 'DMR') then
        x = xx
        y = yy
        ymin = (x - xmin)*rt3
        if (y > ymin) then
          !in the shock region
          gr_V(DENS_VAR,i,j) = 8.
          gr_V(VELX_VAR,i,j) = 7.1447096
          gr_V(VELY_VAR,i,j) = -4.125
          gr_V(PRES_VAR,i,j) = 116.5
        else
          !outside shock
          gr_V(DENS_VAR,i,j) = 1.4
          gr_V(VELX_VAR,i,j) = 0.
          gr_V(VELY_VAR,i,j) = 0.
          gr_V(PRES_VAR,i,j) = 1.
        end if

      elseif (sim_icType == '2DRP-C3') then
        x = xx
        y = yy

        x0 = 0.8
        y0 = 0.8

        !      DENS     VELX     VELY     PRES
        Q1 = (/ 1.5,    0.0,     0.0,     1.5   /)
        Q2 = (/ 0.5323, 1.206,   0.0,     0.3   /)
        Q3 = (/ 0.138,  1.206,   1.206,   0.029 /)
        Q4 = (/ 0.5323, 0.0,     1.206,   0.3   /)

        if (x > x0 .and. y >= y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q1(:)
        elseif (x <= x0 .and. y >= y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q2(:)
        elseif ( x <= x0 .and. y < y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q3(:)
        elseif (x > x0 .and. y < y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q4(:)
        end if

      elseif (sim_icType == '2DRP-C5') then
        x = xx
        y = yy

        x0 = 0.5
        y0 = 0.5

        !      DENS     VELX     VELY     PRES
        Q1 = (/ 1.0,   -0.75,   -0.5,     1.0 /)
        Q2 = (/ 2.0,   -0.75,    0.5,     1.0 /)
        Q3 = (/ 1.0,    0.75,    0.5,     1.0 /)
        Q4 = (/ 3.0,    0.75,   -0.5,     1.0 /)

        if (x > x0 .and. y >= y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q1(:)
        elseif (x <= x0 .and. y >= y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q2(:)
        elseif ( x <= x0 .and. y < y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q3(:)
        elseif (x > x0 .and. y < y0) then
          gr_V(DENS_VAR:PRES_VAR,i,j) = Q4(:)
        end if

      elseif (sim_icType == 'implosion') then

        x = abs(xx) + .5*sqrt(2.)*gr_dx
        y = abs(yy) + .5*sqrt(2.)*gr_dy
        if (x + y > 0.15) then
          gr_V(DENS_VAR,i,j) = 1.
          gr_V(PRES_VAR,i,j) = 1.
        else
          gr_V(DENS_VAR,i,j) = 0.125
          gr_V(PRES_VAR,i,j) = 0.14
        end if
        gr_V(VELX_VAR:VELY_VAR,i,j) = 0.

      elseif (sim_icType == 'sedov') then
        x = xx
        y = yy
        r2 = x**2 + y**2
        if (sqrt(r2) < sqrt(dr2)) then
          gr_V(PRES_VAR,i,j) = (gamm-1.)*E/(PI*dr2)
        else
          gr_V(PRES_VAR,i,j) = 1.e-5
        end if

        gr_V(DENS_VAR,i,j) = 1.
        gr_V(VELX_VAR,i,j) = 0.
        gr_V(VELY_VAR,i,j) = 0.

      elseif (sim_icType == 'KH') then
        x = xx
        y = yy
        call RANDOM_NUMBER(ranx)
        ranx = ranx - 0.5
        call RANDOM_NUMBER(rany)
        rany = rany - 0.5

        if (abs(y) > 0.25) then
          gr_V(VELX_VAR,i,j) = -0.5
          gr_V(DENS_VAR,i,j) = 1.
        else
          gr_V(VELX_VAR,i,j) = 0.5
          gr_V(DENS_VAR,i,j) = 2.
        end if
        gr_V(PRES_VAR,i,j) = 2.5
        gr_V(VELX_VAR,i,j) = gr_V(VELX_VAR,i,j) + small*ranx
        gr_V(VELY_VAR,i,j) = small*rany

      elseif (sim_icType == 'sodx') then
        x = xx
        y = yy
        gr_V(VELY_VAR,i,j) = 0.
        if(x < .5) then
          gr_V(DENS_VAR,i,j) = 1.
          gr_V(VELX_VAR,i,j) = 0.
          gr_V(PRES_VAR,i,j) = 1.
        else
          gr_V(DENS_VAR,i,j) = .125
          gr_V(VELX_VAR,i,j) = 0.
          gr_V(PRES_VAR,i,j) = .1
        end if

      end if    ! icType

      gr_V(GAMC_VAR,i,j) = sim_gamma
      gr_V(GAME_VAR,i,j) = sim_gamma
      gr_V(EINT_VAR,i,j) = gr_V(PRES_VAR,i,j)/(gr_V(GAME_VAR,i,j)-1.)/gr_V(DENS_VAR,i,j)

      call prim2cons(gr_V(:,i,j), gr_U(:,i,j))

    end do
  end do

  call bc_apply(gr_V)

end subroutine sim_initBlock
