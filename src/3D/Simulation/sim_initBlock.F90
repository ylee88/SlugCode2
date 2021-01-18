subroutine sim_initBlock()

#include "definition.h"

  use sim_data
  use grid_data,  only : gr_V, gr_U,          &
                         gr_xbeg, gr_xend,    &
                         gr_i0, gr_imax,      &
                         gr_xCoord, gr_dx,    &
                         gr_yCoord, gr_dy,    &
                         gr_zCoord, gr_dz
  use primconsflux, only : prim2cons
  use bc

  implicit none

  integer :: i, j, k
  integer :: i0, imax, j0, jmax, k0, kmax
  real :: x, y, z, small, r2, beta, T, dr2, dr3, E, gamm
  real :: rt3, xcntr, ycntr, zcntr

  ! for riemann problem
  real, dimension(NSYS_VAR) :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8

  real :: x0, y0, z0
  real :: xx, yy, zz

  small = 0.01
  call RANDOM_SEED()

  i0   = gr_i0(XDIM)
  imax = gr_imax(XDIM)
  j0   = gr_i0(YDIM)
  jmax = gr_imax(YDIM)
  k0   = gr_i0(ZDIM)
  kmax = gr_imax(ZDIM)

  xcntr = (gr_xend(XDIM) + gr_xbeg(XDIM))/2.
  ycntr = (gr_xend(YDIM) + gr_xbeg(YDIM))/2.
  zcntr = (gr_xend(ZDIM) + gr_xbeg(ZDIM))/2.

  !some handy constants
  gamm = sim_gamma
  !these are used for sedov
  dr2 = (3.5*MIN(gr_dx, gr_dy, gr_dz))
  dr3 = (3.5*MIN(gr_dx, gr_dy, gr_dz))**3
  E = 0.979264

  rt3 = sqrt(3.)

  do k = k0, kmax
    zz = gr_zCoord(k)
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
          z = zz - zcntr
          r2 = x**2+y**2
          T = 1. - (gamm-1.)*beta*beta*EXP(1.-r2)/(8.*gamm*PI*PI)

          gr_V(DENS_VAR,i,j,k) = (T)**(1./(gamm-1.))
          gr_V(VELX_VAR,i,j,k) = 1. - y*beta/(2.*PI)*EXP(0.5*(1.-r2))
          gr_V(VELY_VAR,i,j,k) = 1. + x*beta/(2.*PI)*EXP(0.5*(1.-r2))
          gr_V(VELZ_VAR,i,j,k) = 1.
          gr_V(PRES_VAR,i,j,k) = gr_V(DENS_VAR,i,j,k)*T

        elseif (sim_icType == '3DRP') then
          x = xx
          y = yy
          z = zz

          x0 = 0.0
          y0 = 0.0
          z0 = 0.0

          !           DENS           VELX           VELY           VELZ      PRES
          Q1 = (/   0.53125,         0.0,           0.0,           0.0,      0.4 /)
          Q2 = (/     1.0,       0.727606875,       0.0,           0.0,      1.0 /)
          Q3 = (/     0.8,           0.0,           0.0,     -0.727606875,   1.0 /)
          Q4 = (/     1.0,           0.0,       0.727606875,       0.0,      1.0 /)
          Q5 = (/     1.0,           0.0,           0.0,      0.727606875,   1.0 /)
          Q6 = (/     0.8,           0.0,      -0.727606875,       0.0,      1.0 /)
          Q7 = (/ 1.016216216,  -0.401442839,  -0.401442839,  -0.401442839,  1.4 /)
          Q8 = (/     0.8,      -0.727606875,       0.0,           0.0,      1.0 /)

          if (x >= x0 .and. y >= y0 .and. z >= z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q1(:)
          elseif (x < x0 .and. y >= y0 .and. z >= z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q2(:)
          elseif ( x < x0 .and. y < y0 .and. z >= z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q3(:)
          elseif (x >= x0 .and. y < y0 .and. z >= z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q4(:)
          elseif (x >= x0 .and. y >= y0 .and. z < z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q5(:)
          elseif (x < x0 .and. y >= y0 .and. z < z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q6(:)
          elseif ( x < x0 .and. y < y0 .and. z < z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q7(:)
          elseif (x >= x0 .and. y < y0 .and. z < z0) then
            gr_V(DENS_VAR:PRES_VAR,i,j,k) = Q8(:)
          end if
        elseif (sim_icType == 'sedov') then
          x = xx
          y = yy
          z = zz
          r2 = x**2 + y**2 + z**2
          if (sqrt(r2) < dr2) then
            gr_V(PRES_VAR,i,j,k) = (gamm-1.)*3.*E/(4.*PI*dr3)
          else
            gr_V(PRES_VAR,i,j,k) = 1.e-5
          end if

          gr_V(DENS_VAR,i,j,k) = 1.
          gr_V(VELX_VAR,i,j,k) = 0.
          gr_V(VELY_VAR,i,j,k) = 0.
          gr_V(VELZ_VAR,i,j,k) = 0.

        elseif (sim_icType == 'tgv') then
          ! see eqn. 4.13 of 10.1016/j.jcp.2020.110006
          x = xx
          y = yy
          z = zz

          gr_V(DENS_VAR,i,j,k) = 1.
          gr_V(VELX_VAR,i,j,k) = SIN(x)*COS(y)*COS(z)
          gr_V(VELY_VAR,i,j,k) = -COS(x)*SIN(y)*COS(z)
          gr_V(VELZ_VAR,i,j,k) = 0.
          gr_V(PRES_VAR,i,j,k) = 100. + 1./16.*( COS(2.*x) + COS(2.*y) )*( COS(2.*z) + 2. )

        end if    ! icType

        gr_V(GAMC_VAR,i,j,k) = sim_gamma
        gr_V(GAME_VAR,i,j,k) = sim_gamma
        gr_V(EINT_VAR,i,j,k) = gr_V(PRES_VAR,i,j,k)/(gr_V(GAME_VAR,i,j,k)-1.)/gr_V(DENS_VAR,i,j,k)

        call prim2cons(gr_V(:,i,j,k), gr_U(:,i,j,k))

      end do
    end do
  end do

  call bc_apply(gr_V)

end subroutine sim_initBlock
