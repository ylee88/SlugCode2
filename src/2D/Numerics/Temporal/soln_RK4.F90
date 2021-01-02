subroutine soln_RK4(dt)
  ! 5 stages, 4th-order RK method.
  ! ref: 10.1137/S0036142901389025

#include "definition.h"

  use grid_data, only: gr_V,             &
                       gr_U,             &
                       gr_flux,          &
                       gr_dx, gr_dy,     &
                       gr_i0, gr_imax,   &
                       gr_ibeg, gr_iend
  use primconsflux, only: prim2cons,     &
                          prim2flux,     &
                          cons2prim
  use bc, only: bc_apply

  implicit none

  real, intent(IN) :: dt
  real :: dtx, dty
  integer :: m, n, i, j, dir
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: mFlux !(var,i,j,ndim)

  real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: prim
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: cons
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: flux
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM, 4) :: nFlux

  real, dimension(5) :: F
  real, dimension(4,4) :: A

  A(:, :) = 0.
  A(1  , 1) = 0.39175222657189
  A(1:2, 2) = (/ 0.217669096261169, 0.368410593050371  /)
  A(1:3, 3) = (/ 0.0826920866578106, 0.139958502191895, 0.251891774271694 /)
  A(1:4, 4) = (/ 0.0679662836371148, 0.115034698504632, 0.207034898597386, 0.544974750228521 /)

  F = (/ 0.146811876084786, 0.248482909444976, 0.104258830331981, 0.27443890090135, 0.226007483236906 /)

  mFlux = 0.
  nFlux = 0.

  do m = 1, 5

    ! initial data for spatial recon/intp
    prim = gr_V
    do j = gr_i0(YDIM), gr_imax(YDIM)
      do i = gr_i0(XDIM), gr_imax(XDIM)
        call prim2cons(prim(:,i,j), cons(:,i,j))
        do dir = XDIM, NDIM
          call prim2flux(prim(:,i,j), flux(:,i,j,dir), dir)
        end do
      end do
    end do
    ! spatial recon/intp
    call soln_spatial(dt, prim, cons, flux)

    dtx = dt/gr_dx
    dty = dt/gr_dy

    if (m /= 5) then

      nFlux(:,:,:,:,m) = gr_flux(:,:,:,:)      ! save interstage fluxes

      Uk = gr_U
      do n = 1, m
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
          do i = gr_ibeg(XDIM), gr_iend(XDIM)
            ! uk = uk-1 - dt*A0*L(u0) - ... - dt*Ak-1*L(uk-1)
            Uk(DENS_VAR:ENER_VAR,i,j) = Uk(DENS_VAR:ENER_VAR,i,j) - &
              A(n,m)*dtx*(nFlux(DENS_VAR:ENER_VAR,i+1,j,XDIM,n) - nFlux(DENS_VAR:ENER_VAR,i,j,XDIM,n)) - &
              A(n,m)*dty*(nFlux(DENS_VAR:ENER_VAR,i,j+1,YDIM,n) - nFlux(DENS_VAR:ENER_VAR,i,j,YDIM,n))
            call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j), gr_V(DENS_VAR:GAME_VAR,i,j))
          end do
        end do
      end do

    end if

    call bc_apply(gr_V)

    mFlux(:,:,:,:) = mFlux(:,:,:,:) + F(m)*gr_flux(:,:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:,:) = mFlux(DENS_VAR:ENER_VAR,:,:,:)

  return

end subroutine soln_RK4
