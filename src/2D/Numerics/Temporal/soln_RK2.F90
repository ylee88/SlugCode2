subroutine soln_RK2(dt)

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
  integer :: m, i, j, dir
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: mFlux !(var,i,j,ndim)

  real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: prim
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: cons
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: flux

  mFlux = 0.
  dtx = dt/gr_dx
  dty = dt/gr_dy

  do m = 1, 2

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

    Uk = 0.
    if (m == 1) then
      do j = gr_ibeg(YDIM), gr_iend(YDIM)
        do i = gr_ibeg(XDIM), gr_iend(XDIM)
          Uk(DENS_VAR:ENER_VAR,i,j) = gr_U(DENS_VAR:ENER_VAR,i,j)                             - &
            dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,XDIM)) - &
            dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,YDIM))
          call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j), gr_V(DENS_VAR:GAME_VAR,i,j))
        end do
      end do
    end if

    call bc_apply(gr_V)

    mFlux(:,:,:,:) = mFlux(:,:,:,:) + 0.5*gr_flux(:,:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:,:) = mFlux(DENS_VAR:ENER_VAR,:,:,:)

  return

end subroutine soln_RK2
