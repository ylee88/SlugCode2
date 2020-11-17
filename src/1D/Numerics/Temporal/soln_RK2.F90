subroutine soln_RK2(dt)

#include "definition.h"

  use grid_data, only: gr_V,                &
                       gr_U,                &
                       gr_flux,             &
                       gr_dx,               &
                       gr_i0, gr_imax,      &
                       gr_ibeg, gr_iend
  use primconsflux, only: prim2cons,        &
                          prim2flux,        &
                          cons2prim
  use bc, only: bc_apply

  implicit none

  real, intent(IN) :: dt
  real :: dtx
  integer :: m, i, dir
  real, dimension(NSYS_VAR, gr_imax(XDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), NDIM) :: mFlux !(var,i,ndim)

  real, dimension(NUMB_VAR, gr_imax(XDIM)) :: prim
  real, dimension(NSYS_VAR, gr_imax(XDIM)) :: cons
  real, dimension(NSYS_VAR, gr_imax(XDIM), NDIM) :: flux

  mFlux = 0.
  dtx = dt/gr_dx

  do m = 1, 2

    ! initial data for spatial recon/intp
    prim = gr_V
    do i = gr_i0(XDIM), gr_imax(XDIM)
      call prim2cons(prim(:,i), cons(:,i))
      do dir = XDIM, NDIM
        call prim2flux(prim(:,i), flux(:,i,dir), dir)
      end do
    end do
    ! spatial recon/intp
    call soln_spatial(dt, prim, cons, flux)

    Uk = 0.
    if (m == 1) then
      do i = gr_ibeg(XDIM), gr_iend(XDIM)
        Uk(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
          dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,XDIM))
        call cons2prim(Uk(DENS_VAR:ENER_VAR,i), gr_V(DENS_VAR:GAME_VAR,i))
      end do
    end if

    call bc_apply(gr_V)

    mFlux(:,:,:) = mFlux(:,:,:) + 0.5*gr_flux(:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:) = mFlux(DENS_VAR:ENER_VAR,:,:)

  return

end subroutine soln_RK2
