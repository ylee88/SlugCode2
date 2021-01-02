subroutine soln_RK44(dt)
  ! traditional, nonSSP RK4 method


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
  real :: dtx, F, A
  integer :: m, i, dir
  real, dimension(NSYS_VAR, gr_imax(XDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), NDIM) :: mFlux !(var,i,j,k,ndim)

  real, dimension(NUMB_VAR, gr_imax(XDIM)) :: prim
  real, dimension(NSYS_VAR, gr_imax(XDIM)) :: cons
  real, dimension(NSYS_VAR, gr_imax(XDIM), NDIM) :: flux

  mFlux = 0.

  do m = 1, 4

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

    if (m == 1 .OR. m == 2) then
      A = 0.5
    else 
      A = 1.
    end if

    if (m == 1 .OR. m == 4) then
      F = 1./6.
    else
      F = 1./3.
    end if

    dtx = A*dt/gr_dx

    Uk = 0.

    if (m .NE. 4) then
      !update cons variables to mth step only if not 4th step
      do i = gr_ibeg(XDIM), gr_iend(XDIM)
        Uk(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
          dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,XDIM))
          call cons2prim(Uk(DENS_VAR:ENER_VAR,i), gr_V(DENS_VAR:GAME_VAR,i))
      end do
    end if

    call bc_apply(gr_V)

    mFlux(:,:,:) = mFlux(:,:,:) + F*gr_flux(:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:) = mFlux(DENS_VAR:ENER_VAR,:,:)

  return

end subroutine soln_RK44
