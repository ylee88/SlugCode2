subroutine soln_RK3(dt)

#include "definition.h"

  use grid_data, only: gr_V,                &
                       gr_U,                &
                       gr_flux,             &
                       gr_dx, gr_dy, gr_dz, &
                       gr_i0, gr_imax,      &
                       gr_ibeg, gr_iend
  use primconsflux, only: prim2cons,        &
                          prim2flux,        &
                          cons2prim
  use bc, only: bc_apply

  implicit none

  real, intent(IN) :: dt
  real :: dtx, dty, dtz, F
  integer :: m, i, j, k, dir
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM) :: mFlux !(var,i,j,k,ndim)

  real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)) :: prim
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)) :: cons
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM) :: flux

  mFlux = 0.
  dtx = dt/gr_dx
  dty = dt/gr_dy
  dtz = dt/gr_dz

  do m = 1, 3

    ! initial data for spatial recon/intp
    prim = gr_V
    do k = gr_i0(ZDIM), gr_imax(ZDIM)
      do j = gr_i0(YDIM), gr_imax(YDIM)
        do i = gr_i0(XDIM), gr_imax(XDIM)
          call prim2cons(prim(:,i,j,k), cons(:,i,j,k))
          do dir = XDIM, NDIM
            call prim2flux(prim(:,i,j,k), flux(:,i,j,k,dir), dir)
          end do
        end do
      end do
    end do
    ! spatial recon/intp
    call soln_spatial(dt, prim, cons, flux)

    Uk = 0.
    if (m == 1) then
      do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
          do i = gr_ibeg(XDIM), gr_iend(XDIM)
            !getting the U1 state
            !U1 = U0 + k1
            Uk(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k)                             - &
              dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) - &
              dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM)) - &
              dtz*(gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))
            call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j,k), gr_V(DENS_VAR:GAME_VAR,i,j,k))
          end do
        end do
      end do

    elseif (m == 2) then
      do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
          do i = gr_ibeg(XDIM), gr_iend(XDIM)
            !getting the U2 state
            !U2 = U0 + 1/4(k1 +k2)
            !k1 = 6.*mFlux
            Uk(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k) - 0.25*( &
              dtx*( &
              6.*(mFlux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) -mFlux(DENS_VAR:ENER_VAR,i,j,k,XDIM) ) + &
              (gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) &
              ) + dty*( &
              6.*(mFlux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) -mFlux(DENS_VAR:ENER_VAR,i,j,k,YDIM) ) +&
              (gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM))  &
              ) + dtz*( &
              6.*(mFlux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) -mFlux(DENS_VAR:ENER_VAR,i,j,k,ZDIM) ) +&
              (gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))  &
              ))
            call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j,k), gr_V(DENS_VAR:GAME_VAR,i,j,k))
          end do
        end do
      end do
    end if

    call bc_apply(gr_V)

    !F is the factor that multiplies the Km flux
    if (m == 3) then
      F = 2./3.
    else
      F = 1./6.
    end if

    mFlux(:,:,:,:,:) = mFlux(:,:,:,:,:) + F*gr_flux(:,:,:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:,:,:) = mFlux(DENS_VAR:ENER_VAR,:,:,:,:)

  return

end subroutine soln_RK3
