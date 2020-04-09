subroutine soln_RK3(dt)

#include "definition.h"

  use grid_data, only: gr_V, &
                       gr_flux, &
                       gr_imax
  use primconsflux
  use bc, only: bc_apply

  implicit none
  real, intent(IN) :: dt
  real :: dtx,dty,F
  integer :: m, i, j
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: Flux !(var,i,j,ndim)

  Flux = 0.
  dtx = dt/gr_dx
  dty = dt/gr_dy

  do m = 1, 3

    call soln_spatial(dt)

    Uk = 0.
    if (m == 1) then
      do j = gr_ibeg(YDIM), gr_iend(YDIM)
        do i = gr_ibeg(XDIM), gr_iend(XDIM)
          !getting the U1 state
          !U1 = U0 + k1
          Uk(DENS_VAR:ENER_VAR,i,j) = gr_U(DENS_VAR:ENER_VAR,i,j)                             - &
            dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,XDIM)) - &
            dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,YDIM))
          call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j), gr_V(DENS_VAR:GAME_VAR,i,j))
        end do
      end do

    elseif (m == 2) then
      do j = gr_ibeg(YDIM), gr_iend(YDIM)
        do i = gr_ibeg(XDIM), gr_iend(XDIM)
          !getting the U2 state
          !U2 = U0 + 1/4(k1 +k2)
          !k1 = 6.*Flux
          Uk(DENS_VAR:ENER_VAR,i,j) = gr_U(DENS_VAR:ENER_VAR,i,j) - 0.25*( &
            dtx*( &
            6.*(Flux(DENS_VAR:ENER_VAR,i+1,j,XDIM) -Flux(DENS_VAR:ENER_VAR,i,j,XDIM) ) + &
            (gr_flux(DENS_VAR:ENER_VAR,i+1,j,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,XDIM)) &
            ) + dty*( &
            6.*(Flux(DENS_VAR:ENER_VAR,i,j+1,YDIM) -Flux(DENS_VAR:ENER_VAR,i,j,YDIM) ) +&
            (gr_flux(DENS_VAR:ENER_VAR,i,j+1,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,YDIM))  &
            ))
          call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j), gr_V(DENS_VAR:GAME_VAR,i,j))
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

    Flux(:,:,:,:) = Flux(:,:,:,:) + F*gr_flux(:,:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:,:) = Flux(DENS_VAR:ENER_VAR,:,:,:)

  return

end subroutine soln_RK3
