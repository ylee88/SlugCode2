subroutine soln_RK4(dt)

#include "definition.h"

  use grid_data, only: gr_V, &
                       gr_flux, &
                       gr_imax
  use primconsflux
  use bc, only: bc_apply

  implicit none
  real, intent(IN) :: dt
  real :: dtx, dty, F, A
  integer :: m, i, j
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: Flux !(var,i,j,ndim)

  Flux = 0.

  do m = 1, 4

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
    dty = A*dt/gr_dy

    call soln_spatial(dt)

    Uk = 0.

    if (m .NE. 4) then
      !update cons variables to mth step only if not 4th step
      do j = gr_ibeg(YDIM), gr_iend(YDIM)
        do i = gr_ibeg(XDIM), gr_iend(XDIM)
          Uk(DENS_VAR:ENER_VAR,i,j) = gr_U(DENS_VAR:ENER_VAR,i,j) - &
            dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,XDIM)) - &
            dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,YDIM))
          call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j), gr_V(DENS_VAR:GAME_VAR,i,j))
        end do
      end do
    end if

    call bc_apply(gr_V)

    Flux(:,:,:,:) = Flux(:,:,:,:) + F*gr_flux(:,:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:,:) = Flux(DENS_VAR:ENER_VAR,:,:,:)

  return

end subroutine soln_RK4
