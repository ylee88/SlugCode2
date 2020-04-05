subroutine soln_RK2(dt)

#include "definition.h"

  use grid_data, only: gr_V, &
                       gr_flux, &
                       gr_imax
  use primconsflux
  use bc, only: bc_apply

  implicit none
  real, intent(IN) :: dt
  real :: dtx,dty
  integer :: m, i, j
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: Uk
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: Flux !(var,i,j,ndim)

  ! Vk = gr_V
  Flux = 0.

  do m = 1, 2

    call soln_spatial(dt)

    dtx = dt/gr_dx
    dty = dt/gr_dy
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
    ! gr_V = Vk
    Flux(:,:,:,:) = Flux(:,:,:,:) + 0.5*gr_flux(:,:,:,:)

  end do

  gr_flux(DENS_VAR:ENER_VAR,:,:,:) = Flux(DENS_VAR:ENER_VAR,:,:,:)

  return

end subroutine soln_RK2
