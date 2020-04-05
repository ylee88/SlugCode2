subroutine soln_update(dt)

#include "definition.h"

  use grid_data
  use primconsflux, only : cons2prim

  implicit none
  real, intent(IN) :: dt
  integer :: i, j
  real :: dtx, dty

  dtx = dt/gr_dx
  dty = dt/gr_dy

  !! update conservative vars
  do j = gr_ibeg(YDIM), gr_iend(YDIM)
    do i = gr_ibeg(XDIM), gr_iend(XDIM)
      gr_U(DENS_VAR:ENER_VAR,i,j) = gr_U(DENS_VAR:ENER_VAR,i,j) - &
        dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,XDIM)) - &
        dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,YDIM))

      call cons2prim(gr_U(DENS_VAR:ENER_VAR,i,j), gr_V(DENS_VAR:GAME_VAR,i,j))
    end do
  end do


  return
end subroutine soln_update
