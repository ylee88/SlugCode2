subroutine soln_update(dt)

#include "definition.h"

  use grid_data
  use primconsflux, only : cons2prim

  implicit none
  real, intent(IN) :: dt
  integer :: i
  real :: dtx

  dtx = dt/gr_dx

  !! update conservative vars
  do i = gr_ibeg(XDIM), gr_iend(XDIM)
    gr_U(DENS_VAR:ENER_VAR,i) = gr_U(DENS_VAR:ENER_VAR,i) - &
      dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,XDIM))

    call cons2prim(gr_U(DENS_VAR:ENER_VAR,i), gr_V(DENS_VAR:GAME_VAR,i))
  end do

  return
end subroutine soln_update
