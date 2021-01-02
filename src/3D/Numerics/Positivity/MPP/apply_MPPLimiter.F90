subroutine apply_MPPLimiter(dt)

#include "definition.h"

  !!!!!!!!!!!!!!!!!!!!!!!
  ! NOT YET IMPLEMENTED !
  !!!!!!!!!!!!!!!!!!!!!!!

end subroutine apply_MPPLimiter


function get_pres(U) result(pres)

  use sim_data,  only: sim_gamma

  implicit none

  real, dimension(NSYS_VAR), intent(IN) :: U
  real :: ekin, eint, pres

  ekin = 0.5*dot_product(U(MOMX_VAR:MOMZ_VAR), U(MOMX_VAR:MOMZ_VAR))/U(DENS_VAR)
  eint = U(ENER_VAR) - ekin    ! eint = rho*e

  pres = (sim_gamma-1.)*eint

end function get_pres
