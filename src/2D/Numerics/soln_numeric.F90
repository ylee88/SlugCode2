! Let's assume it is FDM
subroutine soln_numeric(dt)

#include "definition.h"

  use num_data, only: num_temporal_method

  implicit none

  real, intent(IN) :: dt


  call num_temporal_method(dt)

  call apply_MPPLimiter(dt)


end subroutine soln_numeric
