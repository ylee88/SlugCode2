subroutine num_init()

  use sim_data
  use num_data

  implicit none

  procedure(temporal) :: soln_RK2, soln_RK3, soln_RK4
  procedure(spatial) :: soln_WENO5

  if (sim_RK) then
    if (sim_Torder == 2) then
      num_temporal_method => soln_RK2
    else if (sim_Torder == 3) then
      num_temporal_method => soln_RK3
    else if (sim_Torder == 4) then
      num_temporal_method => soln_RK4
    else
      call abort_slug("Unrecognized sim_Torder")
    end if
  else
    call abort_slug("non-RK method is not supported yet")
  end if

  if (sim_order == 5) then
    num_radius = 2
    num_spatial_method => soln_WENO5
  else
    call abort_slug("unrecognized sim_order")
  end if



end subroutine num_init
