subroutine num_init()

  use sim_data, only: sim_Torder, &
                      sim_order,  &
                      sim_RK,     &
                      sim_cornerBC
  use num_data

  implicit none

  procedure(temporal) :: soln_RK2, soln_RK3, soln_RK4
  procedure(temporal) :: soln_sfPIF3, soln_sfPIF4
  procedure(spatial) :: soln_WENO5

  if (sim_RK) then
    ! RK doesn't need corner exchanges
    sim_cornerBC = .false.
    if (sim_Torder == 2) then
      num_temporal_method => soln_RK2
    else if (sim_Torder == 3) then
      num_temporal_method => soln_RK3
    else if (sim_Torder == 4) then
      num_temporal_method => soln_RK4
    else
      call abort_slug("unrecognized sim_Torder")
    end if
  else
    ! sfPIF needs corner exchanges
    if (.not. sim_cornerBC) call abort_slug("sfPIF needs sim_cornerBC")
    if (sim_Torder == 3) then
      num_temporal_method => soln_sfPIF3
    else if (sim_Torder == 4) then
      num_temporal_method => soln_sfPIF4
    else
      call abort_slug("unrecognized sim_RK & sim_Torder")
    end if
  end if

  if (sim_order == 5) then
    num_radius = 2
    num_spatial_method => soln_WENO5
  else
    call abort_slug("unrecognized sim_order")
  end if



end subroutine num_init
