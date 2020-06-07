subroutine num_init()

  use sim_data, only: sim_Torder,   &
                      sim_order,    &
                      sim_RK,       &
                      sim_cornerBC, &
                      sim_gpWENO,   &
                      sim_gpRadii
  use grid_data, only: gr_ngc
  use num_data

  implicit none

  procedure(temporal) :: soln_RK2, soln_RK3, soln_RK4
  procedure(temporal) :: soln_sfPIF3, soln_sfPIF4
  procedure(spatial) :: soln_WENO5
  procedure(spatial) :: soln_gpWENO5, soln_gpWENO7

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

  if (sim_gpWENO) then
    call gp_WENOinit()
    num_radius = sim_gpRadii
    if (num_radius == 2) then
      num_spatial_method => soln_gpWENO5
    else if (num_radius == 3) then
      num_spatial_method => soln_gpWENO7
    else
      call abort_slug("unrecognized sim_gpRadii")
    end if
  else
    if (sim_order == 5) then
      num_radius = 2
      num_spatial_method => soln_WENO5
    else
      call abort_slug("unrecognized sim_order")
    end if
  end if

  ! check if there are sufficient guard cells.
  ! sfPIF needs more guard cells b/c num_diffs
  ! they need two more guard cells, as we store `div` in grid.
  if (.not. sim_RK) then
    if (num_radius + 1 + 2 + 2 > gr_ngc) then
      print *, "With sfPIF, at least", num_radius+1+2+2, "guard cells are required. gr_ngc =", gr_ngc
      call abort_slug("Wrong # of guard cells")
    end if
  else
    ! RK needs r+1 gc's
    if (num_radius + 1 > gr_ngc) then
      print *, "At least", num_radius+1, "guard cells are required. gr_ngc =", gr_ngc
      call abort_slug("Wrong # of guard cells")
    end if
  end if

end subroutine num_init
