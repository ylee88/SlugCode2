subroutine num_init()

  use sim_data, only: sim_spatialMethod,  &
                      sim_temporalMethod, &
                      sim_Torder,         &
                      sim_order,          &
                      sim_RK,             &
                      sim_cornerBC,       &
                      sim_gpWENO,         &
                      sim_gpRadii
  use grid_data, only: gr_ngc
  use num_data

  implicit none

  procedure(temporal) :: soln_RK2, soln_RK3, soln_RK4
  procedure(temporal) :: soln_sfPIF3, soln_sfPIF4
  procedure(spatial) :: soln_WENO5, soln_nestedWENO5
  procedure(spatial) :: soln_gpWENO5, soln_gpWENO7


  ! choose temporal method
  select case(sim_temporalMethod)
  ! RK
  case('RK2')
    num_temporal_method => soln_RK2
    sim_Torder = 2
    sim_cornerBC = .false.
    sim_RK = .true.
  case('RK3')
    num_temporal_method => soln_RK3
    sim_Torder = 3
    sim_cornerBC = .false.
    sim_RK = .true.
  case('RK4')
    num_temporal_method => soln_RK4
    sim_Torder = 4
    sim_cornerBC = .false.
    sim_RK = .true.
  ! SF-PIF
  case('SF3')
    num_temporal_method => soln_sfPIF3
    sim_Torder = 3
    if (.not. sim_cornerBC) call abort_slug("sfPIF needs sim_cornerBC")
    sim_RK = .false.
  case('SF4')
    num_temporal_method => soln_sfPIF4
    sim_Torder = 4
    if (.not. sim_cornerBC) call abort_slug("sfPIF needs sim_cornerBC")
    sim_RK = .false.
  case DEFAULT
    call abort_slug("unrecognized sim_temporalMethod: "//sim_temporalMethod)
  end select



  ! choose spatial method
  select case(sim_spatialMethod)
  case('WENO5')
    num_spatial_method => soln_WENO5
    sim_order = 5
    num_radius = 2
  case('nestedWENO5')
    num_spatial_method => soln_nestedWENO5
    sim_order = 5
    num_radius = 2
  case('gpWENO5')
    num_spatial_method => soln_gpWENO5
    call gp_WENOinit()
    sim_order = 5
    sim_gpWENO = .true.
    num_radius = sim_gpRadii
  case('gpWENO7')
    num_spatial_method => soln_gpWENO7
    call gp_WENOinit()
    sim_order = 7
    sim_gpWENO = .true.
    num_radius = sim_gpRadii
  case DEFAULT
    call abort_slug("unrecognized sim_spatialMethod: "//sim_spatialMethod)
  end select


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
