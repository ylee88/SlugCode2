module num_data

#include "definition.h"

  use num_interface
  use sim_data, only: sim_order,  &
                      sim_Torder, &
                      sim_RK

  implicit none


  integer, save :: num_radius

  procedure(spatial),  pointer, save :: num_spatial_method
  procedure(temporal), pointer, save :: num_temporal_method





end module num_data
