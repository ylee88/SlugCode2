module block_data

#include "definition.h"

  implicit none

  integer, save :: bl_iProcs             ! # of procs in each direction
  integer, save :: bl_nProcs             ! total # of procs
  integer, save :: bl_ID                 ! processor ID
  integer, save :: bl_i                  ! block coordinate for a processor

  real, save :: bl_delT                  ! min dt for each block

  integer, allocatable, dimension(:), save :: bl_grid     ! virtual block topology
  integer, allocatable, dimension(:), save :: bl_grid_ext ! extended virtual block topology, useful for BC's
  integer, allocatable, dimension(:), save :: bl_BC

end module block_data
