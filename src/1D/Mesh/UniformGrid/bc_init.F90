subroutine bc_init()

#include "definition.h"

  use block_data, only: bl_grid, bl_grid_ext,  &
                        bl_BC,                 &
                        bl_iProcs,             &
                        bl_i

  use sim_data, only: sim_bcTypex

  use mpi, only: MPI_PROC_NULL

  implicit none

  ! determine extended bl_grid
  ! for normal direction.
  if (sim_bcTypex == 'periodic') then
    bl_grid_ext(0          ) = bl_grid(bl_iProcs)
    bl_grid_ext(bl_iProcs+1) = bl_grid(1        )
  else
    bl_grid_ext(0          ) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1) = MPI_PROC_NULL
  end if

  ! now we figure out the boundary conditions
  ! by determining neighboring block.
  ! Negative neighboring block means global boundary,
  ! while positive block means block-block boundary.

  allocate(bl_BC(2*NDIM)); bl_BC = 0

  bl_BC(1) = bl_grid_ext(bl_i-1) ! left face
  bl_BC(2) = bl_grid_ext(bl_i+1) ! right face


  return

end subroutine bc_init
