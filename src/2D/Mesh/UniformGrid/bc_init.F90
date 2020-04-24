subroutine bc_init()

#include "definition.h"

  use block_data, only: bl_grid, bl_grid_ext, &
                        bl_BC, bl_cornerBC,   &
                        bl_iProcs, bl_jProcs, &
                        bl_i, bl_j

  use sim_data, only: sim_bcTypex, &
                      sim_bcTypey, &
                      sim_cornerBC

  use mpi, only: MPI_PROC_NULL

  implicit none

  ! determine extended bl_grid
  ! for normal direction.
  if (sim_bcTypex == 'periodic') then
    bl_grid_ext(0,           1:bl_jProcs) = bl_grid(bl_iProcs, :)
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs) = bl_grid(1,         :)
  else
    bl_grid_ext(0,           1:bl_jProcs) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs) = MPI_PROC_NULL
  end if

  if (sim_bcTypey == 'periodic') then
    bl_grid_ext(1:bl_iProcs, 0          ) = bl_grid(:, bl_jProcs)
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1) = bl_grid(:,         1)
  else
    bl_grid_ext(1:bl_iProcs, 0          ) = MPI_PROC_NULL
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1) = MPI_PROC_NULL
  end if
  ! and for corners.
  if (      sim_bcTypex == 'periodic' &
      .and. sim_bcTypey == 'periodic' ) then
    bl_grid_ext(0,                     0) = bl_grid(bl_iProcs, bl_jProcs)
    bl_grid_ext(bl_iProcs+1,           0) = bl_grid(1,         bl_jProcs)
    bl_grid_ext(0,           bl_jProcs+1) = bl_grid(bl_iProcs,         1)
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1) = bl_grid(1,                 1)
  else
    bl_grid_ext(0,                     0) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1,           0) = MPI_PROC_NULL
    bl_grid_ext(0,           bl_jProcs+1) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1) = MPI_PROC_NULL
  end if

  ! now we figure out the boundary conditions
  ! by determining neighboring block.
  ! Negative neighboring block means global boundary,
  ! while positive block means block-block boundary.

  ! normal direction,
  !                 4
  !            +----------+
  !            |          |
  !            |          |
  !         1  |          | 2     j
  !            |          |       ^
  !            |          |       |
  !            +----------+       |
  !                 3             +-----> i
  allocate(bl_BC(2*NDIM)); bl_BC = 0

  bl_BC(1) = bl_grid_ext(bl_i-1, bl_j  ) ! left face
  bl_BC(2) = bl_grid_ext(bl_i+1, bl_j  ) ! right face
  bl_BC(3) = bl_grid_ext(bl_i,   bl_j-1) ! bottom face
  bl_BC(4) = bl_grid_ext(bl_i,   bl_j+1) ! top face
  ! corner direction,
  !
  !          1  ---------- 2
  !            |          |
  !            |          |
  !            |          |       j
  !            |          |       ^
  !            |          |       |
  !          4  ----------  3     |
  !                               +-----> i
  if (sim_cornerBC) then
    allocate(bl_cornerBC(4)); bl_cornerBC = 0

    bl_cornerBC(1) = bl_grid_ext(bl_i-1,bl_j+1)
    bl_cornerBC(2) = bl_grid_ext(bl_i+1,bl_j+1)
    bl_cornerBC(3) = bl_grid_ext(bl_i+1,bl_j-1)
    bl_cornerBC(4) = bl_grid_ext(bl_i-1,bl_j-1)
  end if

  return

end subroutine bc_init
