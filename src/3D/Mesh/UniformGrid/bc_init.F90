subroutine bc_init()

#include "definition.h"

  use block_data, only: bl_grid, bl_grid_ext,                        &
                        bl_BC, bl_corner_cornerBC, bl_corner_sideBC, &
                        bl_iProcs, bl_jProcs, bl_kProcs,             &
                        bl_i, bl_j, bl_k

  use sim_data, only: sim_bcTypex, &
                      sim_bcTypey, &
                      sim_bcTypez, &
                      sim_cornerBC

  use mpi, only: MPI_PROC_NULL

  implicit none

  ! determine extended bl_grid
  ! for normal direction.
  if (sim_bcTypex == 'periodic') then
    bl_grid_ext(0,           1:bl_jProcs, 1:bl_kProcs) = bl_grid(bl_iProcs, :, :)
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs, 1:bl_kProcs) = bl_grid(1,         :, :)
  else
    bl_grid_ext(0,           1:bl_jProcs, 1:bl_kProcs) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs, 1:bl_kProcs) = MPI_PROC_NULL
  end if

  if (sim_bcTypey == 'periodic') then
    bl_grid_ext(1:bl_iProcs, 0,           1:bl_kProcs) = bl_grid(:, bl_jProcs, :)
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1, 1:bl_kProcs) = bl_grid(:,         1, :)
  else
    bl_grid_ext(1:bl_iProcs, 0,           1:bl_kProcs) = MPI_PROC_NULL
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1, 1:bl_kProcs) = MPI_PROC_NULL
  end if

  if (sim_bcTypez == 'periodic') then
    bl_grid_ext(1:bl_iProcs, 1:bl_jProcs, 0          ) = bl_grid(:, :, bl_kProcs)
    bl_grid_ext(1:bl_iProcs, 1:bl_jProcs, bl_kProcs+1) = bl_grid(:, :,         1)
  else
    bl_grid_ext(1:bl_iProcs, 1:bl_jProcs, 0          ) = MPI_PROC_NULL
    bl_grid_ext(1:bl_iProcs, 1:bl_jProcs, bl_kProcs+1) = MPI_PROC_NULL
  end if

  ! and for corners.
  if (      sim_bcTypex == 'periodic' &
      .and. sim_bcTypey == 'periodic' &
      .and. sim_bcTypez == 'periodic') then
    ! sides
    bl_grid_ext(          0, bl_jProcs+1, 1:bl_kProcs) = bl_grid(bl_iProcs,         1, :)
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1, 1:bl_kProcs) = bl_grid(        1,         1, :)
    bl_grid_ext(          0,           0, 1:bl_kProcs) = bl_grid(bl_iProcs, bl_jProcs, :)
    bl_grid_ext(bl_iProcs+1,           0, 1:bl_kProcs) = bl_grid(        1, bl_jProcs, :)

    bl_grid_ext(          0, 1:bl_jProcs,           0) = bl_grid(bl_iProcs, :, bl_kProcs)
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs,           0) = bl_grid(        1, :, bl_kProcs)
    bl_grid_ext(          0, 1:bl_jProcs, bl_kProcs+1) = bl_grid(bl_iProcs, :,         1)
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs, bl_kProcs+1) = bl_grid(        1, :,         1)

    bl_grid_ext(1:bl_iProcs,           0,           0) = bl_grid(:, bl_jProcs, bl_kProcs)
    bl_grid_ext(1:bl_iProcs,           0, bl_kProcs+1) = bl_grid(:, bl_jProcs,         1)
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1,           0) = bl_grid(:,         1, bl_kProcs)
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1, bl_kProcs+1) = bl_grid(:,         1,         1)

    ! corners
    bl_grid_ext(          0,           0,           0) = bl_grid(bl_iProcs, bl_jProcs, bl_kProcs)
    bl_grid_ext(bl_iProcs+1,           0,           0) = bl_grid(        1, bl_jProcs, bl_kProcs)
    bl_grid_ext(          0,           0, bl_kProcs+1) = bl_grid(bl_iProcs, bl_jProcs,         1)
    bl_grid_ext(bl_iProcs+1,           0, bl_kProcs+1) = bl_grid(        1, bl_jProcs,         1)
    bl_grid_ext(          0, bl_jProcs+1,           0) = bl_grid(bl_iProcs,         1, bl_kProcs)
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1,           0) = bl_grid(        1,         1, bl_kProcs)
    bl_grid_ext(          0, bl_jProcs+1, bl_kProcs+1) = bl_grid(bl_iProcs,         1,         1)
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1, bl_kProcs+1) = bl_grid(        1,         1,         1)
  else
    ! sides
    bl_grid_ext(          0, bl_jProcs+1, 1:bl_kProcs) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1, 1:bl_kProcs) = MPI_PROC_NULL
    bl_grid_ext(          0,           0, 1:bl_kProcs) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1,           0, 1:bl_kProcs) = MPI_PROC_NULL

    bl_grid_ext(          0, 1:bl_jProcs,           0) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs,           0) = MPI_PROC_NULL
    bl_grid_ext(          0, 1:bl_jProcs, bl_kProcs+1) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, 1:bl_jProcs, bl_kProcs+1) = MPI_PROC_NULL

    bl_grid_ext(1:bl_iProcs,           0,           0) = MPI_PROC_NULL
    bl_grid_ext(1:bl_iProcs,           0, bl_kProcs+1) = MPI_PROC_NULL
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1,           0) = MPI_PROC_NULL
    bl_grid_ext(1:bl_iProcs, bl_jProcs+1, bl_kProcs+1) = MPI_PROC_NULL

    ! corners
    bl_grid_ext(0,                     0,           0) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1,           0,           0) = MPI_PROC_NULL
    bl_grid_ext(0,           bl_jProcs+1,           0) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1,           0) = MPI_PROC_NULL

    bl_grid_ext(0,                     0, bl_kProcs+1) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1,           0, bl_kProcs+1) = MPI_PROC_NULL
    bl_grid_ext(0,           bl_jProcs+1, bl_kProcs+1) = MPI_PROC_NULL
    bl_grid_ext(bl_iProcs+1, bl_jProcs+1, bl_kProcs+1) = MPI_PROC_NULL
  end if

  ! now we figure out the boundary conditions
  ! by determining neighboring block.
  ! Negative neighboring block means global boundary,
  ! while positive block means block-block boundary.

  ! normal direction,
  !                        6
  !                 4    /
  !            +----------+
  !            |          |
  !            |          |
  !         1  |    /     | 2     j    k
  !            |  _/      |       ^  7
  !            |  5       |       | /
  !            +----------+       |/
  !                 3             +-----> i
  allocate(bl_BC(2*NDIM)); bl_BC = 0

  bl_BC(1) = bl_grid_ext(bl_i-1, bl_j  , bl_k  ) ! left face
  bl_BC(2) = bl_grid_ext(bl_i+1, bl_j  , bl_k  ) ! right face
  bl_BC(3) = bl_grid_ext(bl_i,   bl_j-1, bl_k  ) ! bottom face
  bl_BC(4) = bl_grid_ext(bl_i,   bl_j+1, bl_k  ) ! top face
  bl_BC(5) = bl_grid_ext(bl_i,   bl_j,   bl_k-1) ! down face
  bl_BC(6) = bl_grid_ext(bl_i,   bl_j,   bl_k+1) ! up face

  !here we deal with the corner BCs
  ! bl_corner_sideBC :: side-corners, ngc*ngc*gr_n[x,y,z] rectangles
  !             *------12------*
  !            / |            /|
  !           1  |           2 |
  !          /   |          /  |
  !         *----+--11-----*   |
  !         |    7         |   8
  !         |    |         |   |
  !         5    |         6   |
  !         |    *----10---+---*        j    k
  !         |   /          |  /         ^  7
  !         |  3           | 4          | /
  !         | /            |/           |/
  !         *-------9------*            +-----> i

  ! bl_corner_cornerBC :: corner-corners, ngc*ngc*ngc cubes
  !             7--------------8
  !            / |            /|
  !           /  |           / |
  !          /   |          /  |
  !         5----+---------6   |
  !         |    |         |   |
  !         |    |         |   |
  !         |    |         |   |
  !         |    3---------+---4        j    k
  !         |   /          |  /         ^  7
  !         |  /           | /          | /
  !         | /            |/           |/
  !         1--------------2            +-----> i
  if (sim_cornerBC) then
    allocate(bl_corner_sideBC(12)); bl_corner_sideBC = 0
    allocate(bl_corner_cornerBC(8)); bl_corner_cornerBC = 0

    !block-block side-corners
    bl_corner_sideBC(1)   = bl_grid_ext(bl_i-1,bl_j+1,bl_k)
    bl_corner_sideBC(2)   = bl_grid_ext(bl_i+1,bl_j+1,bl_k)
    bl_corner_sideBC(3)   = bl_grid_ext(bl_i-1,bl_j-1,bl_k)
    bl_corner_sideBC(4)   = bl_grid_ext(bl_i+1,bl_j-1,bl_k)
    bl_corner_sideBC(5)   = bl_grid_ext(bl_i-1,bl_j,bl_k-1)
    bl_corner_sideBC(6)   = bl_grid_ext(bl_i+1,bl_j,bl_k-1)
    bl_corner_sideBC(7)   = bl_grid_ext(bl_i-1,bl_j,bl_k+1)
    bl_corner_sideBC(8)   = bl_grid_ext(bl_i+1,bl_j,bl_k+1)
    bl_corner_sideBC(9)   = bl_grid_ext(bl_i,bl_j-1,bl_k-1)
    bl_corner_sideBC(10)  = bl_grid_ext(bl_i,bl_j-1,bl_k+1)
    bl_corner_sideBC(11)  = bl_grid_ext(bl_i,bl_j+1,bl_k-1)
    bl_corner_sideBC(12)  = bl_grid_ext(bl_i,bl_j+1,bl_k+1)

    ! block-block corner-corners
    bl_corner_cornerBC(1) = bl_grid_ext(bl_i-1, bl_j-1, bl_k-1)
    bl_corner_cornerBC(2) = bl_grid_ext(bl_i+1, bl_j-1, bl_k-1)
    bl_corner_cornerBC(3) = bl_grid_ext(bl_i-1, bl_j-1, bl_k+1)
    bl_corner_cornerBC(4) = bl_grid_ext(bl_i+1, bl_j-1, bl_k+1)
    bl_corner_cornerBC(5) = bl_grid_ext(bl_i-1, bl_j+1, bl_k-1)
    bl_corner_cornerBC(6) = bl_grid_ext(bl_i+1, bl_j+1, bl_k-1)
    bl_corner_cornerBC(7) = bl_grid_ext(bl_i-1, bl_j+1, bl_k+1)
    bl_corner_cornerBC(8) = bl_grid_ext(bl_i+1, bl_j+1, bl_k+1)

  end if


  return

end subroutine bc_init
