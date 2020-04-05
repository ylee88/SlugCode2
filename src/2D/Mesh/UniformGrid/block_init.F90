subroutine block_init()

#include "definition.h"

  use block_data
  use sim_data
  use grid_data

  use mpi

  implicit none

  integer :: i, j, nBlockx, nBlocky

  integer :: ierr

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, bl_nProcs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, bl_ID,     ierr)

  if (bl_iProcs*bl_jProcs .ne. bl_nprocs) then
    write(*,*) 'iProcs*jProcs: ', bl_iProcs, bl_jProcs, "=/ ", bl_nProcs
    call abort_slug('iProcs*jProcs must be number of procs')
    stop
  end if

  ! make virtual topology
  allocate(bl_grid(bl_iProcs, bl_jProcs)); bl_grid = 0
  do j = 1, bl_jProcs
    do i = 1, bl_iProcs
      bl_grid(i,j) = (j-1) + (i-1)*bl_jProcs
      if (bl_grid(i,j) == bl_ID) then
        bl_i = i
        bl_j = j
        print *, bl_ID, i, j
      end if
    end do
  end do

  if (MOD(gr_nx, bl_iProcs) + MOD(gr_ny, bl_jProcs) /= 0) then
    call abort_slug("Grid points should be evenly distributed.")
  end if

  !now we assign grid information to block
  nBlockx = gr_nx / bl_iProcs
  nBlocky = gr_ny / bl_jProcs

  ! from here, gr_n[x,y] represents "local" # of grid points
  gr_glb_nx = gr_nx
  gr_glb_ny = gr_ny

  gr_nx = nBlockx
  gr_ny = nBlocky

  ! init extended virtual block topology.
  ! This will be useful to apply BC's.
  !
  !       +----------------------+
  !       |                      |
  !       |    +------------+    |
  !       |    |            |    |
  !       |    |  bl_grid   |    |
  !       |    |            |    |
  !       |    +------------+    |
  !       |      bl_grid_ext     |
  !       +----------------------+
  !
  allocate(bl_grid_ext(0:bl_iProcs+1, 0:bl_jProcs+1)) ! note that the index start from 0.
  bl_grid_ext = -1
  bl_grid_ext(1:bl_iProcs, 1:bl_jProcs) = bl_grid(:,:)

  return

end subroutine block_init
