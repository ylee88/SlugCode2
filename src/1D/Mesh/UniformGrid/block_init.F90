subroutine block_init()

#include "definition.h"

  use block_data
  use sim_data
  use grid_data

  use mpi

  implicit none

  integer :: i, nBlockx

  integer :: ierr

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, bl_nProcs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, bl_ID,     ierr)

  if (bl_iProcs .ne. bl_nprocs) then
    write(*,*) 'iProcs: ', bl_iProcs, "=/ ", bl_nProcs
    call abort_slug('iProcs must be number of procs')
    stop
  end if

  ! make virtual topology
  allocate(bl_grid(bl_iProcs)); bl_grid = 0
  do i = 1, bl_iProcs
    bl_grid(i) = (i-1)
    if (bl_grid(i) == bl_ID) then
      bl_i = i
      print *, bl_grid(i), i
    end if
  end do

  if (MOD(gr_nx, bl_iProcs) /= 0) then
    call abort_slug("Grid points should be evenly distributed.")
  end if

  !now we assign grid information to block
  nBlockx = gr_nx / bl_iProcs

  ! from here, gr_nx represents "local" # of grid points
  gr_glb_nx = gr_nx

  gr_nx = nBlockx

  ! check if there are too much of guard cells.
  if (gr_nx < gr_ngc) call abort_slug("Wrong # of guard cells: gr_nx < gr_ngc")

  ! init extended virtual block topology.
  ! This will be useful to apply BC's.
  !
  !            +------------+
  !               bl_grid
  !       +----------------------+
  !              bl_grid_ext
  !
  allocate(bl_grid_ext(0:bl_iProcs+1)) ! note that the index start from 0.
  bl_grid_ext = -1
  bl_grid_ext(1:bl_iProcs) = bl_grid(:)

  return

end subroutine block_init
