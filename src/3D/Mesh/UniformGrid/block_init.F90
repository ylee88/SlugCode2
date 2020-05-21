subroutine block_init()

#include "definition.h"

  use block_data
  use sim_data
  use grid_data

  use mpi

  implicit none

  integer :: i, j, k, nBlockx, nBlocky, nBlockz

  integer :: ierr

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, bl_nProcs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, bl_ID,     ierr)

  if (bl_iProcs*bl_jProcs*bl_kProcs .ne. bl_nprocs) then
    write(*,*) 'iProcs*jProcs*kProcs: ', bl_iProcs, bl_jProcs, bl_kProcs, "=/ ", bl_nProcs
    call abort_slug('iProcs*jProcs*kProcs must be number of procs')
    stop
  end if

  ! make virtual topology
  allocate(bl_grid(bl_iProcs, bl_jProcs, bl_kProcs)); bl_grid = 0
  do i = 1, bl_iProcs
    do j = 1, bl_jProcs
      do k = 1, bl_kProcs
        bl_grid(i,j,k) = (k-1) + (j-1)*bl_iProcs + (i-1)*bl_jProcs*bl_kProcs
        if (bl_grid(i,j,k) == bl_ID) then
          bl_i = i
          bl_j = j
          bl_k = k
          print *, bl_grid(i,j,k), i, j, k
        end if
      end do
    end do
  end do

  if (MOD(gr_nx, bl_iProcs) + MOD(gr_ny, bl_jProcs) + MOD(gr_nz, bl_kProcs) /= 0) then
    call abort_slug("Grid points should be evenly distributed.")
  end if

  !now we assign grid information to block
  nBlockx = gr_nx / bl_iProcs
  nBlocky = gr_ny / bl_jProcs
  nBlockz = gr_nz / bl_kProcs

  ! from here, gr_n[x,y] represents "local" # of grid points
  gr_glb_nx = gr_nx
  gr_glb_ny = gr_ny
  gr_glb_nz = gr_nz

  gr_nx = nBlockx
  gr_ny = nBlocky
  gr_nz = nBlockz

  ! check if there are too much of guard cells.
  if (gr_nx < gr_ngc) call abort_slug("Wrong # of guard cells: gr_nx < gr_ngc")
  if (gr_ny < gr_ngc) call abort_slug("Wrong # of guard cells: gr_ny < gr_ngc")
  if (gr_nz < gr_ngc) call abort_slug("Wrong # of guard cells: gr_nz < gr_ngc")

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
  allocate(bl_grid_ext(0:bl_iProcs+1, 0:bl_jProcs+1, 0:bl_kProcs+1)) ! note that the index start from 0.
  bl_grid_ext = -1
  bl_grid_ext(1:bl_iProcs, 1:bl_jProcs, 1:bl_kProcs) = bl_grid(:,:,:)

  return

end subroutine block_init
