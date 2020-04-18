subroutine cfl(dt)

#include "definition.h"

  use sim_data,   only: sim_cfl
  use block_data, only: bl_delT
  use grid_data,  only: gr_V,             &
                        gr_ibeg, gr_iend, &
                        gr_dx, gr_dy

  use mpi

  implicit none

  real, intent(OUT) :: dt
  integer :: i, j
  real :: cs, u, v, delT

  integer :: ierr

  bl_delT = 1e30
  ! update conservative vars
  do j = gr_ibeg(YDIM), gr_iend(YDIM)
    do i = gr_ibeg(XDIM), gr_iend(XDIM)

      cs = sqrt(gr_V(GAMC_VAR,i,j)*gr_V(PRES_VAR,i,j)/gr_V(DENS_VAR,i,j))
      u = gr_V(VELX_VAR,i,j); v = gr_V(VELY_VAR,i,j)

      u = abs(u) + cs
      v = abs(v) + cs
      bl_delT = min(gr_dx/u, gr_dy/v, bl_delT)

    end do
  end do

  ! reduce bl_delT
  call MPI_Allreduce(bl_delT, delT, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

  dt = sim_cfl*delT

  return
end subroutine cfl
