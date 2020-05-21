subroutine cfl(dt)

#include "definition.h"

  use sim_data,   only: sim_cfl
  use block_data, only: bl_delT
  use grid_data,  only: gr_V,             &
                        gr_ibeg, gr_iend, &
                        gr_dx, gr_dy, gr_dz

  use mpi

  implicit none

  real, intent(OUT) :: dt
  integer :: i, j, k
  real :: cs, u, v, w, delT

  integer :: ierr

  bl_delT = 1e30
  ! update conservative vars
  do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
    do j = gr_ibeg(YDIM), gr_iend(YDIM)
      do i = gr_ibeg(XDIM), gr_iend(XDIM)

        cs = sqrt(gr_V(GAMC_VAR,i,j,k)*gr_V(PRES_VAR,i,j,k)/gr_V(DENS_VAR,i,j,k))
        u = gr_V(VELX_VAR,i,j,k); v = gr_V(VELY_VAR,i,j,k); w = gr_V(VELZ_VAR,i,j,k)

        u = abs(u) + cs
        v = abs(v) + cs
        w = abs(w) + cs
        bl_delT = min(gr_dx/u, gr_dy/v, gr_dz/w, bl_delT)

      end do
    end do
  end do

  ! reduce bl_delT
  call MPI_Allreduce(bl_delT, delT, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

  dt = sim_cfl*delT

  return
end subroutine cfl
