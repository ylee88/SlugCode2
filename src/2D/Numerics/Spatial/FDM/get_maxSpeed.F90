subroutine get_maxSpeed(prim)

#include "definition.h"


  use eigensystem
  use primconsflux
  use grid_data, only: gr_ibeg,     &
                       gr_iend,     &
                       gr_imax,     &
                       gr_maxSpeed

  use mpi

  implicit none

  real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM)), intent(IN) :: prim

  real, dimension(NUMB_WAVE, NDIM) :: local_maxSpeed
  real, dimension(NUMB_WAVE) :: lambda

  integer i, j, var, dir, ierr

  local_maxSpeed = -1e30
  do j = gr_ibeg(YDIM), gr_iend(YDIM)
    do i = gr_ibeg(XDIM), gr_iend(XDIM)

      do dir = XDIM, NDIM
        call eigenvalues(prim(:,i,j), lambda(:), dir)
        do var = 1, NUMB_WAVE
          local_maxSpeed(var,dir) = MAX(local_maxSpeed(var,dir), ABS(lambda(var)))
        end do
      end do

    end do
  end do

  call MPI_Allreduce(local_maxSpeed, gr_maxSpeed, NUMB_WAVE*NDIM, MPI_DOUBLE_PRECISION, &
                     MPI_MAX, MPI_COMM_WORLD, ierr)

  return
end subroutine get_maxSpeed
