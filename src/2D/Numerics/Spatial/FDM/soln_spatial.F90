subroutine soln_spatial(dt)

#include "definition.h"

  use num_data,  only: num_radius, num_spatial_method
  use grid_data, only: gr_V, &
                       gr_fL, gr_fR, &
                       gr_ibeg, gr_iend

  implicit none

  real, intent(IN) :: dt

  real, dimension(num_radius*2+1, NUMB_VAR) :: stencil
  real, dimension(NSYS_VAR) :: tempL, tempR

  integer :: i, j, dir
  integer :: ip, im, jp, jm
  integer :: s, si, sj
  integer :: ibeg, iend, jbeg, jend


  ibeg = gr_ibeg(XDIM)
  iend = gr_iend(XDIM)

  jbeg = gr_ibeg(YDIM)
  jend = gr_iend(YDIM)

  ! getting maximum speed for flux splitting
  call get_maxSpeed(gr_V)

  do dir = XDIM, NDIM

    select case(dir)
    case(XDIM)
      ip = num_radius
      im = num_radius
      jp = 0
      jm = 0
    case(YDIM)
      ip = 0
      im = 0
      jp = num_radius
      jm = num_radius
    end select

    do j = jbeg-1, jend+1
      do i = ibeg-1, iend+1

        ! make 1D stencil
        s = 0
        do sj = j-jm, j+jp
          do si = i-im, i+ip
            s = s + 1
            stencil(s, :) = gr_V(:, si, sj)
          end do
        end do

        call num_spatial_method(dt, num_radius, stencil, tempL, tempR, dir)
        gr_fL(:,i,j,dir) = tempL(:)
        gr_fR(:,i,j,dir) = tempR(:)

      end do !i
    end do   !j

  end do     !dir

  call soln_getFlux()


  return
end subroutine soln_spatial
