subroutine soln_getFlux()

#include "definition.h"

  use grid_data, only: gr_ibeg, gr_iend, &
                       gr_fL, gr_fR,     &
                       gr_flux

  implicit none

  integer :: i, j, dir
  integer :: im, jm
  integer :: ibeg, iend, jbeg, jend

  ibeg = gr_ibeg(XDIM)
  iend = gr_iend(XDIM)

  jbeg = gr_ibeg(YDIM)
  jend = gr_iend(YDIM)

  ! get flux
  do dir = XDIM, NDIM

    select case(dir)
    case(XDIM)
      im = 1
      jm = 0
    case(YDIM)
      im = 0
      jm = 1
    end select

    do j = jbeg, jend+1
      do i = ibeg, iend+1
        ! flux_imh
        gr_flux(:,i,j,dir) = gr_fL(:,i,j,dir) + gr_fR(:,i-im,j-jm,dir)
      end do
    end do

  end do

end subroutine soln_getFlux
