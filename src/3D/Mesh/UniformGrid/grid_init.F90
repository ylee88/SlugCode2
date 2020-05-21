subroutine grid_init()

#include "definition.h"

  use grid_data
  use sim_data
  use block_data, only : bl_i, bl_j, bl_k

  implicit none

  integer :: i

  ! the first and the last interior cell index
  allocate(gr_ibeg(NDIM)); allocate(gr_iend(NDIM))
  allocate(gr_i0(  NDIM)); allocate(gr_imax(NDIM))

  gr_ibeg(XDIM) = gr_ngc + 1
  gr_ibeg(YDIM) = gr_ngc + 1
  gr_ibeg(ZDIM) = gr_ngc + 1
  gr_i0(XDIM)   = 1
  gr_i0(YDIM)   = 1
  gr_i0(ZDIM)   = 1

  gr_iend(XDIM) = gr_ngc + gr_nx
  gr_iend(YDIM) = gr_ngc + gr_ny
  gr_iend(ZDIM) = gr_ngc + gr_nz
  gr_imax(XDIM) = gr_iend(XDIM) + gr_ngc
  gr_imax(YDIM) = gr_iend(YDIM) + gr_ngc
  gr_imax(ZDIM) = gr_iend(ZDIM) + gr_ngc

  ! allocate cell coordinates
  allocate(gr_xCoord(gr_nx+2*gr_ngc)); gr_xCoord = 0.0
  allocate(gr_yCoord(gr_ny+2*gr_ngc)); gr_yCoord = 0.0
  allocate(gr_zCoord(gr_nz+2*gr_ngc)); gr_zCoord = 0.0

  ! grid delta
  gr_dx = (gr_xend(XDIM) - gr_xbeg(XDIM))/gr_glb_nx
  gr_dy = (gr_xend(YDIM) - gr_xbeg(YDIM))/gr_glb_ny
  gr_dz = (gr_xend(ZDIM) - gr_xbeg(ZDIM))/gr_glb_nz

  do i = gr_i0(XDIM), gr_imax(XDIM)
    ! generate x-coordinate
    gr_xCoord(i) = (real(i-gr_ngc + (bl_i-1)*gr_nx)-0.5)*gr_dx + gr_xbeg(XDIM)
  end do

  do i = gr_i0(YDIM), gr_imax(YDIM)
    ! generate y-coordinates
    gr_yCoord(i) = (real(i-gr_ngc + (bl_j-1)*gr_ny)-0.5)*gr_dy + gr_xbeg(YDIM)
  end do

  do i = gr_i0(ZDIM), gr_imax(ZDIM)
    ! generate z-coordinates
    gr_zCoord(i) = (real(i-gr_ngc + (bl_k-1)*gr_nz)-0.5)*gr_dz + gr_xbeg(ZDIM)
  end do

  ! allocate grid variables
  allocate(gr_U(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM))); gr_U = 0.
  allocate(gr_V(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM))); gr_V = 0.

  ! allocate Riemann states
  allocate(gr_vL(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM)); gr_vL = 0.
  allocate(gr_vR(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM)); gr_vR = 0.

  ! allocate interfacial fluxes
  allocate(gr_fL(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM)); gr_fL = 0.
  allocate(gr_fR(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM)); gr_fR = 0.

  ! allocate grid fluxes
  allocate(gr_flux(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM)); gr_flux = 0.

  ! allocate maximum speed
  allocate(gr_maxSpeed(NUMB_WAVE,NDIM)); gr_maxSpeed = 0.

  return
end subroutine grid_init
