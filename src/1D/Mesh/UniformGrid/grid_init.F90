subroutine grid_init()

#include "definition.h"

  use grid_data
  use sim_data
  use block_data, only : bl_i

  implicit none

  integer :: i

  ! the first and the last interior cell index
  allocate(gr_ibeg(NDIM)); allocate(gr_iend(NDIM))
  allocate(gr_i0(  NDIM)); allocate(gr_imax(NDIM))

  gr_ibeg(XDIM) = gr_ngc + 1
  gr_i0(XDIM)   = 1

  gr_iend(XDIM) = gr_ngc + gr_nx
  gr_imax(XDIM) = gr_iend(XDIM) + gr_ngc

  ! allocate cell coordinates
  allocate(gr_xCoord(gr_nx+2*gr_ngc)); gr_xCoord = 0.0

  ! grid delta
  gr_dx = (gr_xend(XDIM) - gr_xbeg(XDIM))/gr_glb_nx

  do i = gr_i0(XDIM), gr_imax(XDIM)
    ! generate x-coordinate
    gr_xCoord(i) = (real(i-gr_ngc + (bl_i-1)*gr_nx)-0.5)*gr_dx + gr_xbeg(XDIM)
  end do

  ! allocate grid variables
  allocate(gr_U(NSYS_VAR, gr_imax(XDIM))); gr_U = 0.
  allocate(gr_V(NUMB_VAR, gr_imax(XDIM))); gr_V = 0.

  ! allocate Riemann states
  allocate(gr_vL(NUMB_VAR, gr_imax(XDIM), NDIM)); gr_vL = 0.
  allocate(gr_vR(NUMB_VAR, gr_imax(XDIM), NDIM)); gr_vR = 0.

  ! allocate interfacial fluxes
  allocate(gr_fL(NSYS_VAR, gr_imax(XDIM), NDIM)); gr_fL = 0.
  allocate(gr_fR(NSYS_VAR, gr_imax(XDIM), NDIM)); gr_fR = 0.

  ! allocate grid fluxes
  allocate(gr_flux(NSYS_VAR, gr_imax(XDIM), NDIM)); gr_flux = 0.

  ! allocate maximum speed
  allocate(gr_maxSpeed(NUMB_WAVE,NDIM)); gr_maxSpeed = 0.

  return
end subroutine grid_init
