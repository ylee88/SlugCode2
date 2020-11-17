subroutine soln_sfPIF3(dt)

#include "definition.h"

  use grid_data, only: gr_V,               &
                       gr_i0, gr_imax,     &
                       gr_ibeg, gr_iend,   &
                       gr_dx
  use num_data, only: num_radius
  use primconsflux
  use sfPIF
  use num_diff

  implicit none
  real, intent(IN) :: dt
  integer :: i
  integer :: ibeg, iend
  integer :: i0, imax

  real, dimension(NUMB_VAR, gr_imax(XDIM)) :: V
  real, dimension(NSYS_VAR, gr_imax(XDIM)) :: U, F
  real, dimension(NSYS_VAR, gr_imax(XDIM), NDIM) :: Flux   ! solution for sfPIF

  real, dimension(NSYS_VAR, gr_imax(XDIM)) :: div_gr

  real, dimension(NSYS_VAR) :: Ui, Ux, Uxx
  real, dimension(NSYS_VAR) :: Fx
  real, dimension(NSYS_VAR) :: Fxt
  real, dimension(NSYS_VAR) :: Ft, Ftt
  real, dimension(NSYS_VAR) :: div, divx, divt

  ! for readability
  ibeg = gr_ibeg(XDIM)
  iend = gr_iend(XDIM)

  i0   = gr_i0(XDIM)
  imax = gr_imax(XDIM)

  V = gr_V
  ! initial data for spatial recon/intp
  do i = i0, imax
    call prim2cons(V(:,i), U(:,i))
    call prim2flux(V(:,i), F(:,i), XDIM)
  end do

  do i = ibeg-(num_radius+1+2), iend+(num_radius+1+2)
    Fx = diff1(F(:, i-2:i+2), 5, gr_dx)

    div_gr(:,i) = Fx
  end do

  ! initialize solution vector
  Flux = 0.

  ! calculate time averaging fluxes
  do i = ibeg-(num_radius+1), iend+(num_radius+1)

    Ui = U(:,i)
    div = div_gr(:,i)

    ! second order
    Ft = -get_Jv(Ui, div, dt, XDIM)

    Flux(:, i, XDIM) = F(:, i) + dt/2.*Ft(:)

    ! start to building third order term
    Ux = diff1(U(:, i-2:i+2), 5, gr_dx)
    Uxx = diff2(U(:, i-2:i+2), 5, gr_dx)

    ! building divt
    divx = diff1(div_gr(:, i-2:i+2), 5, gr_dx)
    Fxt = -get_Hvw(Ui, Ux, div, dt, XDIM) - get_Jv(Ui, divx, dt, XDIM)

    ! finalizing divt
    divt = Fxt

    ! finalizing third order terms
    Ftt = get_Hvw(Ui, div, div, dt, XDIM) - get_Jv(Ui, divt, dt, XDIM)

    Flux(:, i, XDIM) = Flux(:, i, XDIM) + dt**2/6.*Ftt(:)

  end do


  ! spatial recon/intp
  call soln_spatial(dt, V, U, Flux)

  return

end subroutine soln_sfPIF3
