subroutine soln_sfPIF3(dt)

#include "definition.h"

  use grid_data, only: gr_V,             &
                       gr_i0, gr_imax,   &
                       gr_ibeg, gr_iend, &
                       gr_dx, gr_dy
  use num_data, only: num_radius
  use primconsflux
  use sfPIF
  use num_diff

  implicit none
  real, intent(IN) :: dt
  integer :: i, j
  integer :: ibeg, iend, jbeg, jend
  integer :: i0, imax, j0, jmax

  real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: V
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: U, F, G
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: Flux   ! solution for sfPIF

  real, dimension(NSYS_VAR) :: Ui, Ux, Uxx, Uy, Uyy
  real, dimension(NSYS_VAR) :: Fx, Fxy, Fxx, Gy, Gxy, Gyy
  real, dimension(NSYS_VAR) :: Fxt, Gyt
  real, dimension(NSYS_VAR) :: Ft, Ftt, Gt, Gtt
  real, dimension(NSYS_VAR) :: div, divx, divy, divt

  ! for readability
  ibeg = gr_ibeg(XDIM)
  iend = gr_iend(XDIM)
  jbeg = gr_ibeg(YDIM)
  jend = gr_iend(YDIM)

  i0   = gr_i0(XDIM)
  imax = gr_imax(XDIM)
  j0   = gr_i0(YDIM)
  jmax = gr_imax(YDIM)

  V = gr_V
  ! initial data for spatial recon/intp
  do j = j0, jmax
    do i = i0, imax
      call prim2cons(V(:,i,j), U(:,i,j))
      call prim2flux(V(:,i,j), F(:,i,j), XDIM)
      call prim2flux(V(:,i,j), G(:,i,j), YDIM)
    end do
  end do

  ! initialize solution vector
  Flux = 0.

  ! calculate time averaging fluxes
  do j = jbeg-(num_radius+1), jend+(num_radius+1)
    do i = ibeg-(num_radius+1), iend+(num_radius+1)

      Ui = U(:,i,j)

      Fx = diff1(F(:, i-2:i+2, j), 5, gr_dx)
      Gy = diff1(G(:, i, j-2:j+2), 5, gr_dy)

      div = Fx + Gy

      ! second order
      Ft = -get_Jv(Ui, div, dt, XDIM)
      Gt = -get_Jv(Ui, div, dt, YDIM)

      Flux(:, i, j, XDIM) = F(:, i, j)+ dt/2.*Ft(:)
      Flux(:, i, j, YDIM) = G(:, i, j)+ dt/2.*Gt(:)


      ! start to building third order term
      Ux = diff1(U(:, i-2:i+2,       j), 5, gr_dx)
      Uy = diff1(U(:,       i, j-2:j+2), 5, gr_dy)

      Uxx = diff2(U(:, i-2:i+2,       j), 5, gr_dx)
      Uyy = diff2(U(:,       i, j-2:j+2), 5, gr_dy)

      Fxx = diff2(F(:, i-2:i+2, j), 5, gr_dx)
      Gyy = diff2(G(:, i, j-2:j+2), 5, gr_dy)

      ! cross derivatives
      Fxy = diffxy(F(:, i-2:i+2, j-2:j+2), 5)
      Gxy = diffxy(G(:, i-2:i+2, j-2:j+2), 5)

      ! building divt
      divx = Fxx + Gxy
      divy = Gyy + Fxy

      Fxt = -get_Hvw(Ui, Ux, div, dt, XDIM) - get_Jv(Ui, divx, dt, XDIM)
      Gyt = -get_Hvw(Ui, Uy, div, dt, YDIM) - get_Jv(Ui, divy, dt, YDIM)

      ! finalizing divt
      divt = Fxt + Gyt

      ! finalizing third order terms
      Ftt = get_Hvw(Ui, div, div, dt, XDIM) - get_Jv(Ui, divt, dt, XDIM)
      Gtt = get_Hvw(Ui, div, div, dt, YDIM) - get_Jv(Ui, divt, dt, YDIM)

      Flux(:, i, j, XDIM) = Flux(:, i, j, XDIM) + dt**2/6.*Ftt(:)
      Flux(:, i, j, YDIM) = Flux(:, i, j, YDIM) + dt**2/6.*Gtt(:)

    end do
  end do


  ! spatial recon/intp
  call soln_spatial(dt, V, U, Flux)

  return

end subroutine soln_sfPIF3
