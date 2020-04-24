subroutine soln_sfPIF4(dt)

#include "definition.h"

  use grid_data, only: gr_V,           &
                       gr_i0, gr_imax, &
                       gr_ibeg, gr_iend
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

  real, dimension(NSYS_VAR) :: Ui, Ux, Uxx, Uy, Uyy, Uxy
  real, dimension(NSYS_VAR) :: Fx, Fxy, Fxx, Gy, Gxy, Gyy
  real, dimension(NSYS_VAR) :: Fxxx, Gyyy, Fxxy, Gxxy, Fxyy, Gxyy
  real, dimension(NSYS_VAR) :: Fxt, Gyt, Fxxt, Gxyt, Fxyt, Gyyt, Fxtt, Gytt
  real, dimension(NSYS_VAR) :: Ft, Ftt, Fttt, Gt, Gtt, Gttt
  real, dimension(NSYS_VAR) :: div, divx, divy, divt, divtx, divty, divtt
  real, dimension(NSYS_VAR) :: divxx, divyy, divxy

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

      Fx = diff1(F(:, i-2:i+2, j), 5, XDIM)
      Gy = diff1(G(:, i, j-2:j+2), 5, YDIM)

      div = Fx + Gy

      ! second order
      Ft = -get_Jv(Ui, div, dt, XDIM)
      Gt = -get_Jv(Ui, div, dt, YDIM)

      Flux(:, i, j, XDIM) = F(:, i, j)+ dt/2.*Ft(:)
      Flux(:, i, j, YDIM) = G(:, i, j)+ dt/2.*Gt(:)


      ! start to building third order term
      Ux = diff1(U(:, i-2:i+2,       j), 5, XDIM)
      Uy = diff1(U(:,       i, j-2:j+2), 5, YDIM)

      Uxx = diff2(U(:, i-2:i+2,       j), 5, XDIM)
      Uyy = diff2(U(:,       i, j-2:j+2), 5, YDIM)

      Fxx = diff2(F(:, i-2:i+2, j), 5, XDIM)
      Gyy = diff2(G(:, i, j-2:j+2), 5, YDIM)

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


      ! start to building fourth order term
      Fxxx = diff3(F(:,i-2:i+2, j), 5, XDIM)
      Gyyy = diff3(G(:,i, j-2:j+2), 5, YDIM)

      Fxxy = diffxxy(F(:, i-2:i+2, j-2:j+2), 5)
      Gxxy = diffxxy(G(:, i-2:i+2, j-2:j+2), 5)

      Fxyy = diffxyy(F(:, i-2:i+2, j-2:j+2), 5)
      Gxyy = diffxyy(G(:, i-2:i+2, j-2:j+2), 5)

      Uxy = diffxy(U(:, i-2:i+2, j-2:j+2), 5)

      divxx = Fxxx + Gxxy
      divyy = Fxyy + Gyyy
      divxy = Fxxy + Gxyy

      Fxxt = -get_Dvwx(Ui, Ux, Ux, div, dt, XDIM) - get_Hvw(Ui, Uxx, div, dt, XDIM) -2.*get_Hvw(Ui, Ux, divx, dt, XDIM) &
             -get_Jv(Ui, divxx, dt, XDIM)
      Gxyt = -get_Dvwx(Ui, Ux, Uy, div, dt, YDIM) - get_Hvw(Ui, Uxy, div, dt, YDIM) - get_Hvw(Ui, Uy, divx, dt, YDIM) &
             -get_Hvw(Ui, Ux, divy, dt, YDIM) - get_Jv(Ui, divxy, dt, YDIM)

      Fxyt = -get_Dvwx(Ui, Ux, Uy, div, dt, XDIM) - get_Hvw(Ui, Uxy, div, dt, XDIM) - get_Hvw(Ui, Ux, divy, dt, XDIM) &
             -get_Hvw(Ui, Uy, divx, dt, XDIM) - get_Jv(Ui, divxy, dt, XDIM)
      Gyyt = -get_Dvwx(Ui, Uy, Uy, div, dt, YDIM) - get_Hvw(Ui, Uyy, div, dt, YDIM) -2.*get_Hvw(Ui, Uy, divy, dt, YDIM) &
             -get_Jv(Ui, divyy, dt, YDIM)

      divtx = Fxxt + Gxyt
      divty = Fxyt + Gyyt

      Fxtt = get_Dvwx(Ui, div, Ux, div, dt, XDIM) + 2.*get_Hvw(Ui, div, divx, dt, XDIM)  &
             -get_Hvw(Ui, Ux, divt, dt, XDIM) - get_Jv(Ui, divtx, dt, XDIM)
      Gytt = get_Dvwx(Ui, div, Uy, div, dt, YDIM) + 2.*get_Hvw(Ui, div, divy, dt, YDIM)  &
             -get_Hvw(Ui, Uy, divt, dt, YDIM) - get_Jv(Ui, divty, dt, YDIM)

      ! finalizing divtt
      divtt = Fxtt + Gytt

      Fttt = -get_Dvwx(Ui, div, div, div, dt, XDIM) +3.*get_Hvw(Ui, div, divt, dt, XDIM) - get_Jv(Ui, divtt, dt, XDIM)
      Gttt = -get_Dvwx(Ui, div, div, div, dt, YDIM) +3.*get_Hvw(Ui, div, divt, dt, YDIM) - get_Jv(Ui, divtt, dt, YDIM)

      Flux(:, i, j, XDIM) = Flux(:, i, j, XDIM) + dt**3/24.*Fttt(:)
      Flux(:, i, j, YDIM) = Flux(:, i, j, YDIM) + dt**3/24.*Gttt(:)


    end do
  end do


  ! spatial recon/intp
  call soln_spatial(dt, V, U, Flux)

  return

end subroutine soln_sfPIF4
