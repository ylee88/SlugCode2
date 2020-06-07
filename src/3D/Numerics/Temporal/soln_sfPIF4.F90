subroutine soln_sfPIF4(dt)

#include "definition.h"

  use grid_data, only: gr_V,               &
                       gr_i0, gr_imax,     &
                       gr_ibeg, gr_iend,   &
                       gr_dx, gr_dy, gr_dz
  use num_data, only: num_radius
  use primconsflux
  use sfPIF
  use num_diff

  implicit none
  real, intent(IN) :: dt
  integer :: i, j, k
  integer :: ibeg, iend, jbeg, jend, kbeg, kend
  integer :: i0, imax, j0, jmax, k0, kmax

  real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)) :: V
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)) :: U, F, G, H
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM) :: Flux   ! solution for sfPIF

  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)) :: div_gr

  real, dimension(NSYS_VAR) :: Ui, Ux, Uxx, Uy, Uyy, Uz, Uzz
  real, dimension(NSYS_VAR) :: Fx, Gy, Hz
  real, dimension(NSYS_VAR) :: Fxt, Gyt, Hzt
  real, dimension(NSYS_VAR) :: Ft, Ftt, Fttt, Gt, Gtt, Gttt, Ht, Htt, Httt
  real, dimension(NSYS_VAR) :: div, divx, divy, divz, divt

  real, dimension(NSYS_VAR) :: Uxy, Uyz, Uxz
  real, dimension(NSYS_VAR) :: Fxtt, Gytt, Hztt

  real, dimension(NSYS_VAR) :: divxx, divyy, divzz, divxy, divxz, divyz
  real, dimension(NSYS_VAR) :: Fxxt, Fxyt, Fxzt, Gyxt, Gyyt, Gyzt, Hzxt, Hzyt, Hzzt
  real, dimension(NSYS_VAR) :: divtx, divty, divtz, divtt

  ! for readability
  ibeg = gr_ibeg(XDIM)
  iend = gr_iend(XDIM)
  jbeg = gr_ibeg(YDIM)
  jend = gr_iend(YDIM)
  kbeg = gr_ibeg(ZDIM)
  kend = gr_iend(ZDIM)

  i0   = gr_i0(XDIM)
  imax = gr_imax(XDIM)
  j0   = gr_i0(YDIM)
  jmax = gr_imax(YDIM)
  k0   = gr_i0(ZDIM)
  kmax = gr_imax(ZDIM)

  V = gr_V
  ! initial data for spatial recon/intp
  do k = k0, kmax
    do j = j0, jmax
      do i = i0, imax
        call prim2cons(V(:,i,j,k), U(:,i,j,k))
        call prim2flux(V(:,i,j,k), F(:,i,j,k), XDIM)
        call prim2flux(V(:,i,j,k), G(:,i,j,k), YDIM)
        call prim2flux(V(:,i,j,k), H(:,i,j,k), ZDIM)
      end do
    end do
  end do

  do k = kbeg-(num_radius+1+2), kend+(num_radius+1+2)
    do j = jbeg-(num_radius+1+2), jend+(num_radius+1+2)
      do i = ibeg-(num_radius+1+2), iend+(num_radius+1+2)
        Fx = diff1(F(:, i-2:i+2, j, k), 5, gr_dx)
        Gy = diff1(G(:, i, j-2:j+2, k), 5, gr_dy)
        Hz = diff1(H(:, i, j, k-2:k+2), 5, gr_dz)

        div_gr(:,i,j,k) = Fx + Gy + Hz
      end do
    end do
  end do

  ! initialize solution vector
  Flux = 0.

  ! calculate time averaging fluxes
  do k = kbeg-(num_radius+1), kend+(num_radius+1)
    do j = jbeg-(num_radius+1), jend+(num_radius+1)
      do i = ibeg-(num_radius+1), iend+(num_radius+1)

        Ui = U(:,i,j,k)
        div = div_gr(:,i,j,k)

        ! second order
        Ft = -get_Jv(Ui, div, dt, XDIM)
        Gt = -get_Jv(Ui, div, dt, YDIM)
        Ht = -get_Jv(Ui, div, dt, ZDIM)

        Flux(:, i, j, k, XDIM) = F(:, i, j, k) + dt/2.*Ft(:)
        Flux(:, i, j, k, YDIM) = G(:, i, j, k) + dt/2.*Gt(:)
        Flux(:, i, j, k, ZDIM) = H(:, i, j, k) + dt/2.*Ht(:)

        ! start to building third order term
        Ux = diff1(U(:, i-2:i+2,       j,       k), 5, gr_dx)
        Uy = diff1(U(:,       i, j-2:j+2,       k), 5, gr_dy)
        Uz = diff1(U(:,       i,       j, k-2:k+2), 5, gr_dz)

        Uxx = diff2(U(:, i-2:i+2,       j,       k), 5, gr_dx)
        Uyy = diff2(U(:,       i, j-2:j+2,       k), 5, gr_dy)
        Uzz = diff2(U(:,       i,       j, k-2:k+2), 5, gr_dz)

        ! building divt
        divx = diff1(div_gr(:, i-2:i+2,       j,       k), 5, gr_dx)
        divy = diff1(div_gr(:,       i, j-2:j+2,       k), 5, gr_dy)
        divz = diff1(div_gr(:,       i,       j, k-2:k+2), 5, gr_dz)

        Fxt = -get_Hvw(Ui, Ux, div, dt, XDIM) - get_Jv(Ui, divx, dt, XDIM)
        Gyt = -get_Hvw(Ui, Uy, div, dt, YDIM) - get_Jv(Ui, divy, dt, YDIM)
        Hzt = -get_Hvw(Ui, Uz, div, dt, ZDIM) - get_Jv(Ui, divz, dt, ZDIM)

        ! finalizing divt
        divt = Fxt + Gyt + Hzt

        ! finalizing third order terms
        Ftt = get_Hvw(Ui, div, div, dt, XDIM) - get_Jv(Ui, divt, dt, XDIM)
        Gtt = get_Hvw(Ui, div, div, dt, YDIM) - get_Jv(Ui, divt, dt, YDIM)
        Htt = get_Hvw(Ui, div, div, dt, ZDIM) - get_Jv(Ui, divt, dt, ZDIM)

        Flux(:, i, j, k, XDIM) = Flux(:, i, j, k, XDIM) + dt**2/6.*Ftt(:)
        Flux(:, i, j, k, YDIM) = Flux(:, i, j, k, YDIM) + dt**2/6.*Gtt(:)
        Flux(:, i, j, k, ZDIM) = Flux(:, i, j, k, ZDIM) + dt**2/6.*Htt(:)


        ! start to building fourth order term
        Uxx = diff2(U(:, i-2:i+2, j, k), 5, gr_dx)
        Uyy = diff2(U(:, i, j-2:j+2, k), 5, gr_dy)
        Uzz = diff2(U(:, i, j, k-2:k+2), 5, gr_dz)

        Uxy = diffxy(U(:, i-2:i+2, j-2:j+2, k), 5, gr_dx, gr_dy)
        Uyz = diffxy(U(:, i, j-2:j+2, k-2:k+2), 5, gr_dy, gr_dz)
        Uxz = diffxy(U(:, i-2:i+2, j, k-2:k+2), 5, gr_dx, gr_dz)

        divxx = diff2(div_gr(:, i-2:i+2,       j,       k), 5, gr_dx)
        divyy = diff2(div_gr(:,       i, j-2:j+2,       k), 5, gr_dy)
        divzz = diff2(div_gr(:,       i,       j, k-2:k+2), 5, gr_dz)

        divxy = diffxy(div_gr(:, i-2:i+2, j-2:j+2,       k), 5, gr_dx, gr_dy)
        divyz = diffxy(div_gr(:,       i, j-2:j+2, k-2:k+2), 5, gr_dy, gr_dz)
        divxz = diffxy(div_gr(:, i-2:i+2,       j, k-2:k+2), 5, gr_dx, gr_dz)

        Fxxt = -get_Dvwx(Ui, Ux, Ux, div, dt, XDIM) - get_Hvw(Ui, Uxx, div, dt, XDIM) -2.*get_Hvw(Ui, Ux, divx, dt, XDIM) &
               -get_Jv(Ui, divxx, dt, XDIM)
        Fxyt = -get_Dvwx(Ui, Ux, Uy, div, dt, XDIM) - get_Hvw(Ui, Uxy, div, dt, XDIM) - get_Hvw(Ui, Ux, divy, dt, XDIM) &
               -get_Hvw(Ui, Uy, divx, dt, XDIM) - get_Jv(Ui, divxy, dt, XDIM)
        Fxzt = -get_Dvwx(Ui, Ux, Uz, div, dt, XDIM) - get_Hvw(Ui, Uxz, div, dt, XDIM) - get_Hvw(Ui, Ux, divz, dt, XDIM) &
               -get_Hvw(Ui, Uz, divx, dt, XDIM) - get_Jv(Ui, divxz, dt, XDIM)

        Gyxt = -get_Dvwx(Ui, Uy, Ux, div, dt, YDIM) - get_Hvw(Ui, Uxy, div, dt, YDIM) - get_Hvw(Ui, Uy, divx, dt, YDIM) &
               -get_Hvw(Ui, Ux, divy, dt, YDIM) - get_Jv(Ui, divxy, dt, YDIM)
        Gyyt = -get_Dvwx(Ui, Uy, Uy, div, dt, YDIM) - get_Hvw(Ui, Uyy, div, dt, YDIM) -2.*get_Hvw(Ui, Uy, divy, dt, YDIM) &
               -get_Jv(Ui, divyy, dt, YDIM)
        Gyzt = -get_Dvwx(Ui, Uy, Uz, div, dt, YDIM) - get_Hvw(Ui, Uyz, div, dt, YDIM) - get_Hvw(Ui, Uy, divz, dt, YDIM) &
               -get_Hvw(Ui, Uz, divy, dt, YDIM) - get_Jv(Ui, divyz, dt, YDIM)

        Hzxt = -get_Dvwx(Ui, Uz, Ux, div, dt, ZDIM) - get_Hvw(Ui, Uxz, div, dt, ZDIM) - get_Hvw(Ui, Uz, divx, dt, ZDIM) &
               -get_Hvw(Ui, Ux, divz, dt, ZDIM) - get_Jv(Ui, divxz, dt, ZDIM)
        Hzyt = -get_Dvwx(Ui, Uz, Uy, div, dt, ZDIM) - get_Hvw(Ui, Uyz, div, dt, ZDIM) - get_Hvw(Ui, Uz, divy, dt, ZDIM) &
               -get_Hvw(Ui, Uy, divz, dt, ZDIM) - get_Jv(Ui, divyz, dt, ZDIM)
        Hzzt = -get_Dvwx(Ui, Uz, Uz, div, dt, ZDIM) - get_Hvw(Ui, Uzz, div, dt, ZDIM) -2.*get_Hvw(Ui, Uz, divz, dt, ZDIM) &
               -get_Jv(Ui, divzz, dt, ZDIM)

        divtx = Fxxt + Gyxt + Hzxt
        divty = Fxyt + Gyyt + Hzyt
        divtz = Fxzt + Gyzt + Hzzt

        Fxtt = get_Dvwx(Ui, div, Ux, div, dt, XDIM) + 2.*get_Hvw(Ui, div, divx, dt, XDIM)  &
               -get_Hvw(Ui, Ux, divt, dt, XDIM) - get_Jv(Ui, divtx, dt, XDIM)
        Gytt = get_Dvwx(Ui, div, Uy, div, dt, YDIM) + 2.*get_Hvw(Ui, div, divy, dt, YDIM)  &
               -get_Hvw(Ui, Uy, divt, dt, YDIM) - get_Jv(Ui, divty, dt, YDIM)
        Hztt = get_Dvwx(Ui, div, Uz, div, dt, ZDIM) + 2.*get_Hvw(Ui, div, divz, dt, ZDIM)  &
               -get_Hvw(Ui, Uz, divt, dt, ZDIM) - get_Jv(Ui, divtz, dt, ZDIM)

        ! finalizing divtt
        divtt = Fxtt + Gytt + Hztt

        Fttt = -get_Dvwx(Ui, div, div, div, dt, XDIM) +3.*get_Hvw(Ui, div, divt, dt, XDIM) - get_Jv(Ui, divtt, dt, XDIM)
        Gttt = -get_Dvwx(Ui, div, div, div, dt, YDIM) +3.*get_Hvw(Ui, div, divt, dt, YDIM) - get_Jv(Ui, divtt, dt, YDIM)
        Httt = -get_Dvwx(Ui, div, div, div, dt, ZDIM) +3.*get_Hvw(Ui, div, divt, dt, ZDIM) - get_Jv(Ui, divtt, dt, ZDIM)

        Flux(:, i, j, k, XDIM) = Flux(:, i, j, k, XDIM) + dt**3/24.*Fttt(:)
        Flux(:, i, j, k, YDIM) = Flux(:, i, j, k, YDIM) + dt**3/24.*Gttt(:)
        Flux(:, i, j, k, ZDIM) = Flux(:, i, j, k, ZDIM) + dt**3/24.*Httt(:)

      end do
    end do
  end do


  ! spatial recon/intp
  call soln_spatial(dt, V, U, Flux)

  return

end subroutine soln_sfPIF4
