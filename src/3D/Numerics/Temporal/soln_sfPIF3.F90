subroutine soln_sfPIF3(dt)

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
  real, dimension(NSYS_VAR) :: Ft, Ftt, Gt, Gtt, Ht, Htt
  real, dimension(NSYS_VAR) :: div, divx, divy, divz, divt

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

      end do
    end do
  end do


  ! spatial recon/intp
  call soln_spatial(dt, V, U, Flux)

  return

end subroutine soln_sfPIF3
