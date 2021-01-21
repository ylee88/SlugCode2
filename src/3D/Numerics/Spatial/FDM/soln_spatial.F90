subroutine soln_spatial(dt, prim, cons, flux)

#include "definition.h"

  use num_data,     only: num_radius, num_spatial_method
  use eigensystem,  only: left_eigenvectors, right_eigenvectors
  use grid_data,    only: gr_flux,          &  ! output
                          gr_ibeg, gr_iend, &
                          gr_i0, gr_imax,   &
                          gr_maxSpeed

  implicit none

  real, intent(IN) :: dt

  real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)), intent(IN) :: prim
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)), intent(IN) :: cons
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM), NDIM), intent(IN) :: flux

  real, dimension(NUMB_VAR, 2*num_radius+2) :: Vstncl
  real, dimension(NSYS_VAR, 2*num_radius+2) :: Fstncl, Ustncl
  real, dimension(2*num_radius+1, NSYS_VAR) :: Wstncl_L, Wstncl_R
  real, dimension(NSYS_VAR, NUMB_WAVE) :: leig, reig
  real, dimension(NUMB_VAR) :: Vimh
  real, dimension(NSYS_VAR) :: tempL, tempR

  real :: alpha

  integer :: i, j, k, dir, var
  integer :: ip, im, jp, jm, kp, km
  integer :: s, si, sj, sk
  integer :: ibeg, iend, jbeg, jend, kbeg, kend
  integer :: i0, imax, j0, jmax, k0, kmax

  logical :: conservative

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

  ! we need conservative eigenvectors
  conservative = .true.

  ! getting maximum speed for flux splitting
  call get_maxSpeed(prim)

  do dir = XDIM, NDIM

    alpha = MAXVAL(gr_maxSpeed(:,dir))

    select case(dir)
    case(XDIM)
      ip = num_radius
      im = num_radius+1  ! we need left biased stencil; i-(r+1):i+r
      jp = 0
      jm = 0
      kp = 0
      km = 0
    case(YDIM)
      ip = 0
      im = 0
      jp = num_radius
      jm = num_radius+1  ! we need left biased stencil; i-(r+1):i+r
      kp = 0
      km = 0
    case(ZDIM)
      ip = 0
      im = 0
      jp = 0
      jm = 0
      kp = num_radius
      km = num_radius+1  ! we need left biased stencil; i-(r+1):i+r
    end select

    do k = kbeg, kend+1
      do j = jbeg, jend+1
        do i = ibeg, iend+1

          ! make 1D stencil
          s = 0
          do sk = k-km, k+kp
            do sj = j-jm, j+jp
              do si = i-im, i+ip
                s = s + 1
                Vstncl(:, s) = prim(:, si, sj, sk)
                Ustncl(:, s) = cons(:, si, sj, sk)
                Fstncl(:, s) = flux(:, si, sj, sk, dir)
              end do
            end do
          end do

          ! below block is a flux splitting; doi.org/10.1016/j.jcp.2010.04.013
          ! we project stencil to i-h
          Vimh = 0.5*(Vstncl(:,num_radius+1) + Vstncl(:,num_radius+2))
          call left_eigenvectors (Vimh(:), conservative, leig, dir)
          call right_eigenvectors(Vimh(:), conservative, reig, dir)
          do s = 1, 2*num_radius+1
            do var = 1, NSYS_VAR
              Wstncl_L(s, var) = 0.5*dot_product( leig(:,var), Fstncl(:,  s) + alpha*Ustncl(:,  s) )
              Wstncl_R(s, var) = 0.5*dot_product( leig(:,var), Fstncl(:,s+1) - alpha*Ustncl(:,s+1) )
            end do
          end do

          ! (projected)flux reconstruction
          !             |
          !    tempL >> | << tmpR
          !             |
          !       ------+------
          !            i-h
          call num_spatial_method(dt, num_radius, Wstncl_L, Wstncl_R, tempL, tempR)

          ! back projection. result is Fimh
          do var = 1, NSYS_VAR
            gr_flux(var, i, j, k, dir) = dot_product( tempL(:) + tempR(:), reig(var,:))
          end do

        end do !i
      end do   !j
    end do     !k

  end do     !dir

  return
end subroutine soln_spatial
