subroutine soln_WENO(dt, radius, V, reconL, reconR, dir)
  ! INPUT:
  !   dt    : infinitesimal time interval
  !   radius: radius of stencil
  !   V     : (2R+1, NUMB_VAR) primitive vector
  !   dir   : direction
  ! OUTPUT:
  !   reconL: (NSYS_VAR) vector at i-h
  !   reconR: (NSYS_VAR) vector at i+h

#include "definition.h"

  use WENO,     only: betas
  use sim_data, only: sim_WENO, sim_WENeps, sim_mval
  use char_limiting

  implicit none

  real, intent(IN) :: dt
  integer, intent(IN) :: radius, dir
  real, dimension(2*radius+1, NUMB_VAR), intent(IN)    :: V
  real, dimension(NSYS_VAR),             intent(INOUT) :: reconL, reconR

  real, dimension(2*radius+1, NSYS_VAR) :: stencil_L, stencil_R
  real, dimension(radius+1) :: lin_w
  real, dimension(radius+1) :: smth_ind_L, smth_ind_R
  real, dimension(radius+1, radius+1) :: ENO_coeff
  real, dimension(radius+1, 2) :: ENO_intp, nonLin_w

  real, dimension(NSYS_VAR) :: tempL, tempR

  integer :: var, s, r
  real    :: w_norm

  ! (left) -->
  ! <--(right)
  select case(radius+1)
  case(3)
    lin_w(:)       = (/ .3,   .6,  .1  /)
    ENO_coeff(:,1) = (/ -1.,  5.,   2. /)
    ENO_coeff(:,2) = (/  2.,  5.,  -1. /)
    ENO_coeff(:,3) = (/ 11., -7.,   2. /)

    ENO_coeff = ENO_coeff/6.
  case default
    call abort_slug("unsupported WENO radius")
  end select


  ! char limiting for FVM
  ! flux splitting for FDM
  call char_proj(radius, V(:,:), stencil_L(:,:), stencil_R(:,:), dir)

  do var = 1, NSYS_VAR

    !get smoothness indicators
    call betas(stencil_L(:,var), radius, smth_ind_L)
    call betas(stencil_R(:,var), radius, smth_ind_R)

    !compute non-linear weights
    do s = 1, radius+1
      r = (radius+1)+1-s
      if (sim_WENO == '5') then
        !WENO-JS
        nonLin_w(s,1) = lin_w(s)/(sim_WENeps + smth_ind_L(s))**sim_mval  ! left : i->imh
        nonLin_w(s,2) = lin_w(r)/(sim_WENeps + smth_ind_R(s))**sim_mval  ! right: i->iph
      else if (sim_WENO == 'Z') then
        !WENO-Z
        nonLin_w(s,1) = lin_w(s)*(1.+&   ! left : i->imh
             abs(smth_ind_L(radius+1)-smth_ind_L(1))/(sim_WENeps+smth_ind_L(s)))**sim_mval
        nonLin_w(s,2) = lin_w(r)*(1.+&   ! right: i->iph
             abs(smth_ind_R(radius+1)-smth_ind_R(1))/(sim_WENeps+smth_ind_R(s)))**sim_mval
      else
        call abort_slug("unrecognized sim_WENO")
      end if
    end do

    w_norm = SUM(nonLin_w(:,1))
    nonLin_w(:,1) = nonLin_w(:,1)/w_norm
    w_norm = SUM(nonLin_w(:,2))
    nonLin_w(:,2) = nonLin_w(:,2)/w_norm

    !calculate ENO interpolations
    do s = 1, radius+1
      r = (radius+1)+1 - s
      ENO_intp(s,1) = dot_product(stencil_L(s:s+2   ,var), ENO_coeff(:,s))
      ENO_intp(s,2) = dot_product(stencil_R(s+2:s:-1,var), ENO_coeff(:,r))
    end do

    !take convex combination of ENO statets
    tempL(var) = dot_product(nonLin_w(:,1), ENO_intp(:,1))
    tempR(var) = dot_product(nonLin_w(:,2), ENO_intp(:,2))

  end do

  ! project back from char field
  call char_back_proj(radius, V(:,:), tempL(:), tempR(:), &
                      reconL(:), reconR(:), dir)

  return
end subroutine soln_WENO
