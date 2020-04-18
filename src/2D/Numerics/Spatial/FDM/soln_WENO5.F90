subroutine soln_WENO5(dt, radius, Wl, Wr, reconL, reconR)
  ! INPUT:
  !   dt    : infinitesimal time interval
  !   radius: radius of stencil
  !   Wl    : left biased (projected)stencil;  i-(r+1):i+(r-1)
  !   Wr    : right biased (projected)stencil; i-r:i+r
  ! OUTPUT:
  !   reconL: reconstructed vector, left of i-h
  !   reconR: reconstructed vector, right of i-h
  ! EXAMPLE:
  !   |         |         | reconL>>|<<reconR |         |         |
  !   |         |         |         |         |         |         |
  !   |         |         |         |    i    |         |         |
  !   +---------+---------+---------+---------+---------+---------+
  !   |<---------------------  Wl  -------------------->|
  !             |<---------------------  Wr  -------------------->|

#include "definition.h"

  use WENO,     only: betas
  use sim_data, only: sim_WENO, sim_WENeps, sim_mval

  implicit none

  real, intent(IN) :: dt
  integer, intent(IN) :: radius
  real, dimension(5, NSYS_VAR), intent(IN)  :: Wl, Wr
  real, dimension(NSYS_VAR),    intent(OUT) :: reconL, reconR

  real, dimension(3) :: lin_w
  real, dimension(3) :: smth_ind_L, smth_ind_R
  real, dimension(3, 3) :: ENO_coeff
  real, dimension(3, 2) :: ENO_intp, nonLin_w

  integer :: var, s, r
  real    :: w_norm

  ! (left) -->
  ! <--(right)
  lin_w(:)       = (/ .3,   .6,  .1  /)
  ENO_coeff(:,1) = (/ -1.,  5.,   2. /)
  ENO_coeff(:,2) = (/  2.,  5.,  -1. /)
  ENO_coeff(:,3) = (/ 11., -7.,   2. /)

  ENO_coeff = ENO_coeff/6.


  do var = 1, NSYS_VAR

    !get smoothness indicators
    call betas(Wl(:,var), radius, smth_ind_L)
    call betas(Wr(:,var), radius, smth_ind_R)

    ! compute non-linear weights
    do s = 1, 3
      r = 4 - s
      if (sim_WENO == '5') then
        !WENO-JS
        nonLin_w(s,1) = lin_w(r)/(sim_WENeps + smth_ind_L(s))**sim_mval  ! left : i->iph
        nonLin_w(s,2) = lin_w(s)/(sim_WENeps + smth_ind_R(s))**sim_mval  ! right: i->imh
      else if (sim_WENO == 'Z') then
        !WENO-Z
        nonLin_w(s,1) = lin_w(r)*(1.+&   ! left : i->iph
             abs(smth_ind_L(3)-smth_ind_L(1))/(sim_WENeps+smth_ind_L(s)))**sim_mval
        nonLin_w(s,2) = lin_w(s)*(1.+&   ! right: i->imh
             abs(smth_ind_R(3)-smth_ind_R(1))/(sim_WENeps+smth_ind_R(s)))**sim_mval
      else
        call abort_slug("unrecognized sim_WENO")
      end if
    end do

    ! normalize non-linear weights
    w_norm = SUM(nonLin_w(:,1))
    nonLin_w(:,1) = nonLin_w(:,1)/w_norm
    w_norm = SUM(nonLin_w(:,2))
    nonLin_w(:,2) = nonLin_w(:,2)/w_norm

    ! calculate ENO interpolations
    do s = 1, 3
      r = 4 - s
      ENO_intp(s,1) = dot_product(Wl(s+2:s:-1,var), ENO_coeff(:,r))
      ENO_intp(s,2) = dot_product(Wr(s:s+2   ,var), ENO_coeff(:,s))
    end do

    ! take convex combination of ENO statets
    !           |
    !  reconL>> | <<reconR
    !      -----+-----
    !          imh
    reconL(var) = dot_product(nonLin_w(:,1), ENO_intp(:,1))
    reconR(var) = dot_product(nonLin_w(:,2), ENO_intp(:,2))

  end do  ! end var


  return
end subroutine soln_WENO5
