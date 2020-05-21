subroutine soln_gpWENO7(dt, radius, Wl, Wr, reconL, reconR)
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

  use gp_data, only: gp_linW, gp_Zk
  use WENO,    only: gp_betas
  use sim_data, only: sim_WENeps, sim_mval

  implicit none

  ! this will determine the order of current subroutine.
  integer, parameter :: R = 3

  real, intent(IN) :: dt
  integer, intent(IN) :: radius
  real, dimension(2*R+1, NSYS_VAR), intent(IN)  :: Wl, Wr
  real, dimension(NSYS_VAR),    intent(OUT) :: reconL, reconR

  real, dimension(R+1) :: smth_ind_L, smth_ind_R

  real, dimension(R+1, 2) :: wbar, weights
  real, dimension(R+1)    :: vMk, vPk

  integer :: var, k, N, M
  real    :: sum_wbar

  M = R+1
  N = 2*R+1

  do var = 1, NSYS_VAR

    call gp_betas(Wl(:,var), R, smth_ind_L)
    call gp_betas(Wr(:,var), R, smth_ind_R)

    do k = 1, M
      wbar(k, 1) = gp_linW(2,k)/(sim_WENeps + smth_ind_L(k))**sim_mval  ! left  : i->iph
      wbar(k, 2) = gp_linW(1,k)/(sim_WENeps + smth_ind_R(k))**sim_mval  ! right : i->imh
    end do

    !normalize weights
    sum_wbar = SUM(wbar(:,1))
    weights(:,1) = wbar(:,1)/sum_wbar
    sum_wbar = SUM(wbar(:,2))
    weights(:,2) = wbar(:,2)/sum_wbar

    !get GP predictions on ENO stencils
    do k = 1, M
      !characteristics moving (+) to interface
      vPk(k) = dot_product(gp_Zk(2,1:M,k), Wl(k:k+R,var))
      !characteristics moving (-) to interface
      vMk(k) = dot_product(gp_Zk(1,1:M,k), Wr(k:k+R,var))
    end do

    ! take convex combination of ENO statets
    !           |
    !  reconL>> | <<reconR
    !      -----+-----
    !          imh
    reconL(var) = dot_product(weights(:,1), vPk(:))
    reconR(var) = dot_product(weights(:,2), vMk(:))
  end do

end subroutine soln_gpWENO7
