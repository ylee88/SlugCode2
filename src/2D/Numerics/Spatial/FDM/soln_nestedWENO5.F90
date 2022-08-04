subroutine soln_nestedWENO5(dt, radius, Wl, Wr, reconL, reconR)
  ! DESCRIPTION:
  !   Nested, multi-resolution WENO scheme. see https://doi.org/10.1016/j.jcp.2020.110006
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

  use WENO,     only: Leg_betas

  implicit none

  real, intent(IN) :: dt
  integer, intent(IN) :: radius
  real, dimension(:, :), intent(IN)  :: Wl, Wr
  real, dimension(NSYS_VAR),    intent(OUT) :: reconL, reconR

  real, dimension(3) :: coeff_q1L, coeff_q1R
  real, dimension(5) :: coeff_q2L, coeff_q2R
  real, dimension(5) :: tmpCoeffL, tmpCoeffR

  real, dimension(2) :: lin_w
  real, dimension(2) :: smth_indL, smth_indR
  real, dimension(2) :: nonLin_wL, nonLin_wR

  real :: q0L, q0R, q1L, q1R, q2L, q2R
  real :: p11L, p11R, p12L, p12R, p22L, p22R


  integer :: var

  integer :: i

  i = 3

  lin_w(1) = 0.05     ! gamma_01, gamma_12
  lin_w(2) = 0.95     ! gamma_11, gamma_22



  do var = 1, NSYS_VAR
    call get_LegCoeff3(Wl(i-1:i+1,var), coeff_q1L(:))
    call get_LegCoeff3(Wr(i-1:i+1,var), coeff_q1R(:))

    call get_LegCoeff5(Wl(i-2:i+2,var), coeff_q2L(:))
    call get_LegCoeff5(Wr(i-2:i+2,var), coeff_q2R(:))


    q0L = Wl(i, var)
    q0R = Wr(i, var)

    q1L = dot_product(coeff_q1L(:), (/1.,  1./2., 1./6./))    ! left stencil -> +1/2
    q1R = dot_product(coeff_q1R(:), (/1., -1./2., 1./6./))    ! right stencil -> -1/2

    p11L = 1./lin_w(2)*q1L - lin_w(1)/lin_w(2)*q0L
    p11R = 1./lin_w(2)*q1R - lin_w(1)/lin_w(2)*q0R

    smth_indL(1) = MIN( (Wl(i,var) - Wl(i-1,var))**2, (Wl(i+1,var) - Wl(i,var))**2 )    ! beta_01
    smth_indR(1) = MIN( (Wr(i,var) - Wr(i-1,var))**2, (Wr(i+1,var) - Wr(i,var))**2 )

    call Leg_betas(coeff_q1L(:)/lin_w(2), 3, smth_indL(2))      ! beta_11
    call Leg_betas(coeff_q1R(:)/lin_w(2), 3, smth_indR(2))

    call get_nestedNonlinWeights(lin_w(:), smth_indL(:), nonLin_wL(:))
    call get_nestedNonlinWeights(lin_w(:), smth_indR(:), nonLin_wR(:))

    ! third order reconstruction
    p12L = nonLin_wL(1)*q0L + nonLin_wL(2)*p11L
    p12R = nonLin_wR(1)*q0R + nonLin_wR(2)*p11R



    ! start to building fifth order
    q2L = dot_product(coeff_q2L(:), (/ 1.,  1./2., 1./6.,  1./20., 1./70. /))    ! left stencil -> +1/2
    q2R = dot_product(coeff_q2R(:), (/ 1., -1./2., 1./6., -1./20., 1./70. /))    ! right stencil -> -1/2

    p22L = 1./lin_w(2)*q2L - lin_w(1)/lin_w(2)*p12L
    p22R = 1./lin_w(2)*q2R - lin_w(1)/lin_w(2)*p12R

    call Leg_betas(nonLin_wL(2)/lin_w(2)*coeff_q1L(:), 3, smth_indL(1))     ! beta_12
    call Leg_betas(nonLin_wR(2)/lin_w(2)*coeff_q1R(:), 3, smth_indR(1))

    ! building beta_22
    tmpCoeffL(:) = coeff_q2L(:)
    tmpCoeffR(:) = coeff_q2R(:)

    tmpCoeffL(2) = coeff_q2L(2) - lin_w(1)*nonLin_wL(2)/lin_w(2)*coeff_q1L(2)
    tmpCoeffL(3) = coeff_q2L(3) - lin_w(1)*nonLin_wL(2)/lin_w(2)*coeff_q1L(3)
    tmpCoeffL(:) = tmpCoeffL(:)/lin_w(2)

    tmpCoeffR(2) = coeff_q2R(2) - lin_w(1)*nonLin_wR(2)/lin_w(2)*coeff_q1R(2)
    tmpCoeffR(3) = coeff_q2R(3) - lin_w(1)*nonLin_wR(2)/lin_w(2)*coeff_q1R(3)
    tmpCoeffR(:) = tmpCoeffR(:)/lin_w(2)

    call Leg_betas(tmpCoeffL, 5, smth_indL(2))     ! beta_22
    call Leg_betas(tmpCoeffR, 5, smth_indR(2))

    call get_nestedNonlinWeights(lin_w(:), smth_indL(:), nonLin_wL(:))
    call get_nestedNonlinWeights(lin_w(:), smth_indR(:), nonLin_wR(:))

    ! fifth order reconstruction
    reconL(var) = nonLin_wL(1)*p12L + nonLin_wL(2)*p22L
    reconR(var) = nonLin_wR(1)*p12R + nonLin_wR(2)*p22R



  end do



  return
end subroutine soln_nestedWENO5


subroutine get_LegCoeff3(U, coeff)

  implicit none

  real, dimension(3), intent(IN) :: U
  real, dimension(3), intent(OUT) :: coeff

  integer, PARAMETER :: i = 2

  coeff(1) = U(i)
  coeff(2) = 0.5*(U(i+1) - U(i-1))
  coeff(3) = 0.5*(U(i-1) - 2.*U(i) + U(i+1))


end subroutine get_LegCoeff3



subroutine get_LegCoeff5(U, coeff)

  implicit none

  real, dimension(5), intent(IN) :: U
  real, dimension(5), intent(OUT) :: coeff

  integer, PARAMETER :: i = 3

  coeff(1) = U(i)
  coeff(2) = (11.*U(i-2) - 82.*U(i-1)           + 82.*U(i+1) - 11.*U(i+2))/120.
  coeff(3) = (-3.*U(i-2) + 40.*U(i-1) -74.*U(i) + 40.*U(i+1) -  3.*U(i+2))/56.
  coeff(4) = (   -U(i-2) +  2.*U(i-1)           -  2.*U(i+1) +     U(i+2))/12.
  coeff(5) = (    U(i-2) -  4.*U(i-1) + 6.*U(i) -  4.*U(i+1) +     U(i+2))/24.

end subroutine get_LegCoeff5


subroutine get_nestedNonlinWeights(lin_w, smth_ind, nonLin_w)

  use sim_data, only: sim_WENeps

  implicit none

  real, dimension(2), intent(IN) :: lin_w, smth_ind
  real, dimension(2), intent(OUT) :: nonLin_w

  real, dimension(2) :: wbar
  real :: tau, w_norm

  integer :: n

  tau = (smth_ind(2) - smth_ind(1))**2
  do n = 1, 2
    wbar(n) = lin_w(n)*(1. + tau/(sim_WENeps + smth_ind(n)))
  end do

  w_norm = sum(wbar)

  do n = 1, 2
    nonLin_w(n) = wbar(n)/w_norm
  end do

end subroutine get_nestedNonlinWeights
