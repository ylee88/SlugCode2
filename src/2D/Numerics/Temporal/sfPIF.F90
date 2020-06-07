module sfPIF

#include "definition.h"

  use primconsflux, only: cons2flux
  use sim_data, only: sim_Torder  ! this will be needed for getting eps

  implicit none

contains



  function get_Jv(U, V, dt, dir) result(r)
    implicit none

    real, dimension(NSYS_VAR), intent(IN) :: U, V

    integer, intent(IN) :: dir
    real, intent(IN) :: dt
    real, dimension(NSYS_VAR) :: r

    real, dimension(NSYS_VAR) :: Fr, Fl
    real :: b_jac, eps

    b_jac = 4.8062E-06
    eps =  get_eps(dt, V, b_jac)

    call cons2flux(U+eps*V, Fr, dir)
    call cons2flux(U-eps*V, Fl, dir)

    r = (Fr - Fl)/(2.*eps)

  end function get_Jv

  function get_Hvw(U, V, W, dt, dir) result(r)
    implicit none
    real, dimension(NSYS_VAR), intent(IN) :: U, V, W
    real, intent(IN) :: dt
    integer, intent(IN) :: dir
    real, dimension(NSYS_VAR) :: r

    real, dimension(NSYS_VAR) :: Fm2, Fm1, Fp1, Fp2, VpW, VmW

    real :: b, eps1, eps2

    b = 1.6065E-4
    eps1 =  get_eps(dt, V, b)
    eps2 =  get_eps(dt, W, b)

    VpW = eps1*V+eps2*W
    VmW = eps1*V-eps2*W

    call cons2flux(U + VpW, Fm2, dir)
    call cons2flux(U + VmW, Fm1, dir)
    call cons2flux(U - VmW, Fp1, dir)
    call cons2flux(U - VpW, Fp2, dir)

    r = (Fm2 - Fm1 - Fp1 + Fp2)/(4.*eps1*eps2)

  end function get_Hvw

  function get_Dvwx(U, V, W, X, dt, dir) result(r)
    implicit none
    real, dimension(NSYS_VAR), intent(IN) :: U, V, W, X
    real, intent(IN) :: dt
    integer, intent(IN) :: dir
    real, dimension(NSYS_VAR) :: r

    real, dimension(NSYS_VAR) :: Fm4, Fm3, Fm2, Fm1, Fp1, Fp2, Fp3, Fp4
    real, dimension(NSYS_VAR) :: C1, C2, C3, C4

    real :: b, eps1, eps2, eps3

    b = 0.0005673365502470651   ! see note
    eps1 = get_eps(dt, V, b)
    eps2 = get_eps(dt, W, b)
    eps3 = get_eps(dt, X, b)

    C1 = eps1*V+eps2*W+eps3*X
    C2 = eps1*V-eps2*W+eps3*X
    C3 = eps1*V-eps2*W-eps3*X
    C4 = eps1*V+eps2*W-eps3*X

    call cons2flux(U+C1, Fm4, dir)
    call cons2flux(U+C2, Fm3, dir)
    call cons2flux(U-C3, Fm2, dir)
    call cons2flux(U-C4, Fm1, dir)

    call cons2flux(U+C4, Fp1, dir)
    call cons2flux(U+C3, Fp2, dir)
    call cons2flux(U-C2, Fp3, dir)
    call cons2flux(U-C1, Fp4, dir)

    r = (Fm4 - Fm3 - Fm2 + Fm1 - Fp1 + Fp2 + Fp3 - Fp4)/(8.*eps1*eps2*eps3)

  end function get_Dvwx


  function get_eps(dt, V, b) result(r)
    implicit none

    real, dimension(NSYS_VAR), intent(IN) :: V
    real, intent(IN) :: dt, b

    real :: r
    real :: eps, tmp

    eps = sqrt(b)/NORM2(V)
    tmp = min(dt**(sim_Torder-1), eps**2)
    r = sqrt(tmp)

  end function get_eps


end module sfPIF
