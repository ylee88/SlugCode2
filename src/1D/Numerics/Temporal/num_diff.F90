module num_diff

#include "definition.h"

  use grid_data, only: gr_dx

contains


  !!!!! experimental
  function diff1_weno(q, Nx, d) result(f)

    use WENO
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx), intent(IN) :: q
    real, intent(IN) :: d

    real, dimension(NSYS_VAR) :: f
    real, dimension(3) :: smth_ind, linW, nonLinW
    real, dimension(5) :: C, cl, cm, cr

    real :: delta, sumW
    integer :: var, k

    delta = d

    C  = (/  1., -8.,  0.,  8., -1. /)/(12.*delta)
    cl = (/  1., -4.,  3.,  0.,  0. /)/(2.*delta)
    cm = (/  0., -1.,  0.,  1.,  0. /)/(2.*delta)
    cr = (/  0.,  0., -3.,  4., -1. /)/(2.*delta)

    linW = (/ 1., 4., 1. /)/6.

    do var = 1, NSYS_VAR
      call betas(q(var, :), 2, smth_ind)
    end do

    do k = 1, 3
      nonLinW(k) = linW(k)/(1.E-36 + smth_ind(k))**1
    end do

    sumW = SUM(nonLinW)
    nonLinW(:) = nonLinW(:)/sumW

    C(:) = nonLinW(1)*cl(:) + nonLinW(2)*cm + nonLinW(3)*cr

    do var = 1, NSYS_VAR
      f(var) = dot_product(q(var, :), C)
    end do

  end function diff1_weno
  !!!!!


  function diff1(q, Nx, d) result(f)
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx), intent(IN) :: q
    real, intent(IN) :: d

    real, dimension(NSYS_VAR) :: f

    real :: delta
    integer :: i

    delta = d

    if(Nx == 3) then
      i = 2

      f = -q(:,i-1) + q(:,i+1)
      f = f/(2.*delta)

    else if(Nx == 5) then
      i = 3

      f = q(:,i-2) - 8.*q(:,i-1) + 8.*q(:,i+1) - q(:,i+2)
      f = f/(12.*delta)

    else if(Nx == 7) then
      i = 4

      f = -q(:,i-3) + 9.*q(:,i-2) - 45.*q(:,i-1) + 45.*q(:,i+1) - 9.*q(:,i+2) + q(:,i+3)
      f = f/(60.*delta)

    end if

    return
  end function diff1

  function diff2(q, Nx, d) result(f)
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx), intent(IN) :: q
    real, intent(IN) :: d

    real, dimension(NSYS_VAR) :: f

    real :: delta
    integer :: i

    delta = d

    if(Nx == 3) then
      i = 2

      f = q(:,i-1) - 2.*q(:,i) + q(:,i+1)
      f = f/delta**2

    else if(Nx == 5) then
      i = 3

      f = -q(:,i-2) + 16.*q(:,i-1) - 30.*q(:,i) + 16.*q(:,i+1) - q(:,i+2)
      f = f/(12.*delta**2)

    else if(Nx == 7) then

      i = 4
      f = 2.*q(:,i-3) - 27.*q(:,i-2) + 270.*q(:,i-1) - 490.*q(:,i) + 270.*q(:,i+1) - 27.*q(:,i+2) + 2.*q(:,i+3)
      f = f/(180.*delta**2)

    end if

    return
  end function diff2

  function diff3(q, Nx, d) result(f)
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx), intent(IN) :: q
    real, intent(IN) :: d
    real, dimension(NSYS_VAR) :: f

    real :: delta
    integer :: i

    delta = d

    if(Nx == 5) then
      i = 3

      f = -q(:,i-2) + 2.*q(:,i-1) -2.*q(:,i+1) + q(:,i+2)
      f = f/(2.*delta**3)
    else if(Nx == 7) then
      i = 4

      f = q(:,i-3) - 8.*q(:,i-2) + 13.*q(:,i-1) - 13.*q(:,i+1) + 8.*q(:,i+2) - q(:,i+3)
      f = f/(8.*delta**3)
    end if

  end function diff3



end module num_diff
