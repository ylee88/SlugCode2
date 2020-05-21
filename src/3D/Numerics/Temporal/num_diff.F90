module num_diff

#include "definition.h"

  use grid_data, only: gr_dx, gr_dy, gr_dz

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



  function diffxy(q, Nx, d1, d2) result(f)
    ! it takes Nx x Nx 2D stencil
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx, Nx), intent(IN) :: q
    real, intent(IN) :: d1, d2
    real, dimension(NSYS_VAR) :: f

    real :: dxdy
    integer :: i, j

    dxdy = d1*d2

    if(Nx == 3) then
      i = 2
      j = 2

      f = q(:,i+1,j+1) - q(:,i-1,j+1) - q(:,i+1,j-1) + q(:,i-1,j-1)
      f = f/(4.*dxdy)

    else if(Nx == 5) then
      i = 3
      j = 3

      f =     q(:,i-2,j-2) - 8.*q(:,i-2,j-1) + 8.*q(:,i-2,j+1) - q(:,i-2,j+2)   &
        -8.*( q(:,i-1,j-2) - 8.*q(:,i-1,j-1) + 8.*q(:,i-1,j+1) - q(:,i-1,j+2) ) &
        +8.*( q(:,i+1,j-2) - 8.*q(:,i+1,j-1) + 8.*q(:,i+1,j+1) - q(:,i+1,j+2) ) &
           -( q(:,i+2,j-2) - 8.*q(:,i+2,j-1) + 8.*q(:,i+2,j+1) - q(:,i+2,j+2) )

      f = f/(144.*dxdy)

    else if (Nx == 7) then
      i = 4
      j = 4

      f =  -1.*( -q(:,i-3,j-3) + 9.*q(:,i-3,j-2) - 45.*q(:,i-3,j-1) + 45.*q(:,i-3,j+1) - 9.*q(:,i-3,j+2) + q(:,i-3,j+3) ) &
           +9.*( -q(:,i-2,j-3) + 9.*q(:,i-2,j-2) - 45.*q(:,i-2,j-1) + 45.*q(:,i-2,j+1) - 9.*q(:,i-2,j+2) + q(:,i-2,j+3) ) &
           -45.*( -q(:,i-1,j-3) + 9.*q(:,i-1,j-2) - 45.*q(:,i-1,j-1) + 45.*q(:,i-1,j+1) - 9.*q(:,i-1,j+2) + q(:,i-1,j+3) ) &
           +45.*( -q(:,i+1,j-3) + 9.*q(:,i+1,j-2) - 45.*q(:,i+1,j-1) + 45.*q(:,i+1,j+1) - 9.*q(:,i+1,j+2) + q(:,i+1,j+3) ) &
           -9.*( -q(:,i+2,j-3) + 9.*q(:,i+2,j-2) - 45.*q(:,i+2,j-1) + 45.*q(:,i+2,j+1) - 9.*q(:,i+2,j+2) + q(:,i+2,j+3) ) &
           +1.*( -q(:,i+3,j-3) + 9.*q(:,i+3,j-2) - 45.*q(:,i+3,j-1) + 45.*q(:,i+3,j+1) - 9.*q(:,i+3,j+2) + q(:,i+3,j+3) )

      f = f/(3600.*dxdy)

    end if

  end function diffxy

  function diffxxy(q, Nx, d1, d2) result(f)
    ! it takes Nx x Nx 2D stencil
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx, Nx), intent(IN) :: q
    real, intent(IN) :: d1, d2
    real, dimension(NSYS_VAR) :: f

    real :: dx2dy
    integer :: i, j

    dx2dy = d1**2*d2

    if(Nx == 3) then
      i = 2
      j = 2

      f =  -(q(:,i-1,j-1) + q(:,i+1,j-1)) &
           +(q(:,i-1,j+1) + q(:,i+1,j+1))
      f = f/(2*dx2dy)

    else if(Nx == 5) then
      i = 3
      j = 3

      f =     ( -q(:,i-2,j-2) + 16.*q(:,i-1,j-2) - 30.*q(:,i,j-2) + 16.*q(:,i+1,j-2) - q(:,i+2,j-2) ) &
          -8.*( -q(:,i-2,j-1) + 16.*q(:,i-1,j-1) - 30.*q(:,i,j-1) + 16.*q(:,i+1,j-1) - q(:,i+2,j-1) ) &
          +8.*( -q(:,i-2,j+1) + 16.*q(:,i-1,j+1) - 30.*q(:,i,j+1) + 16.*q(:,i+1,j+1) - q(:,i+2,j+1) ) &
             -( -q(:,i-2,j+2) + 16.*q(:,i-1,j+2) - 30.*q(:,i,j+2) + 16.*q(:,i+1,j+2) - q(:,i+2,j+2) )

      f = f/(144.*dx2dy)

    end if

  end function diffxxy


  function diffxyy(q, Nx, d1, d2) result(f)
    ! it takes Nx x Nx 2D stencil
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx, Nx), intent(IN) :: q
    real, intent(IN) :: d1, d2
    real, dimension(NSYS_VAR) :: f

    real :: dxdy2
    integer :: i, j

    dxdy2 = d1*d2**2

    if(Nx == 3) then
      i = 2
      j = 2

      f =  (-q(:,i-1,j-1) + q(:,i-1,j+1)) &
          +(-q(:,i+1,j-1) + q(:,i+1,j+1))
      f = f/(2*dxdy2)

    else if(Nx == 5) then
      i = 3
      j = 3

      f =     ( -q(:,i-2,j-2) + 16.*q(:,i-2,j-1) - 30.*q(:,i-2,j) + 16.*q(:,i-2,j+1) - q(:,i-2,j+2) ) &
          -8.*( -q(:,i-1,j-2) + 16.*q(:,i-1,j-1) - 30.*q(:,i-1,j) + 16.*q(:,i-1,j+1) - q(:,i-1,j+2) ) &
          +8.*( -q(:,i+1,j-2) + 16.*q(:,i+1,j-1) - 30.*q(:,i+1,j) + 16.*q(:,i+1,j+1) - q(:,i+1,j+2) ) &
             -( -q(:,i+2,j-2) + 16.*q(:,i+2,j-1) - 30.*q(:,i+2,j) + 16.*q(:,i+2,j+1) - q(:,i+2,j+2) )

      f = f/(144.*dxdy2)

    end if

  end function diffxyy


  function diffxyz(q, Nx) result(f)
    ! it takes Nx x Nx 2D stencil
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx, Nx, Nx), intent(IN) :: q
    real, dimension(NSYS_VAR) :: f1, f2, f3, f4, f

    integer :: i, j, k

    if(Nx == 5) then
      i = 3
      j = 3
      k = 3

      f1 = diffxy(q(:, i-2:i+2, j-2:j+2, k-2), 5, gr_dx, gr_dy)
      f2 = diffxy(q(:, i-2:i+2, j-2:j+2, k-1), 5, gr_dx, gr_dy)
      f3 = diffxy(q(:, i-2:i+2, j-2:j+2, k+1), 5, gr_dx, gr_dy)
      f4 = diffxy(q(:, i-2:i+2, j-2:j+2, k+2), 5, gr_dx, gr_dy)

      f = f1 - 8.*f2 + 8.*f3 - f4
      f = f/12.*gr_dz
    else
      call abort_slug("[diffxyz] unrecognized Nx")
    end if
  end function diffxyz


end module num_diff
