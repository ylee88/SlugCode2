module num_diff

#include "definition.h"

  use grid_data, only: gr_dx, gr_dy

contains


  function diff1(q, Nx, dir) result(f)
    implicit none

    integer, intent(IN) :: Nx, dir
    real, dimension(NSYS_VAR, Nx), intent(IN) :: q

    real, dimension(NSYS_VAR) :: f

    real :: delta
    integer :: i

    select case(dir)
    case(XDIM)
      delta = gr_dx
    case(YDIM)
      delta = gr_dy
    case DEFAULT
      delta = 0.
      call abort_slug("[diff1] Wrong dir value")
    end select

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

  function diff2(q, Nx, dir) result(f)
    implicit none

    integer, intent(IN) :: Nx, dir
    real, dimension(NSYS_VAR, Nx), intent(IN) :: q

    real, dimension(NSYS_VAR) :: f

    real :: delta
    integer :: i

    select case(dir)
    case(XDIM)
      delta = gr_dx
    case(YDIM)
      delta = gr_dy
    case DEFAULT
      delta = 0
      call abort_slug("[diff2] Wrong dir value")
    end select

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

  function diff3(q, Nx, dir) result(f)
    implicit none

    integer, intent(IN) :: Nx, dir
    real, dimension(NSYS_VAR, Nx), intent(IN) :: q
    real, dimension(NSYS_VAR) :: f

    real :: delta
    integer :: i

    select case(dir)
    case(XDIM)
      delta = gr_dx
    case(YDIM)
      delta = gr_dy
    case DEFAULT
      delta = 0.
      call abort_slug("[diff3] Wrong dir value")
    end select

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



  function diffxy(q, Nx) result(f)
    ! it takes Nx x Nx 2D stencil
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx, Nx), intent(IN) :: q
    real, dimension(NSYS_VAR) :: f

    real :: dxdy
    integer :: i, j

    dxdy = gr_dx*gr_dy

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

  function diffxxy(q, Nx) result(f)
    ! it takes Nx x Nx 2D stencil
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx, Nx), intent(IN) :: q
    real, dimension(NSYS_VAR) :: f

    real :: dx1, dx2, dx2dy
    integer :: i, j

    dx2dy = gr_dx**2*gr_dy

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


  function diffxyy(q, Nx) result(f)
    ! it takes Nx x Nx 2D stencil
    implicit none

    integer, intent(IN) :: Nx
    real, dimension(NSYS_VAR, Nx, Nx), intent(IN) :: q
    real, dimension(NSYS_VAR) :: f

    real :: dxdy2
    integer :: i, j

    dxdy2 = gr_dx*gr_dy**2

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


end module num_diff
