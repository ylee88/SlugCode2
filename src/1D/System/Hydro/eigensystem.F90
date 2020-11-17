module eigensystem

#include "definition.h"

  use grid_data

contains

  subroutine eigenvalues(V,lambda,dir)
    implicit none

    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NUMB_WAVE),intent(OUT) :: lambda
    integer, intent(IN) :: dir

    real :: a, u
    integer :: VEL

    select case(dir)
    case(XDIM)
      VEL = VELX_VAR
    case default
      VEL = 100
      call abort_slug("[eigval] unrecognized dir")
    end select

    ! sound speed
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VEL)

    lambda(SHOCKLEFT) = u - a
    lambda(CTENTROPY) = u
    lambda(SHOCKRGHT) = u + a

    return
  end subroutine eigenvalues


  subroutine right_eigenvectors(V,conservative,reig, dir)
    implicit none
    integer, intent(IN) :: dir
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: reig

    real :: d, P, a, hda, hdai, v1,  H, ekin, g
    integer :: VEL1_VAR

    VEL1_VAR = 1
    if (dir == XDIM) then
      VEL1_VAR = VELX_VAR
    else
      call abort_slug("[reigvec] invalid dir. stop")
    end if

    ! sound speed, and others
    d = V(DENS_VAR)
    P = V(PRES_VAR)
    g = V(GAMC_VAR) - 1.
    a = SQRT(V(GAMC_VAR)*P/d)
    hdai = 0.5*d/a
    hda = 0.5*d*a
    v1 = V(VEL1_VAR)
    ekin = 0.5*(v1*v1)
    H = ekin + a**2/g

    if (conservative) then
      !! Conservative eigenvector
      reig(DENS_VAR,SHOCKLEFT) = 1.
      reig(VEL1_VAR,SHOCKLEFT) = v1 - a
      reig(PRES_VAR,SHOCKLEFT) = H - a*v1

      reig(DENS_VAR,CTENTROPY) = 1.
      reig(VEL1_VAR,CTENTROPY) = v1
      reig(PRES_VAR,CTENTROPY) = ekin

      reig(DENS_VAR,SHOCKRGHT) = 1.
      reig(VEL1_VAR,SHOCKRGHT) = v1 + a
      reig(PRES_VAR,SHOCKRGHT) = H + v1*a
    else
      !! Primitive eigenvector
      reig(DENS_VAR,SHOCKLEFT) = -hdai
      reig(VEL1_VAR,SHOCKLEFT) = 0.5
      reig(PRES_VAR,SHOCKLEFT) = -hda

      reig(DENS_VAR,CTENTROPY) = 1.
      reig(VEL1_VAR,CTENTROPY) = 0.
      reig(PRES_VAR,CTENTROPY) = 0.

      reig(DENS_VAR,SHOCKRGHT) = hdai
      reig(VEL1_VAR,SHOCKRGHT) = 0.5
      reig(PRES_VAR,SHOCKRGHT) = hda
    endif

    return
  end subroutine right_eigenvectors


  subroutine left_eigenvectors(V,conservative,leig,dir)
    implicit none
    integer, intent(IN) :: dir
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: leig

    real :: d, P, a, Na, a2inv, v1, ekin, g
    integer :: VEL1_VAR
    VEL1_VAR = 1
    if (dir == XDIM) then
      VEL1_VAR = VELX_VAR
    else
      call abort_slug("[reigvec] invalid dir. stop")
    end if

    ! sound speed, and others
    d = V(DENS_VAR)
    P = V(PRES_VAR)
    g = V(GAME_VAR) - 1.
    a = SQRT(V(GAMC_VAR)*P/d)
    Na = 0.5/(a*a)
    a2inv = 1./(a*a)
    v1 = V(VEL1_VAR)
    ekin = 0.5*(v1*v1)

    if (conservative) then
      !! Conservative eigenvector
      leig(DENS_VAR,SHOCKLEFT) = g*ekin+v1*a
      leig(VEL1_VAR,SHOCKLEFT) = -g*v1-a
      leig(PRES_VAR,SHOCKLEFT) = g

      leig(DENS_VAR:PRES_VAR,SHOCKLEFT) = Na*leig(DENS_VAR:PRES_VAR,SHOCKLEFT)

      leig(DENS_VAR,CTENTROPY) = 1.-Na*g*2.*ekin
      leig(VEL1_VAR,CTENTROPY) = g*v1*a2inv
      leig(PRES_VAR,CTENTROPY) = -g*a2inv

      leig(DENS_VAR,SHOCKRGHT) = g*ekin-v1*a
      leig(VEL1_VAR,SHOCKRGHT) = -g*v1+a
      leig(PRES_VAR,SHOCKRGHT) = g

      leig(DENS_VAR:PRES_VAR,SHOCKRGHT) = Na*leig(DENS_VAR:PRES_VAR,SHOCKRGHT)
    else
      !! Primitive eigenvector
      leig(DENS_VAR,SHOCKLEFT) = 0.0
      leig(VEL1_VAR,SHOCKLEFT) = 1.0
      leig(PRES_VAR,SHOCKLEFT) = -1./(d*a)

      leig(DENS_VAR,CTENTROPY) = 1.
      leig(VEL1_VAR,CTENTROPY) = 0.
      leig(PRES_VAR,CTENTROPY) = -a2inv

      leig(DENS_VAR,SHOCKRGHT) = 0.0
      leig(VEL1_VAR,SHOCKRGHT) = 1.0
      leig(PRES_VAR,SHOCKRGHT) = 1./(d*a)
    endif

    return
  end subroutine left_eigenvectors



end module eigensystem
