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

    if (dir == XDIM) then
      VEL = VELX_VAR
    else
      VEL = VELY_VAR
    end if

    ! sound speed
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    u = V(VEL)

    lambda(SHOCKLEFT) = u - a
    lambda(SLOWWLEFT) = u
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

    real :: a, d, ekin, hdai, H, p, v1, v2, g
    integer :: VEL1_VAR, VEL2_VAR

    VEL1_VAR = 1; VEL2_VAR = 2
    if (dir == XDIM) then
      VEL1_VAR = VELX_VAR
      VEL2_VAR = VELY_VAR
    else if (dir == YDIM) then
      VEL1_VAR = VELY_VAR
      VEL2_VAR = VELX_VAR
    end if

    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    g = V(GAMC_VAR) - 1.
    v1 = V(VEL1_VAR)
    v2 = V(VEL2_VAR)
    d = V(DENS_VAR)
    p = V(PRES_VAR)
    ekin = 0.5*(v1**2+v2**2)
    H = ekin + a**2/g
    hdai = 0.5/(a*d)

    if (conservative) then
      !! Conservative eigenvector
      reig(DENS_VAR,SHOCKLEFT) = 1.
      reig(VEL1_VAR,SHOCKLEFT) = v1 - a
      reig(VEL2_VAR,SHOCKLEFT) = v2
      reig(PRES_VAR,SHOCKLEFT) = H - a*v1

      reig(DENS_VAR,SLOWWLEFT) = 0.
      reig(VEL1_VAR,SLOWWLEFT) = 0.
      reig(VEL2_VAR,SLOWWLEFT) = 1.
      reig(PRES_VAR,SLOWWLEFT) = v2

      reig(DENS_VAR,CTENTROPY) = 1.
      reig(VEL1_VAR,CTENTROPY) = v1
      reig(VEL2_VAR,CTENTROPY) = v2
      reig(PRES_VAR,CTENTROPY) = ekin

      reig(DENS_VAR,SHOCKRGHT) = 1.
      reig(VEL1_VAR,SHOCKRGHT) = v1 + a
      reig(VEL2_VAR,SHOCKRGHT) = v2
      reig(PRES_VAR,SHOCKRGHT) = H + v1*a

    else
      !! Primitive eigenvector
      reig(DENS_VAR,SHOCKLEFT) = 1.
      reig(VEL1_VAR,SHOCKLEFT) = -hdai
      reig(VEL2_VAR,SHOCKLEFT) = 0.
      reig(PRES_VAR,SHOCKLEFT) = 0.5

      reig(DENS_VAR,SLOWWLEFT) = 0.
      reig(VEL1_VAR,SLOWWLEFT) = 0.
      reig(VEL2_VAR,SLOWWLEFT) = 1
      reig(PRES_VAR,SLOWWLEFT) = 0.

      reig(DENS_VAR,CTENTROPY) = 1.
      reig(VEL1_VAR,CTENTROPY) = 0.
      reig(VEL2_VAR,CTENTROPY) = 0.
      reig(PRES_VAR,CTENTROPY) = 0.

      reig(DENS_VAR,SHOCKRGHT) = 1.
      reig(VEL1_VAR,SHOCKRGHT) = hdai
      reig(VEL2_VAR,SHOCKRGHT) = 0.
      reig(PRES_VAR,SHOCKRGHT) = 0.5

    endif

    return
  end subroutine right_eigenvectors


  subroutine left_eigenvectors(V,conservative,leig,dir)
    implicit none
    integer, intent(IN) :: dir
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: leig

    real :: a, d, g, ekin, ha2i, agi,  p, v1, v2, ad
    integer :: VEL1_VAR, VEL2_VAR
    VEL1_VAR = 1; VEL2_VAR = 2
    if (dir == XDIM) then
      VEL1_VAR = VELX_VAR
      VEL2_VAR = VELY_VAR
    else if (dir == YDIM) then
      VEL1_VAR = VELY_VAR
      VEL2_VAR = VELX_VAR
    end if

    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    v1 = V(VEL1_VAR)
    v2 = V(VEL2_VAR)
    d = V(DENS_VAR)
    p = V(PRES_VAR)
    g = V(GAMC_VAR) - 1.
    ad = a*d
    ha2i = 0.5/(a*a)
    agi = a/g
    ekin = 0.5*(v1**2+v2**2)

    if (conservative) then
      !!$       !! Conservative eigenvector
      leig(DENS_VAR,SHOCKLEFT) = ekin + agi*v1
      leig(VEL1_VAR,SHOCKLEFT) = -agi - v1
      leig(VEL2_VAR,SHOCKLEFT) = -v2
      leig(PRES_VAR,SHOCKLEFT) = 1.

      leig(DENS_VAR,SLOWWLEFT) = -2.*v2*a*agi
      leig(VEL1_VAR,SLOWWLEFT) = 0.
      leig(VEL2_VAR,SLOWWLEFT) = 2*a*agi
      leig(PRES_VAR,SLOWWLEFT) = 0.

      leig(DENS_VAR,CTENTROPY) = 2*a*agi - 2.*ekin
      leig(VEL1_VAR,CTENTROPY) = 2.*v1
      leig(VEL2_VAR,CTENTROPY) = 2.*v2
      leig(PRES_VAR,CTENTROPY) = -2.

      leig(DENS_VAR,SHOCKRGHT) = ekin - v1*agi
      leig(VEL1_VAR,SHOCKRGHT) = -v1 + agi
      leig(VEL2_VAR,SHOCKRGHT) = -v2
      leig(PRES_VAR,SHOCKRGHT) = 1.

      leig(:,:) = leig(:,:)*g/(2.*a*a)
    else
      !! Primitive eigenvector
      leig(DENS_VAR,SHOCKLEFT) = 0.
      leig(VEL1_VAR,SHOCKLEFT) = -ad
      leig(VEL2_VAR,SHOCKLEFT) = 0.
      leig(PRES_VAR,SHOCKLEFT) = 1.

      leig(DENS_VAR,SLOWWLEFT) = 0.
      leig(VEL1_VAR,SLOWWLEFT) = 0.
      leig(VEL2_VAR,SLOWWLEFT) = 1.
      leig(PRES_VAR,SLOWWLEFT) = 0.

      leig(DENS_VAR,CTENTROPY) = 1.
      leig(VEL1_VAR,CTENTROPY) = 0.
      leig(VEL2_VAR,CTENTROPY) = 0.
      leig(PRES_VAR,CTENTROPY) = -ha2i

      leig(DENS_VAR,SHOCKRGHT) = 0.
      leig(VEL1_VAR,SHOCKRGHT) = ad
      leig(VEL2_VAR,SHOCKRGHT) = 0.
      leig(PRES_VAR,SHOCKRGHT) = 1.

    endif

    return
  end subroutine left_eigenvectors



end module eigensystem
