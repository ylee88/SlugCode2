module char_limiting
  ! this module is unique to FDM

#include "definition.h"

  use eigensystem
  use primconsflux

  implicit none

contains

  subroutine char_proj(R, V, stencil_L, stencil_R, dir)
    ! for FDM, this is a flux splitting.
    ! input: V(var, s)
    ! output; F_left(s, var), F_right(s, var)
    use grid_data, only: gr_maxSpeed

    implicit none

    integer, intent(IN) :: R, dir
    real, dimension(NUMB_VAR, 2*R+1), intent(IN) :: V
    real, dimension(2*R+1, NSYS_VAR), intent(OUT) :: stencil_L, stencil_R

    real, dimension(NUMB_VAR) :: Viph, Vimh
    real, dimension(NSYS_VAR) :: Flux, Cons
    real, dimension(NSYS_VAR, NUMB_WAVE) :: leig_L, leig_R

    logical :: conservative
    integer :: cntr, s, var

    cntr = R+1

    do var = DENS_VAR, NUMB_VAR
      Viph(var) = 0.5*(V(var, cntr) + V(var, cntr+1))
      Vimh(var) = 0.5*(V(var, cntr) + V(var, cntr-1))
    end do

    conservative = .true.
    call left_eigenvectors(Viph(:), conservative, leig_R, dir)
    call left_eigenvectors(Vimh(:), conservative, leig_L, dir)

    do s = 1, 2*R+1
      call prim2flux(V(:,s), Flux(:), dir)
      call prim2cons(V(:,s), Cons(:))
      do var = 1, NUMB_WAVE
        stencil_R(s, var) = 0.5*dot_product( leig_R(:,var), Flux(:) + gr_maxSpeed(var,dir)*Cons(:) )
        stencil_L(s, var) = 0.5*dot_product( leig_L(:,var), Flux(:) - gr_maxSpeed(var,dir)*Cons(:) )
      end do
    end do

    return
  end subroutine char_proj


  subroutine char_back_proj(R, V, charL, charR, fluxL, fluxR, dir)

    implicit none

    integer, intent(IN) :: R, dir
    real, dimension(NUMB_VAR, 2*R+1), intent(IN) :: V
    real, dimension(NSYS_VAR), intent(IN) :: charL, charR
    real, dimension(NSYS_VAR), intent(OUT) :: fluxL, fluxR

    real, dimension(NUMB_VAR) :: Viph, Vimh
    real, dimension(NSYS_VAR, NUMB_WAVE) :: reig_L, reig_R

    logical :: conservative
    integer :: cntr, var

    cntr = R+1

    do var = DENS_VAR, NUMB_VAR
      Viph(var) = 0.5*(V(var, cntr) + V(var, cntr+1))
      Vimh(var) = 0.5*(V(var, cntr) + V(var, cntr-1))
    end do

    conservative = .true.
    call right_eigenvectors(Viph(:), conservative, reig_R, dir)
    call right_eigenvectors(Vimh(:), conservative, reig_L, dir)

    do var = 1, NUMB_WAVE
      fluxR(var) = dot_product( charR(:), reig_R(var,:) )
      fluxL(var) = dot_product( charL(:), reig_L(var,:) )
    end do

    return
  end subroutine char_back_proj



end module char_limiting
