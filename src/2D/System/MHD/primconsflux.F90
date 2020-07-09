!> @file 2D/System/MHD/primconsflux.F90
!! @brief This file contains methods to calculate primitive variables,
!! conservative varibales, and numerical flux.

!> @brief Methods to calculate primitive variables, conservative variables,
!! and numerical flux.
!! @version beta
!! @warning This module is still in development.
module primconsflux

#include "definition.h"

    use grid_data
    use sim_data, only : sim_gamma, sim_smallPres

contains

    !> Calculate conservative variables based on primitive variables.
    !! @param V Primitive variables.
    !! @param U Conservative variables.
    !! @version beta
    !! @date 2020-07-09
    subroutine prim2cons(V, U)
        implicit none
        real, dimension(NUMB_VAR), intent(IN) :: V
        real, dimension(NSYS_VAR), intent(OUT) :: U
    end subroutine prim2cons

    !> Calculate primitive variables based on conservative variables.
    !! @param U Conservative variables.
    !! @param V Primitive variables.
    !! @version beta
    !! @date 2020-07-09
    subroutine cons2prim(U, V)
        implicit none
        real, dimension(NSYS_VAR), intent(IN) :: U
        real, dimension(NUMB_VAR), intent(OUT) :: V
    end subroutine cons2prim

    !> Calculate numerical flux based on primitive variables.
    !! @param V Primitive variables.
    !! @param flux Numerical flux.
    !! @param dir Direction.
    !! @version beta
    !! @date 2020-07-09
    subroutine prim2flux(V, flux, dir)
        implicit none
        real, dimension(NUMB_VAR), intent(IN) :: V
        real, dimension(NSYS_VAR), intent(OUT) :: flux
        integer, intent(IN) :: dir
    end subroutine prim2flux

    !> Calculate numerical flux based on conservative variables.
    !! @param U Conservative variables.
    !! @param flux Numerical flux.
    !! @param dir Direction.
    !! @version beta
    !! @date 2020-07-09
    subroutine cons2flux(U, flux, dir)
        implicit none
        real, dimension(NSYS_VAR), intent(IN) :: U
        real, dimension(NSYS_VAR), intent(OUT) :: flux
        integer, intent(IN) :: dir
    end subroutine cons2flux

end module primconsflux
