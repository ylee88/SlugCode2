!> @file 2D/System/MHD/eigensystem.F90
!! @brief This file contains methods to solve for the eigensystem of 2D MHD
!! problems.

!> @brief Methods to solve to the eigensystem of 2D MHD problems.
!! @version beta
!! @warning This module is still in development.
module eigensystem

#include "definition.h"

    use grid_data

contains

    !> Calculate eigenvalues based on primitive variables.
    !! @param V Primitive variables.
    !! @param lambda Eigenvalues.
    !! @param dir Direction.
    !! @version beta
    !! @date 2020-07-09
    subroutine eigenvalues(V, lambda, dir)
        implicit none
        real, dimension(NUMB_VAR), intent(IN) :: V
        real, dimension(NUMB_WAVE), intent(OUT) :: lambda
        integer, intent(IN) :: dir
    end subroutine eigenvalues

    !> Calculate right eigenvectors based on primitive variables.
    !! @param V Primitive variables.
    !! @param conservative Flag for conservative eigenvectors.
    !! @param reig Right eigenvectors.
    !! @param dir Direction.
    !! @version beta
    !! @date 2020-07-09
    subroutine right_eigenvectors(V, conservative, reig, dir)
        implicit none
        real, dimension(NUMB_VAR), intent(IN) :: V
        logical :: conservative
        real, dimension(NSYS_VAR, NUMB_WAVE), intent(OUT) :: reig
        integer, intent(IN) :: dir
    end subroutine right_eigenvectors

    !> Calculate left eigenvectors based on primitive variables.
    !! @param V Primitive variables.
    !! @param conservative Flag for conservative eigenvectors.
    !! @param leig Left eigenvectors.
    !! @param dir Direction.
    !! @version beta
    !! @date 2020-07-09
    subroutine left_eigenvectors(V, conservative, leig, dir)
        implicit none
        real, dimension(NUMB_VAR), intent(IN) :: V
        logical :: conservative
        real, dimension(NSYS_VAR, NUMB_WAVE), intent(OUT) :: leig
        integer, intent(IN) :: dir
    end subroutine left_eigenvectors

end module eigensystem
