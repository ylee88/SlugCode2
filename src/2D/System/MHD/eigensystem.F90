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
    !! @note Reference: K.G. Powell (1994) <https://ntrs.nasa.gov/search.jsp?R=19940028527>
    !! @version beta
    !! @date 2020-07-10
    subroutine eigenvalues(V, lambda, dir)
        implicit none
        real, dimension(NUMB_VAR), intent(IN) :: V
        real, dimension(NUMB_WAVE), intent(OUT) :: lambda
        integer, intent(IN) :: dir

        real :: a ! sound speed
        real :: u ! entropy wave speed
        real :: b ! normal magnetic field
        real :: B2
        real :: ca
        real :: cfs_aux
        real :: cf
        real :: cs
        integer :: VEL ! velocity index
        integer :: MAG ! magnetic field index

        if (dir == XDIM) then
            VEL = VELX_VAR
            MAG = MAGX_VAR
        else if(dir == YDIM) then
            VEL = VELY_VAR
            MAG = MAGY_VAR
        else
        end if

        a = sqrt(V(GAMC_VAR) * V(PRES_VAR) / V(DENS_VAR))
        u = V(VEL)
        b = V(MAG)
        B2 = V(MAGX_VAR) ** 2.0 + V(MAGY_VAR) ** 2.0 + V(MAGZ_VAR) ** 2.0
        ca = b / sqrt(V(DENS_VAR))
        cfs_aux = sqrt(((V(GAMC_VAR) * V(PRES_VAR) + B2) / V(DENS_VAR)) ** 2.0 - (2.0 * a * ca) ** 2.0)
        cf = 0.5 * (a ** 2.0 + B2 / V(DENS_VAR) + cfs_aux)
        cs = 0.5 * (a ** 2.0 + B2 / V(DENS_VAR) - cfs_aux)

        lambda(FAST_LEFT) = u - cf
        lambda(ALFVEN_LEFT) = u - ca
        lambda(SLOW_LEFT) = u - cs
        lambda(ENTROPY) = u
        lambda(SLOW_RIGHT) = u + cs
        lambda(ALFVEN_RIGHT) = u + ca
        lambda(FAST_RIGHT) = u + cf
        lambda(DIVERGENCE) = u

        return
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
