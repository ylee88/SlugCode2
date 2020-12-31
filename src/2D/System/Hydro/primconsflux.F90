module primconsflux

#include "definition.h"

  use grid_data
  use sim_data, only : sim_gamma, sim_smallPres
  ! use eos, only : eos_cell

contains

  subroutine prim2cons(V,U)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: U

    real :: ekin, eint

    U(DENS_VAR) = V(DENS_VAR)
    U(MOMX_VAR) = V(DENS_VAR)*V(VELX_VAR)
    U(MOMY_VAR) = V(DENS_VAR)*V(VELY_VAR)
    ekin = 0.5*V(DENS_VAR)*(V(VELX_VAR)**2+V(VELY_VAR)**2)
    ! eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    eint = V(EINT_VAR)
    U(ENER_VAR) = ekin + eint

  end subroutine prim2cons


  subroutine cons2prim(U,V)
    implicit none
    real, dimension(NSYS_VAR), intent(IN)  :: U
    real, dimension(NUMB_VAR), intent(OUT) :: V
    real :: eint, ekin, pres

    V(DENS_VAR) = U(DENS_VAR)
    V(VELX_VAR) = U(MOMX_VAR)/U(DENS_VAR)
    V(VELY_VAR) = U(MOMY_VAR)/U(DENS_VAR)
    ekin = 0.5*V(DENS_VAR)*(V(VELX_VAR)**2 + V(VELY_VAR)**2)
    eint = max(U(ENER_VAR) - ekin, sim_smallPres) !eint=rho*e
    ! get pressure by calling eos
    ! if (eint < 1.E-11) print*, 'small eint detected', eint
    pres = max((sim_gamma-1.)*eint, 1.E-12)
    ! if (pres < 1.1E-12) print *, eint
    V(PRES_VAR) = pres
    V(EINT_VAR) = eint
    V(GAMC_VAR) = sim_gamma
    V(GAME_VAR) = sim_gamma

  end subroutine cons2prim

  subroutine prim2flux(V,Flux,dir)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: Flux
    integer, intent(IN) :: dir

    real :: ekin,eint,ener
    integer :: VEL1_VAR, VEL2_VAR

    select case(dir)
    case(XDIM)
      VEL1_VAR = VELX_VAR
      VEL2_VAR = VELY_VAR
    case(YDIM)
      VEL1_VAR = VELY_VAR
      VEL2_VAR = VELX_VAR
    case DEFAULT
      VEL1_VAR = 0
      VEL2_VAR = 0
      call abort_slug("[prim2flux] Wrong dir value")
    end select


    ekin = 0.5*V(DENS_VAR)*(V(VELX_VAR)**2 + V(VELY_VAR)**2)
    eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    ener = ekin + eint

    Flux(DENS_VAR) = V(DENS_VAR)*V(VEL1_VAR)
    Flux(VEL1_VAR) = Flux(DENS_VAR)*V(VEL1_VAR) + V(PRES_VAR)
    Flux(VEL2_VAR) = Flux(DENS_VAR)*V(VEL2_VAR)
    Flux(ENER_VAR) = V(VEL1_VAR)*(ener + V(PRES_VAR))

  end subroutine prim2flux

  subroutine cons2flux(U, Flux, dir)
    implicit none

    real, dimension(NSYS_VAR), intent(IN) :: U
    real, dimension(NUMB_WAVE), intent(OUT) :: Flux
    integer, intent(IN) :: dir

    real :: pres
    integer :: VEL1_VAR, VEL2_VAR
    real :: vel1, vel2

    select case(dir)
    case(XDIM)
      VEL1_VAR = VELX_VAR
      VEL2_VAR = VELY_VAR
    case(YDIM)
      VEL1_VAR = VELY_VAR
      VEL2_VAR = VELX_VAR
    case DEFAULT
      VEL1_VAR = 0
      VEL2_VAR = 0
      call abort_slug("[cons2flux] Wrong dir value")
    end select

    vel1 = U(VEL1_VAR)/U(DENS_VAR)
    vel2 = U(VEL2_VAR)/U(DENS_VAR)

    pres = (sim_gamma-1.)*( U(ENER_VAR) - 0.5*(U(VEL1_VAR)*vel1 + U(VEL2_VAR)*vel2) )

    ! Flux(DENS_VAR) = U(DENS_VAR)*vel1
    Flux(DENS_VAR) = U(VEL1_VAR)  ! this has better symmetry
    FLUX(VEL1_VAR) = FLUX(DENS_VAR)*vel1 + pres
    FLUX(VEL2_VAR) = FLUX(DENS_VAR)*vel2
    FLUX(ENER_VAR) = vel1*(U(ENER_VAR) + pres)

  end subroutine cons2flux

end module primconsflux

