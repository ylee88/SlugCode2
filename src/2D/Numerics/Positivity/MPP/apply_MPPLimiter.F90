subroutine apply_MPPLimiter(dt)

#include "definition.h"

  use primconsflux
  use sim_data, only: sim_gamma
  use eigensystem, only: eigenvalues
  use grid_data, only: gr_dx, gr_dy,     &
                       gr_i0, gr_imax,   &
                       gr_ibeg, gr_iend, &
                       gr_maxSpeed,      &
                       gr_V,             &   ! primitive variables
                       gr_U,             &   ! conserved variables
                       gr_flux               ! high order fluxes

  real, intent(IN) :: dt

  real, EXTERNAL :: get_pres
  real, PARAMETER :: eps = 1.E-12

  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: HOflux
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: LFflux
  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM), NDIM) :: flux

  real, dimension(gr_imax(XDIM), gr_imax(YDIM), NDIM) :: weights     !! output

  real, dimension(gr_imax(XDIM), gr_imax(YDIM), 4) :: amin
  real, dimension(2,2,2,2) :: rescale

  real, dimension(NSYS_VAR, gr_imax(XDIM), gr_imax(YDIM)) :: LFU

  real, dimension(NSYS_VAR) :: local_U, Utmp

  real, dimension(4) :: df, aa, rescale2
  real, dimension(NDIM) :: maxSpd

  real :: pres
  real :: root1, root2, rate
  real :: edens, epres, dfdens_sum
  real :: dtx, dty

  integer :: i, j, m, var, n
  integer :: i1, i2, i3, i4
  integer :: dir, im, jm
  integer :: ibeg, iend, jbeg, jend

  ! for readability
  ibeg = gr_ibeg(XDIM)
  iend = gr_iend(XDIM)
  jbeg = gr_ibeg(YDIM)
  jend = gr_iend(YDIM)

  dtx = dt/gr_dx
  dty = dt/gr_dy

  ! copy high-order flux
  HOflux = gr_flux

  ! calculate Lax-Friedrichs flux
  do dir = XDIM, NDIM

    maxSpd(dir) = maxval(abs(gr_maxSpeed(:,dir)))

    select case(dir)
    case(XDIM)
      im = 1
      jm = 0
    case(YDIM)
      im = 0
      jm = 1
    case DEFAULT
      call abort_slug('[apply_MPPLimiter] wrong dir')
    end select

    ! initialize flux values at cell centers
    do j = jbeg-1, jend
      do i = ibeg-1, iend
        call prim2flux(gr_V(:,i,j), flux(:,i,j, dir), dir)
      end do
    end do

    ! calculate first order Lax-Friedrichs fluxes
    do j = jbeg, jend+1
      do i = ibeg, iend+1
        ! imh
        LFflux(:,i,j,dir) = 0.5*(flux(:,i-im,j-jm,dir) + flux(:,i,j,dir)) &
                           -0.5*maxSpd(dir)*(gr_U(:,i,j) - gr_U(:,i-im,j-jm))
      end do
    end do

  end do


  ! initialize amin to large values
  ! in order to handling boundaries
  amin = 1.E12


  edens = eps
  epres = eps
  do j = jbeg, jend
    do i = ibeg, iend

      LFU(DENS_VAR:ENER_VAR,i,j) = gr_U(DENS_VAR:ENER_VAR,i,j) - &
              dtx*(LFflux(DENS_VAR:ENER_VAR,i+1,j,XDIM) - LFflux(DENS_VAR:ENER_VAR,i,j,XDIM)) - &
              dty*(LFflux(DENS_VAR:ENER_VAR,i,j+1,YDIM) - LFflux(DENS_VAR:ENER_VAR,i,j,YDIM))
      pres = get_pres(LFU(:,i,j))

      edens = MIN(edens, LFU(DENS_VAR,i,j))
      epres = MIN(epres, pres)

    end do
  end do


  ! limiting dens
  do j = jbeg, jend
    do i = ibeg, iend

      ddens = edens - LFU(DENS_VAR,i,j)

      ! (left, right) x (NDIM)
      df(1) =  dtx*(HOflux(DENS_VAR,  i,  j,XDIM) - LFflux(DENS_VAR,  i,  j,XDIM))
      df(2) = -dtx*(HOflux(DENS_VAR,i+1,  j,XDIM) - LFflux(DENS_VAR,i+1,  j,XDIM))
      df(3) =  dty*(HOflux(DENS_VAR,  i,  j,YDIM) - LFflux(DENS_VAR,  i,  j,YDIM))
      df(4) = -dty*(HOflux(DENS_VAR,  i,j+1,YDIM) - LFflux(DENS_VAR,  i,j+1,YDIM))

      dfdens_sum = 0.
      do m = 1, 4
        if (df(m) < 0.) dfdens_sum = dfdens_sum + df(m)
      end do

      do m = 1, 4
        if (df(m) < 0.) then
          amin(i,j,m) = MIN(ddens/(dfdens_sum-eps), 1.)
        else
          amin(i,j,m) = 1.
        end if
      end do

    end do
  end do


  ! limiting pres
  do j = jbeg, jend
    do i = ibeg, iend

      rescale(:,:,:,:) = 1.

      do i4 = 1, 2
        do i3 = 1, 2
          do i2 = 1, 2
            do i1 = 1, 2

              aa(1) = REAL(i1-1)*amin(i,j,1)
              aa(2) = REAL(i2-1)*amin(i,j,2)
              aa(3) = REAL(i3-1)*amin(i,j,3)
              aa(4) = REAL(i4-1)*amin(i,j,4)

              do var = 1, NSYS_VAR
                df(1) = aa(1)*HOflux(var,  i,  j,XDIM) + (1.-aa(1))*LFflux(var,  i,  j,XDIM)
                df(2) = aa(2)*HOflux(var,i+1,  j,XDIM) + (1.-aa(2))*LFflux(var,i+1,  j,XDIM)
                df(3) = aa(3)*HOflux(var,  i,  j,YDIM) + (1.-aa(3))*LFflux(var,  i,  j,YDIM)
                df(4) = aa(4)*HOflux(var,  i,j+1,YDIM) + (1.-aa(4))*LFflux(var,  i,j+1,YDIM)

                local_U(var) = gr_U(var,i,j) + dtx*(df(1) - df(2)) + dty*(df(3) - df(4))
              end do

              pres = get_pres(local_U)

              if (pres < epres) then
                root1 = 0.
                root2 = 1.

                Utmp = 0.

                ! bisection
                do n = 1, 10
                  rate = 0.5*(root1 + root2)

                  Utmp(:) = rate*local_U(:) + (1.-rate)*LFU(:,i,j)
                  pres = get_pres(Utmp)

                  if (pres < 0.) then
                    root2 = rate
                  else
                    root1 = rate
                  end if
                end do

                if (pres > 0.) then
                  rate = rate
                else
                  rate = root1
                end if

                rescale(i1, i2, i3, i4) = rate

              end if

            end do    ! i1
          end do      ! i2
        end do        ! i3
      end do          ! i4


      rescale2 = 1.
      do i3 = 1, 2
        do i2 = 1, 2
          do i1 = 1, 2
            rescale2(1) = MIN(rescale2(1), rescale(2, i1, i2, i3))
            rescale2(2) = MIN(rescale2(2), rescale(i1, 2, i2, i3))
            rescale2(3) = MIN(rescale2(3), rescale(i1, i2, 2, i3))
            rescale2(4) = MIN(rescale2(4), rescale(i1, i2, i3, 2))
          end do
        end do
      end do

      do n = 1, 4
        amin(i, j, n) = rescale2(n)*amin(i, j, n)
      end do

    end do    ! j
  end do      ! i


  do j = jbeg, jend+1
    do i = ibeg, iend+1
      weights(i, j, XDIM) = MIN(amin(i, j, 1), amin(i-1, j, 2))
      weights(i, j, YDIM) = MIN(amin(i, j, 3), amin(i, j-1, 4))
    end do
  end do


  do dir = XDIM, NDIM
    do j = jbeg, jend+1
      do i = ibeg, iend+1
        gr_flux(:, i, j, dir) = weights(i, j, dir)*(HOflux(:, i, j, dir) - LFflux(:, i, j, dir)) + LFflux(:, i, j, dir)

      end do
    end do
  end do


  ! !!! DEBUG
  ! do j = jbeg, jend
  !   do i = ibeg, iend
  !      if (weights(i, j, dir) < 0.) print *, 'theta < 0', weights(i,j,dir)
  !
  !       if (i == 52 .and. j == 9) then
  !         print *, weights(i, j, dir)
  !         print *, 'dens = ', 0. +  &
  !       dtx*(gr_flux(DENS_VAR,i+1,j,XDIM) - gr_flux(DENS_VAR,i,j,XDIM)) - &
  !       dty*(gr_flux(DENS_VAR,i,j+1,YDIM) - gr_flux(DENS_VAR,i,j,YDIM))
  !       end if
  !
  !   end do
  ! end do

  ! print *, sum(weights)



end subroutine apply_MPPLimiter


function get_pres(U) result(pres)

  use sim_data,  only: sim_gamma

  implicit none

  real, dimension(NSYS_VAR), intent(IN) :: U
  real :: ekin, eint, pres

  ekin = 0.5*dot_product(U(MOMX_VAR:MOMY_VAR), U(MOMX_VAR:MOMY_VAR))/U(DENS_VAR)
  eint = U(ENER_VAR) - ekin    ! eint = rho*e

  pres = (sim_gamma-1.)*eint

end function get_pres
