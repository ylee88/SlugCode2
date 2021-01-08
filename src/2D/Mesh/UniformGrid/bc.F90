module bc

#include "definition.h"

  use grid_data, only: gr_i0,     &
                       gr_imax,   &
                       gr_ibeg,   &
                       gr_iend,   &
                       gr_nx,     &
                       gr_ny,     &
                       gr_ngc,    &
                       gr_xCoord, &
                       gr_yCoord, &
                       gr_xend

  use block_data, only: bl_BC,       &
                        bl_cornerBC, &
                        bl_ID,       &
                        bl_i,        &
                        bl_j,        &
                        bl_iProcs,   &
                        bl_jProcs

  use sim_data,  only: sim_bcTypex,  &
                       sim_bcTypey,  &
                       sim_cornerBC, &
                       sim_t
  use mpi

  implicit none

  contains

  subroutine bc_apply(V)
    ! by calling this function, all processs are ensured to be synced.

    implicit none

    real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM)), intent(INOUT) :: V

    call bc_normal(V)
    if (sim_cornerBC) call bc_corner(V)

  end subroutine bc_apply

  subroutine bc_normal(V)
    implicit none

    real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM)), intent(INOUT) :: V

    real, dimension(NUMB_VAR, gr_ngc, gr_ny) :: loc_buffL, loc_buffR
    real, dimension(NUMB_VAR, gr_nx, gr_ngc) :: loc_buffB, loc_buffT

    real, dimension(NUMB_VAR, gr_ngc, gr_ny) :: rcv_buffL, rcv_buffR
    real, dimension(NUMB_VAR, gr_nx, gr_ngc) :: rcv_buffB, rcv_buffT

    integer :: i, j
    integer :: i0, imax, ibeg, iend
    integer :: j0, jmax, jbeg, jend

    integer :: gc
    integer :: send, recv, stag, rtag, ierr
    integer, dimension(MPI_STATUS_SIZE) :: stat

    !! for DMR
    real :: xmin, ymax, rt3, xx, yy
    rt3 = sqrt(3.)
    ymax = gr_xend(YDIM)
    xmin = 1./6. + 10.*sim_t/(.5*rt3) + ymax/rt3
    !!

    ! for readability
    gc = gr_ngc - 1

    i0 = gr_i0(XDIM)
    imax = gr_imax(XDIM)
    j0 = gr_i0(YDIM)
    jmax = gr_imax(YDIM)

    ibeg = gr_ibeg(XDIM)
    iend = gr_iend(XDIM)
    jbeg = gr_ibeg(YDIM)
    jend = gr_iend(YDIM)

    ! prepare to send data
    ! init local buffers with inner domains
    loc_buffL(:, :, :) = V(:,    ibeg:ibeg+gc,    jbeg:jend   )
    loc_buffR(:, :, :) = V(:, iend-gc:iend,       jbeg:jend   )
    loc_buffB(:, :, :) = V(:,    ibeg:iend,       jbeg:jbeg+gc)
    loc_buffT(:, :, :) = V(:,    ibeg:iend,    jend-gc:jend   )

    stag = 0
    rtag = 0

    ! excange along non-domain boundary
    ! send to left & receive from right
    send = bl_BC(1)
    recv = bl_BC(2)
    call MPI_Sendrecv(loc_buffL, NUMB_VAR*gr_ngc*gr_ny, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffR, NUMB_VAR*gr_ngc*gr_ny, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to right & receive from left
    send = bl_BC(2)
    recv = bl_BC(1)
    call MPI_Sendrecv(loc_buffR, NUMB_VAR*gr_ngc*gr_ny, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffL, NUMB_VAR*gr_ngc*gr_ny, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to bottom & receive from top
    send = bl_BC(3)
    recv = bl_BC(4)
    call MPI_Sendrecv(loc_buffB, NUMB_VAR*gr_nx*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffT, NUMB_VAR*gr_nx*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to top & receive from bottom
    send = bl_BC(4)
    recv = bl_BC(3)
    call MPI_Sendrecv(loc_buffT, NUMB_VAR*gr_nx*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffB, NUMB_VAR*gr_nx*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)

    ! now apply domain BCs
    if (bl_i == 1) then          ! left
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffL(:,i,:) = loc_buffL(:,1,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffL(:,i,:) = loc_buffL(:,gr_ngc-i+1,:)
        end do
        rcv_buffL(VELX_VAR,:,:) = -rcv_buffL(VELX_VAR,:,:)
      else if (sim_bcTypeX == "DMR") then
        ! do nothing; just copy ghost cell's inform
        rcv_buffL(:,:,:) = V(:, i0:i0+gc, jbeg:jend)
      else if (sim_bcTypeX == "SHOCKVORTEX") then
        ! do nothing; just copy ghost cell's inform
        rcv_buffL(:,:,:) = V(:, i0:i0+gc, jbeg:jend)
      end if
    end if

    if (bl_i == bl_iProcs) then  ! right
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffR(:,i,:) = loc_buffR(:,gr_ngc,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffR(:,i,:) = loc_buffR(:,gr_ngc-i+1,:)
        end do
        rcv_buffR(VELX_VAR,:,:) = -rcv_buffR(VELX_VAR,:,:)
      else if (sim_bcTypeX == "DMR") then
        ! do outflow
        do i = 1, gr_ngc
          rcv_buffR(:,i,:) = loc_buffR(:,gr_ngc,:)
        end do
      else if (sim_bcTypeX == "SHOCKVORTEX") then
        ! do outflow
        do i = 1, gr_ngc
          rcv_buffR(:,i,:) = loc_buffR(:,gr_ngc,:)
        end do
      end if
    end if

    if (bl_j == 1) then          ! bottom
      if (sim_bcTypeY == "outflow") then
        do j = 1, gr_ngc
          rcv_buffB(:,:,j) = loc_buffB(:,:,1)
        end do
      else if (sim_bcTypeY == "reflect") then
        do j = 1, gr_ngc
          rcv_buffB(:,:,j) = loc_buffB(:,:,gr_ngc-j+1)
        end do
        rcv_buffB(VELY_VAR,:,:) = -rcv_buffB(VELY_VAR,:,:)
      else if (sim_bcTypeY == "DMR") then
        rcv_buffB(:,:,:) = loc_buffB(:,:,:)   ! copy EOS data
        do i = 1, gr_nx
          xx = gr_xCoord(i+gr_ngc)
          if (xx > 1./6.) then
            ! reflect
            do j = 1, gr_ngc
              rcv_buffB(:,i,j) = loc_buffB(:,i,gr_ngc-j+1)
            end do
            rcv_buffB(VELY_VAR,i,:) = -rcv_buffB(VELY_VAR,i,:)
          else
            ! inside shock
            rcv_buffB(DENS_VAR,i,:) = 8.
            rcv_buffB(VELX_VAR,i,:) = 7.1447096
            rcv_buffB(VELY_VAR,i,:) = -4.125
            rcv_buffB(PRES_VAR,i,:) = 116.5
          end if
        end do ! end DMR
      else if (sim_bcTypeY == "astrojet") then
        ! rcv_buffB(:,:,:) = loc_buffB(:,:,:)   ! copy EOS data
        do j = 1, gr_ngc
          do i = 1, gr_nx
            xx = gr_xCoord(i+gr_ngc)
            if (xx < 0.05 .and. xx > -0.05) then
              ! do nothing
              rcv_buffB(:,i,j) = V(:,i+gr_ngc,j)
            else
              ! do outflow
              rcv_buffB(:,i,j) = loc_buffB(:,i,1)
            end if
          end do
        end do  ! end astrojet
      end if
    end if

    if (bl_j == bl_jProcs) then  ! top
      if (sim_bcTypey == "outflow") then
        do j = 1, gr_ngc
          rcv_buffT(:,:,j) = loc_buffT(:,:,gr_ngc)
        end do
      else if (sim_bcTypey == "reflect") then
        do j = 1, gr_ngc
          rcv_buffT(:,:,j) = loc_buffT(:,:,gr_ngc-j+1)
        end do
        rcv_buffT(VELY_VAR,:,:) = -rcv_buffT(VELY_VAR,:,:)
      else if (sim_bcTypeY == "DMR") then
        rcv_buffT(:,:,:) = loc_buffT(:,:,:)   ! copy EOS data
        do i = 1, gr_nx
          do j = 1, gr_ngc
            xx = gr_xCoord(i+gr_ngc)
            yy = gr_yCoord(j+jend)
            if (yy-ymax > (xx-xmin)*rt3) then
              ! still in shock, use left state
              rcv_buffT(DENS_VAR,i,j) = 8.
              rcv_buffT(VELX_VAR,i,j) = 7.1447096
              rcv_buffT(VELY_VAR,i,j) = -4.125
              rcv_buffT(PRES_VAR,i,j) = 116.5
            else
              ! outside of shock
              rcv_buffT(DENS_VAR,i,j) = 1.4
              rcv_buffT(VELX_VAR,i,j) = 0.
              rcv_buffT(VELY_VAR,i,j) = 0.
              rcv_buffT(PRES_VAR,i,j) = 1.
            end if
          end do
        end do ! end DMR
      else if (sim_bcTypey == "astrojet") then
        ! do outflow
        do j = 1, gr_ngc
          rcv_buffT(:,:,j) = loc_buffT(:,:,gr_ngc)
        end do
      end if
    end if

    ! filling halo
    V(:,      i0:i0+gc,   jbeg:jend ) = rcv_buffL
    V(:, imax-gc:imax,    jbeg:jend ) = rcv_buffR
    V(:,    ibeg:iend,      j0:j0+gc) = rcv_buffB
    V(:,    ibeg:iend, jmax-gc:jmax ) = rcv_buffT

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
  end subroutine bc_normal

  subroutine bc_corner(V)
    implicit none
    real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM)), intent(INOUT) :: V

    real, dimension(NUMB_VAR, gr_ngc, gr_ngc) :: loc_buffTL, loc_buffTR
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc) :: loc_buffBL, loc_buffBR

    real, dimension(NUMB_VAR, gr_ngc, gr_ngc) :: rcv_buffTL, rcv_buffTR
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc) :: rcv_buffBL, rcv_buffBR

    integer :: i, j
    integer :: i0, imax, ibeg, iend
    integer :: j0, jmax, jbeg, jend

    integer :: gc
    integer :: send, recv, stag, rtag, ierr
    integer, dimension(MPI_STATUS_SIZE) :: stat

    !! for DMR
    real :: xmin, ymax, rt3, xx, yy
    rt3 = sqrt(3.)
    ymax = gr_xend(YDIM)
    xmin = 1./6. + 10.*sim_t/(.5*rt3) + ymax/rt3
    !!

    ! for readability
    gc = gr_ngc - 1

    i0 = gr_i0(XDIM)
    imax = gr_imax(XDIM)
    j0 = gr_i0(YDIM)
    jmax = gr_imax(YDIM)

    ibeg = gr_ibeg(XDIM)
    iend = gr_iend(XDIM)
    jbeg = gr_ibeg(YDIM)
    jend = gr_iend(YDIM)

    ! prepare to send data
    ! init local buffers with inner domains
    loc_buffTL(:, :, :) = V(:,    ibeg:ibeg+gc, jend-gc:jend   )
    loc_buffTR(:, :, :) = V(:, iend-gc:iend,    jend-gc:jend   )
    loc_buffBL(:, :, :) = V(:,    ibeg:ibeg+gc,    jbeg:jbeg+gc)
    loc_buffBR(:, :, :) = V(:, iend-gc:iend,       jbeg:jbeg+gc)

    stag = 0
    rtag = 0

    ! excange along non-domain boundary
    ! send to TL & receive from BR
    send = bl_cornerBC(1)
    recv = bl_cornerBC(3)
    call MPI_Sendrecv(loc_buffTL, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBR, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BR & receive from TL
    send = bl_cornerBC(3)
    recv = bl_cornerBC(1)
    call MPI_Sendrecv(loc_buffBR, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTL, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TR & receive from BL
    send = bl_cornerBC(2)
    recv = bl_cornerBC(4)
    call MPI_Sendrecv(loc_buffTR, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBL, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BL & receive from TR
    send = bl_cornerBC(4)
    recv = bl_cornerBC(2)
    call MPI_Sendrecv(loc_buffBL, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTR, NUMB_VAR*gr_ngc*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)

    ! now apply domain BCs
    if (bl_i == 1) then          ! left boundary
      loc_buffTL(:,:,:) = V(:, ibeg:ibeg+gc, jmax-gc:jmax )
      loc_buffBL(:,:,:) = V(:, ibeg:ibeg+gc,      j0:j0+gc)
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffTL(:,i,:) = loc_buffTL(:,1,:)
          rcv_buffBL(:,i,:) = loc_buffBL(:,1,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffTL(:,i,:) = loc_buffTL(:,gr_ngc-i+1,:)
          rcv_buffBL(:,i,:) = loc_buffBL(:,gr_ngc-i+1,:)
        end do
        rcv_buffTL(VELX_VAR,:,:) = -rcv_buffTL(VELX_VAR,:,:)
        rcv_buffBL(VELX_VAR,:,:) = -rcv_buffBL(VELX_VAR,:,:)
      else if (sim_bcTypeX == "DMR") then
        rcv_buffTL(:,:,:) = loc_buffTL(:,:,:)   ! copy EOS data
        rcv_buffBL(:,:,:) = loc_buffBL(:,:,:)   ! copy EOS data
        ! inflow
        rcv_buffTL(DENS_VAR,:,:) = 8.
        rcv_buffTL(VELX_VAR,:,:) = 7.1447096
        rcv_buffTL(VELY_VAR,:,:) = -4.125
        rcv_buffTL(PRES_VAR,:,:) = 116.5

        rcv_buffBL(DENS_VAR,:,:) = 8.
        rcv_buffBL(VELX_VAR,:,:) = 7.1447096
        rcv_buffBL(VELY_VAR,:,:) = -4.125
        rcv_buffBL(PRES_VAR,:,:) = 116.5
      else if (sim_bcTypeX == "SHOCKVORTEX") then
        ! inflow, do nothing
        rcv_buffTL(:,:,:) = V(:, i0:i0+gc, jmax-gc:jmax)
        rcv_buffBL(:,:,:) = V(:, i0:i0+gc, j0:j0+gc)
      end if  ! sim_bcTypeX
    end if    ! bl_i

    if (bl_i == bl_iProcs) then  ! right boundary
      loc_buffTR(:,:,:) = V(:, iend-gc:iend, jmax-gc:jmax )
      loc_buffBR(:,:,:) = V(:, iend-gc:iend,      j0:j0+gc)
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffTR(:,i,:) = loc_buffTR(:,gr_ngc,:)
          rcv_buffBR(:,i,:) = loc_buffBR(:,gr_ngc,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffTR(:,i,:) = loc_buffTR(:,gr_ngc-i+1,:)
          rcv_buffBR(:,i,:) = loc_buffBR(:,gr_ngc-i+1,:)
        end do
        rcv_buffTR(VELX_VAR,:,:) = -rcv_buffTR(VELX_VAR,:,:)
        rcv_buffBR(VELX_VAR,:,:) = -rcv_buffBR(VELX_VAR,:,:)
      else if (sim_bcTypeX == "DMR") then
        ! outflow
        do i = 1, gr_ngc
          rcv_buffTR(:,i,:) = loc_buffTR(:,gr_ngc,:)
          rcv_buffBR(:,i,:) = loc_buffBR(:,gr_ngc,:)
        end do
      else if (sim_bcTypeX == "SHOCKVORTEX") then
        ! outflow
        do i = 1, gr_ngc
          rcv_buffTR(:,i,:) = loc_buffTR(:,gr_ngc,:)
          rcv_buffBR(:,i,:) = loc_buffBR(:,gr_ngc,:)
        end do
      end if  ! sim_bcTypeX
    end if    ! bl_i

    if (bl_j == 1) then          ! bottom boundary
      loc_buffBL(:,:,:) = V(:,      i0:i0+gc, jbeg:jbeg+gc)
      loc_buffBR(:,:,:) = V(:, imax-gc:imax,  jbeg:jbeg+gc)
      if (sim_bcTypeY == "outflow") then
        do j = 1, gr_ngc
          rcv_buffBL(:,:,j) = loc_buffBL(:,:,1)
          rcv_buffBR(:,:,j) = loc_buffBR(:,:,1)
        end do
      else if (sim_bcTypeY == "reflect") then
        do j = 1, gr_ngc
          rcv_buffBL(:,:,j) = loc_buffBL(:,:,gr_ngc-j+1)
          rcv_buffBR(:,:,j) = loc_buffBR(:,:,gr_ngc-j+1)
        end do
        rcv_buffBL(VELY_VAR,:,:) = -rcv_buffBL(VELY_VAR,:,:)
        rcv_buffBR(VELY_VAR,:,:) = -rcv_buffBR(VELY_VAR,:,:)
      else if (sim_bcTypeY == "DMR") then
        rcv_buffBL(:,:,:) = loc_buffBL(:,:,:)   ! copy EOS data
        rcv_buffBR(:,:,:) = loc_buffBR(:,:,:)   ! copy EOS data
        do i = 1, gr_ngc
          ! BL
          if(gr_xCoord(i) > 1./6.) then
            ! do reflect
            do j = 1, gr_ngc
              rcv_buffBL(:,i,j) = loc_buffBL(:,i,gr_ngc-j+1)
            end do
            rcv_buffBL(VELY_VAR,:,:) = -rcv_buffBL(VELY_VAR,:,:)
          else
            ! inside shock
            rcv_buffBL(DENS_VAR,i-gr_ngc,:) = 8.
            rcv_buffBL(VELX_VAR,i-gr_ngc,:) = 7.1447096
            rcv_buffBL(VELY_VAR,i-gr_ngc,:) = -4.125
            rcv_buffBL(PRES_VAR,i-gr_ngc,:) = 116.5
          end if
          ! BR
          if(gr_xCoord(iend+i) > 1./6.) then
            ! do reflect
            do j = 1, gr_ngc
              rcv_buffBR(:,i,j) = loc_buffBR(:,i,gr_ngc-j+1)
            end do
            rcv_buffBR(VELY_VAR,:,:) = -rcv_buffBR(VELY_VAR,:,:)
          else
            ! inside shock
            rcv_buffBR(DENS_VAR,i-gr_ngc,:) = 8.
            rcv_buffBR(VELX_VAR,i-gr_ngc,:) = 7.1447096
            rcv_buffBR(VELY_VAR,i-gr_ngc,:) = -4.125
            rcv_buffBR(PRES_VAR,i-gr_ngc,:) = 116.5
          end if
        end do  ! DMR
      else if (sim_bcTypeY == "astrojet") then
        ! do outflow
        do j = 1, gr_ngc
          rcv_buffBL(:,:,j) = loc_buffBL(:,:,1)
          rcv_buffBR(:,:,j) = loc_buffBR(:,:,1)
        end do
      end if  ! sim_bcTypeY
    end if    ! bl_j

    if (bl_j == bl_jProcs) then  ! top boundary
      loc_buffTL(:,:,:) = V(:,      i0:i0+gc, jend-gc:jend)
      loc_buffTR(:,:,:) = V(:, imax-gc:imax,  jend-gc:jend)
      if (sim_bcTypeY == "outflow") then
        do j = 1, gr_ngc
          rcv_buffTL(:,:,j) = loc_buffTL(:,:,gr_ngc)
          rcv_buffTR(:,:,j) = loc_buffTR(:,:,gr_ngc)
        end do
      else if (sim_bcTypeY == "reflect") then
        do j = 1, gr_ngc
          rcv_buffTL(:,:,j) = loc_buffTL(:,:,gr_ngc-j+1)
          rcv_buffTR(:,:,j) = loc_buffTR(:,:,gr_ngc-j+1)
        end do
        rcv_buffTL(VELY_VAR,:,:) = -rcv_buffTL(VELY_VAR,:,:)
        rcv_buffTR(VELY_VAR,:,:) = -rcv_buffTR(VELY_VAR,:,:)
      else if (sim_bcTypeY == "DMR") then
        rcv_buffTL(:,:,:) = loc_buffTL(:,:,:)   ! copy EOS data
        rcv_buffTR(:,:,:) = loc_buffTR(:,:,:)   ! copy EOS data
        do i = 1, gr_ngc
          do j = 1, gr_ngc
            ! TL
            xx = gr_xCoord(i)
            yy = gr_yCoord(j+jend)
            if (yy-ymax > (xx-xmin)*rt3) then
              ! still in shock, use left state
              rcv_buffTL(DENS_VAR,i-gr_ngc,j) = 8.
              rcv_buffTL(VELX_VAR,i-gr_ngc,j) = 7.1447096
              rcv_buffTL(VELY_VAR,i-gr_ngc,j) = -4.125
              rcv_buffTL(PRES_VAR,i-gr_ngc,j) = 116.5
            else
              ! outside of shock
              rcv_buffTL(DENS_VAR,i-gr_ngc,j) = 1.4
              rcv_buffTL(VELX_VAR,i-gr_ngc,j) = 0.
              rcv_buffTL(VELY_VAR,i-gr_ngc,j) = 0.
              rcv_buffTL(PRES_VAR,i-gr_ngc,j) = 1.
            end if
            ! TR
            xx = gr_xCoord(i+iend)
            yy = gr_yCoord(j+jend)
            if (yy-ymax > (xx-xmin)*rt3) then
              ! still in shock, use left state
              rcv_buffTR(DENS_VAR,i-gr_ngc,j) = 8.
              rcv_buffTR(VELX_VAR,i-gr_ngc,j) = 7.1447096
              rcv_buffTR(VELY_VAR,i-gr_ngc,j) = -4.125
              rcv_buffTR(PRES_VAR,i-gr_ngc,j) = 116.5
            else
              ! outside of shock
              rcv_buffTR(DENS_VAR,i-gr_ngc,j) = 1.4
              rcv_buffTR(VELX_VAR,i-gr_ngc,j) = 0.
              rcv_buffTR(VELY_VAR,i-gr_ngc,j) = 0.
              rcv_buffTR(PRES_VAR,i-gr_ngc,j) = 1.
            end if
          end do
        end do ! end DMR
      else if (sim_bcTypeY == "astrojet") then
        ! do outflow
        do j = 1, gr_ngc
          rcv_buffTL(:,:,j) = loc_buffTL(:,:,gr_ngc)
          rcv_buffTR(:,:,j) = loc_buffTR(:,:,gr_ngc)
        end do
      end if  ! sim_bcTypeY
    end if    ! bl_j

    ! filling halo
    V(:,      i0:i0+gc,      j0:j0+gc) = rcv_buffBL
    V(:, imax-gc:imax,       j0:j0+gc) = rcv_buffBR
    V(:,      i0:i0+gc, jmax-gc:jmax ) = rcv_buffTL
    V(:, imax-gc:imax,  jmax-gc:jmax ) = rcv_buffTR

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
  end subroutine bc_corner


end module bc
