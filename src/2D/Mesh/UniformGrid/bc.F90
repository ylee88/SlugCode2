module bc

#include "definition.h"

  use grid_data, only: gr_i0,   &
                       gr_imax, &
                       gr_ibeg, &
                       gr_iend, &
                       gr_nx,   &
                       gr_ny,   &
                       gr_ngc

  use block_data, only: bl_BC,       &
                        bl_cornerBC, &
                        bl_ID,       &
                        bl_i,        &
                        bl_j,        &
                        bl_iProcs,   &
                        bl_jProcs

  use sim_data,  only: sim_bcTypex,  &
                       sim_bcTypey,  &
                       sim_cornerBC
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
    integer :: send, recv, stag, rtag, stat, ierr

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
      if (sim_bcTypex == "outflow") then
        do i = 1, gr_ngc
          rcv_buffL(:,i,:) = loc_buffL(:,1,:)
        end do
      else if (sim_bcTypex == "reflect") then
        do i = 1, gr_ngc
          rcv_buffL(:,i,:) = loc_buffL(:,gr_ngc-i+1,:)
        end do
        rcv_buffL(VELX_VAR,:,:) = -rcv_buffL(VELX_VAR,:,:)
      end if
    end if

    if (bl_i == bl_iProcs) then  ! right
      if (sim_bcTypex == "outflow") then
        do i = 1, gr_ngc
          rcv_buffR(:,i,:) = loc_buffR(:,gr_ngc,:)
        end do
      else if (sim_bcTypex == "reflect") then
        do i = 1, gr_ngc
          rcv_buffR(:,i,:) = loc_buffR(:,gr_ngc-i+1,:)
        end do
        rcv_buffR(VELX_VAR,:,:) = -rcv_buffR(VELX_VAR,:,:)
      end if
    end if

    if (bl_j == 1) then          ! bottom
      if (sim_bcTypey == "outflow") then
        do j = 1, gr_ngc
          rcv_buffB(:,:,j) = loc_buffB(:,:,1)
        end do
      else if (sim_bcTypey == "reflect") then
        do j = 1, gr_ngc
          rcv_buffB(:,:,j) = loc_buffB(:,:,gr_ngc-j+1)
        end do
        rcv_buffB(VELY_VAR,:,:) = -rcv_buffB(VELY_VAR,:,:)
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
    integer :: send, recv, stag, rtag, stat, ierr

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
