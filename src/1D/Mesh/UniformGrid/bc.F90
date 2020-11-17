module bc

#include "definition.h"

  use grid_data, only: gr_i0,   &
                       gr_imax, &
                       gr_ibeg, &
                       gr_iend, &
                       gr_nx,   &
                       gr_ngc

  use block_data, only: bl_BC,              &
                        bl_ID,              &
                        bl_i,               &
                        bl_iProcs

  use sim_data,  only: sim_bcTypex

  use mpi

  implicit none

  contains

  subroutine bc_apply(V)
    ! by calling this function, all processs are ensured to be synced.

    implicit none

    real, dimension(NUMB_VAR, gr_imax(XDIM)), intent(INOUT) :: V

    call bc_normal(V)

  end subroutine bc_apply

  subroutine bc_normal(V)
    implicit none

    real, dimension(NUMB_VAR, gr_imax(XDIM)), intent(INOUT) :: V

    real, dimension(NUMB_VAR, gr_ngc) :: loc_buffL, loc_buffR

    real, dimension(NUMB_VAR, gr_ngc) :: rcv_buffL, rcv_buffR

    integer :: i
    integer :: i0, imax, ibeg, iend

    integer :: gc
    integer :: send, recv, stag, rtag, ierr
    integer, dimension(MPI_STATUS_SIZE) :: stat

    ! for readability
    gc = gr_ngc - 1

    i0 = gr_i0(XDIM)
    imax = gr_imax(XDIM)

    ibeg = gr_ibeg(XDIM)
    iend = gr_iend(XDIM)

    ! prepare to send data
    ! init local buffers with inner domains
    loc_buffL(:, :) = V(:,    ibeg:ibeg+gc)
    loc_buffR(:, :) = V(:, iend-gc:iend   )

    stag = 0
    rtag = 0

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    ! excange along non-domain boundary
    ! send to left & receive from right
    send = bl_BC(1)
    recv = bl_BC(2)
    call MPI_Sendrecv(loc_buffL, NUMB_VAR*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffR, NUMB_VAR*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to right & receive from left
    send = bl_BC(2)
    recv = bl_BC(1)
    call MPI_Sendrecv(loc_buffR, NUMB_VAR*gr_ngc, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffL, NUMB_VAR*gr_ngc, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)

    ! now apply domain BCs
    if (bl_i == 1) then          ! left
      if (sim_bcTypex == "outflow") then
        do i = 1, gr_ngc
          rcv_buffL(:,i) = loc_buffL(:,1)
        end do
      else if (sim_bcTypex == "reflect") then
        do i = 1, gr_ngc
          rcv_buffL(:,i) = loc_buffL(:,gr_ngc-i+1)
        end do
        rcv_buffL(VELX_VAR,:) = -rcv_buffL(VELX_VAR,:)
      end if
    end if

    if (bl_i == bl_iProcs) then  ! right
      if (sim_bcTypex == "outflow") then
        do i = 1, gr_ngc
          rcv_buffR(:,i) = loc_buffR(:,gr_ngc)
        end do
      else if (sim_bcTypex == "reflect") then
        do i = 1, gr_ngc
          rcv_buffR(:,i) = loc_buffR(:,gr_ngc-i+1)
        end do
        rcv_buffR(VELX_VAR,:) = -rcv_buffR(VELX_VAR,:)
      end if
    end if

    ! filling halo
    V(:,      i0:i0+gc) = rcv_buffL
    V(:, imax-gc:imax ) = rcv_buffR
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
  end subroutine bc_normal

end module bc
