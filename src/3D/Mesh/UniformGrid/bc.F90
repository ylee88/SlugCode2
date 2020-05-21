module bc

#include "definition.h"

  use grid_data, only: gr_i0,   &
                       gr_imax, &
                       gr_ibeg, &
                       gr_iend, &
                       gr_nx,   &
                       gr_ny,   &
                       gr_nz,   &
                       gr_ngc

  use block_data, only: bl_BC,              &
                        bl_corner_cornerBC, &
                        bl_corner_sideBC,   &
                        bl_ID,              &
                        bl_i,               &
                        bl_j,               &
                        bl_k,               &
                        bl_iProcs,          &
                        bl_jProcs,          &
                        bl_kProcs

  use sim_data,  only: sim_bcTypex,  &
                       sim_bcTypey,  &
                       sim_bcTypez,  &
                       sim_cornerBC
  use mpi

  implicit none

  contains

  subroutine bc_apply(V)
    ! by calling this function, all processs are ensured to be synced.

    implicit none

    real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)), intent(INOUT) :: V

    call bc_normal(V)
    if (sim_cornerBC) call bc_corner(V)

  end subroutine bc_apply

  subroutine bc_normal(V)
    implicit none

    real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)), intent(INOUT) :: V

    real, dimension(NUMB_VAR, gr_ngc, gr_ny,  gr_nz ) :: loc_buffL, loc_buffR
    real, dimension(NUMB_VAR, gr_nx,  gr_ngc, gr_nz ) :: loc_buffB, loc_buffT
    real, dimension(NUMB_VAR, gr_nx,  gr_ny,  gr_ngc) :: loc_buffD, loc_buffU

    real, dimension(NUMB_VAR, gr_ngc, gr_ny,  gr_nz ) :: rcv_buffL, rcv_buffR
    real, dimension(NUMB_VAR, gr_nx,  gr_ngc, gr_nz ) :: rcv_buffB, rcv_buffT
    real, dimension(NUMB_VAR, gr_nx,  gr_ny,  gr_ngc) :: rcv_buffD, rcv_buffU

    integer :: i, j, k
    integer :: i0, imax, ibeg, iend
    integer :: j0, jmax, jbeg, jend
    integer :: k0, kmax, kbeg, kend

    integer :: gc
    integer :: send, recv, stag, rtag, stat, ierr

    ! for readability
    gc = gr_ngc - 1

    i0 = gr_i0(XDIM)
    imax = gr_imax(XDIM)
    j0 = gr_i0(YDIM)
    jmax = gr_imax(YDIM)
    k0 = gr_i0(ZDIM)
    kmax = gr_imax(ZDIM)

    ibeg = gr_ibeg(XDIM)
    iend = gr_iend(XDIM)
    jbeg = gr_ibeg(YDIM)
    jend = gr_iend(YDIM)
    kbeg = gr_ibeg(ZDIM)
    kend = gr_iend(ZDIM)

    ! prepare to send data
    ! init local buffers with inner domains
    loc_buffL(:, :, :, :) = V(:,    ibeg:ibeg+gc,    jbeg:jend,       kbeg:kend   )
    loc_buffR(:, :, :, :) = V(:, iend-gc:iend,       jbeg:jend,       kbeg:kend   )
    loc_buffB(:, :, :, :) = V(:,    ibeg:iend,       jbeg:jbeg+gc,    kbeg:kend   )
    loc_buffT(:, :, :, :) = V(:,    ibeg:iend,    jend-gc:jend,       kbeg:kend   )
    loc_buffD(:, :, :, :) = V(:,    ibeg:iend,       jbeg:jend,       kbeg:kbeg+gc)
    loc_buffU(:, :, :, :) = V(:,    ibeg:iend,       jbeg:jend,    kend-gc:kend   )

    stag = 0
    rtag = 0

    ! excange along non-domain boundary
    ! send to left & receive from right
    send = bl_BC(1)
    recv = bl_BC(2)
    call MPI_Sendrecv(loc_buffL, NUMB_VAR*gr_ngc*gr_ny*gr_nz, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffR, NUMB_VAR*gr_ngc*gr_ny*gr_nz, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to right & receive from left
    send = bl_BC(2)
    recv = bl_BC(1)
    call MPI_Sendrecv(loc_buffR, NUMB_VAR*gr_ngc*gr_ny*gr_nz, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffL, NUMB_VAR*gr_ngc*gr_ny*gr_nz, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to bottom & receive from top
    send = bl_BC(3)
    recv = bl_BC(4)
    call MPI_Sendrecv(loc_buffB, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffT, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to top & receive from bottom
    send = bl_BC(4)
    recv = bl_BC(3)
    call MPI_Sendrecv(loc_buffT, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffB, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to down & receive from up
    send = bl_BC(5)
    recv = bl_BC(6)
    call MPI_Sendrecv(loc_buffD, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffU, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to down & receive from up
    send = bl_BC(6)
    recv = bl_BC(5)
    call MPI_Sendrecv(loc_buffU, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffD, NUMB_VAR*gr_nx*gr_ngc*gr_nz, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)

    ! now apply domain BCs
    if (bl_i == 1) then          ! left
      if (sim_bcTypex == "outflow") then
        do i = 1, gr_ngc
          rcv_buffL(:,i,:,:) = loc_buffL(:,1,:,:)
        end do
      else if (sim_bcTypex == "reflect") then
        do i = 1, gr_ngc
          rcv_buffL(:,i,:,:) = loc_buffL(:,gr_ngc-i+1,:,:)
        end do
        rcv_buffL(VELX_VAR,:,:,:) = -rcv_buffL(VELX_VAR,:,:,:)
      end if
    end if

    if (bl_i == bl_iProcs) then  ! right
      if (sim_bcTypex == "outflow") then
        do i = 1, gr_ngc
          rcv_buffR(:,i,:,:) = loc_buffR(:,gr_ngc,:,:)
        end do
      else if (sim_bcTypex == "reflect") then
        do i = 1, gr_ngc
          rcv_buffR(:,i,:,:) = loc_buffR(:,gr_ngc-i+1,:,:)
        end do
        rcv_buffR(VELX_VAR,:,:,:) = -rcv_buffR(VELX_VAR,:,:,:)
      end if
    end if

    if (bl_j == 1) then          ! bottom
      if (sim_bcTypey == "outflow") then
        do j = 1, gr_ngc
          rcv_buffB(:,:,j,:) = loc_buffB(:,:,1,:)
        end do
      else if (sim_bcTypey == "reflect") then
        do j = 1, gr_ngc
          rcv_buffB(:,:,j,:) = loc_buffB(:,:,gr_ngc-j+1,:)
        end do
        rcv_buffB(VELY_VAR,:,:,:) = -rcv_buffB(VELY_VAR,:,:,:)
      end if
    end if

    if (bl_j == bl_jProcs) then  ! top
      if (sim_bcTypey == "outflow") then
        do j = 1, gr_ngc
          rcv_buffT(:,:,j,:) = loc_buffT(:,:,gr_ngc,:)
        end do
      else if (sim_bcTypey == "reflect") then
        do j = 1, gr_ngc
          rcv_buffT(:,:,j,:) = loc_buffT(:,:,gr_ngc-j+1,:)
        end do
        rcv_buffT(VELY_VAR,:,:,:) = -rcv_buffT(VELY_VAR,:,:,:)
      end if
    end if

    if (bl_k == 1) then  ! down
      if (sim_bcTypez == "outflow") then
        do k = 1, gr_ngc
          rcv_buffD(:,:,:,k) = loc_buffD(:,:,:,1)
        end do
      else if (sim_bcTypey == "reflect") then
        do k = 1, gr_ngc
          rcv_buffD(:,:,:,k) = loc_buffD(:,:,:,gr_ngc-k+1)
        end do
        rcv_buffD(VELZ_VAR,:,:,:) = -rcv_buffD(VELZ_VAR,:,:,:)
      end if
    end if

    if (bl_k == bl_kProcs) then  ! up
      if (sim_bcTypez == "outflow") then
        do k = 1, gr_ngc
          rcv_buffU(:,:,:,k) = loc_buffU(:,:,:,gr_ngc)
        end do
      else if (sim_bcTypey == "reflect") then
        do k = 1, gr_ngc
          rcv_buffU(:,:,:,k) = loc_buffU(:,:,:,gr_ngc-k+1)
        end do
        rcv_buffU(VELZ_VAR,:,:,:) = -rcv_buffU(VELZ_VAR,:,:,:)
      end if
    end if

    ! filling halo
    V(:,      i0:i0+gc,   jbeg:jend,     kbeg:kend ) = rcv_buffL
    V(:, imax-gc:imax,    jbeg:jend,     kbeg:kend ) = rcv_buffR
    V(:,    ibeg:iend,      j0:j0+gc,    kbeg:kend ) = rcv_buffB
    V(:,    ibeg:iend, jmax-gc:jmax,     kbeg:kend ) = rcv_buffT
    V(:,    ibeg:iend,    jbeg:jend,       k0:k0+gc) = rcv_buffD
    V(:,    ibeg:iend,    jbeg:jend,  kmax-gc:kmax ) = rcv_buffU
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
  end subroutine bc_normal

  subroutine bc_corner(V)
    implicit none
    real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)), intent(INOUT) :: V

    ! 12 rectangle buffers for side-coners
    real, dimension(NUMB_VAR, gr_nx , gr_ngc, gr_ngc) :: loc_buffBS, loc_buffBN, loc_buffTS, loc_buffTN
    real, dimension(NUMB_VAR, gr_ngc, gr_ny , gr_ngc) :: loc_buffDW, loc_buffDE, loc_buffUW, loc_buffUE
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc, gr_nz ) :: loc_buffBW, loc_buffBE, loc_buffTW, loc_buffTE
    ! receive buffer
    real, dimension(NUMB_VAR, gr_nx , gr_ngc, gr_ngc) :: rcv_buffBS, rcv_buffBN, rcv_buffTS, rcv_buffTN
    real, dimension(NUMB_VAR, gr_ngc, gr_ny , gr_ngc) :: rcv_buffDW, rcv_buffDE, rcv_buffUW, rcv_buffUE
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc, gr_nz ) :: rcv_buffBW, rcv_buffBE, rcv_buffTW, rcv_buffTE

    ! 8 cubic buffers for corner-corners
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc, gr_ngc) :: loc_buffBSW, loc_buffBSE, loc_buffBNW, loc_buffBNE
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc, gr_ngc) :: loc_buffTSW, loc_buffTSE, loc_buffTNW, loc_buffTNE
    ! receive buffers
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc, gr_ngc) :: rcv_buffBSW, rcv_buffBSE, rcv_buffBNW, rcv_buffBNE
    real, dimension(NUMB_VAR, gr_ngc, gr_ngc, gr_ngc) :: rcv_buffTSW, rcv_buffTSE, rcv_buffTNW, rcv_buffTNE

    integer :: i, j, k
    integer :: i0, imax, ibeg, iend
    integer :: j0, jmax, jbeg, jend
    integer :: k0, kmax, kbeg, kend

    integer :: gc
    integer :: send, recv, stag, rtag, stat, ierr, ndata

    ! for readability
    gc = gr_ngc - 1

    i0 = gr_i0(XDIM)
    imax = gr_imax(XDIM)
    j0 = gr_i0(YDIM)
    jmax = gr_imax(YDIM)
    k0 = gr_i0(ZDIM)
    kmax = gr_imax(ZDIM)

    ibeg = gr_ibeg(XDIM)
    iend = gr_iend(XDIM)
    jbeg = gr_ibeg(YDIM)
    jend = gr_iend(YDIM)
    kbeg = gr_ibeg(ZDIM)
    kend = gr_iend(ZDIM)

    ! prepare to send data
    ! init local buffers with inner domains

    ! 12 rectangle buffers for side-corners
    ! bl_corner_sideBC :: side-corners, ngc*ngc*gr_n[x,y,z] rectangles
    !             *------12------*
    !            / |            /|
    !           1  |           2 |
    !          /   |          /  |
    !         *----+--11-----*   |
    !         |    7         |   8
    !         |    |         |   |
    !         5    |         6   |
    !         |    *----10---+---*        j    k
    !         |   /          |  /         ^  7
    !         |  3           | 4          | /
    !         | /            |/           |/
    !         *-------9------*            +-----> i
    loc_buffBW(:,:,:,:) = V(:,    ibeg:ibeg+gc,    jbeg:jbeg+gc,    kbeg:kend   ) !3
    loc_buffBE(:,:,:,:) = V(:, iend-gc:iend,       jbeg:jbeg+gc,    kbeg:kend   ) !4
    loc_buffBS(:,:,:,:) = V(:,    ibeg:iend,       jbeg:jbeg+gc,    kbeg:kbeg+gc) !9
    loc_buffBN(:,:,:,:) = V(:,    ibeg:iend,       jbeg:jbeg+gc, kend-gc:kend   ) !10

    loc_buffTW(:,:,:,:) = V(:,    ibeg:ibeg+gc, jend-gc:jend,       kbeg:kend   ) !1
    loc_buffTE(:,:,:,:) = V(:, iend-gc:iend,    jend-gc:jend,       kbeg:kend   ) !2
    loc_buffTS(:,:,:,:) = V(:,    ibeg:iend,    jend-gc:jend,       kbeg:kbeg+gc  ) !11
    loc_buffTN(:,:,:,:) = V(:,    ibeg:iend,    jend-gc:jend,    kend-gc:kend     ) !12

    loc_buffDW(:,:,:,:) = V(:,    ibeg:ibeg+gc,    jbeg:jend,       kbeg:kbeg+gc  ) !5
    loc_buffDE(:,:,:,:) = V(:, iend-gc:iend,       jbeg:jend,       kbeg:kbeg+gc  ) !6
    loc_buffUW(:,:,:,:) = V(:,    ibeg:ibeg+gc,    jbeg:jend,    kend-gc:kend     ) !7
    loc_buffUE(:,:,:,:) = V(:, iend-gc:iend,       jbeg:jend,    kend-gc:kend     ) !8

    ! 8 cubic buffers for corner-corners
    ! bl_corner_cornerBC :: corner-corners, ngc*ngc*ngc cubes
    !             7--------------8
    !            / |            /|
    !           /  |           / |
    !          /   |          /  |
    !         5----+---------6   |
    !         |    |         |   |
    !         |    |         |   |
    !         |    |         |   |
    !         |    3---------+---4        j    k
    !         |   /          |  /         ^  7
    !         |  /           | /          | /
    !         | /            |/           |/
    !         1--------------2            +-----> i
    loc_buffBSW(:,:,:,:) = V(:, ibeg   :ibeg+gc, jbeg:jbeg+gc, kbeg   :kbeg+gc) !1
    loc_buffBSE(:,:,:,:) = V(:, iend-gc:iend,    jbeg:jbeg+gc, kbeg   :kbeg+gc) !2
    loc_buffBNW(:,:,:,:) = V(:, ibeg   :ibeg+gc, jbeg:jbeg+gc, kend-gc:kend   ) !3
    loc_buffBNE(:,:,:,:) = V(:, iend-gc:iend,    jbeg:jbeg+gc, kend-gc:kend   ) !4

    loc_buffTSW(:,:,:,:) = V(:, ibeg   :ibeg+gc, jend-gc:jend, kbeg   :kbeg+gc) !5
    loc_buffTSE(:,:,:,:) = V(:, iend-gc:iend,    jend-gc:jend, kbeg   :kbeg+gc) !6
    loc_buffTNW(:,:,:,:) = V(:, ibeg   :ibeg+gc, jend-gc:jend, kend-gc:kend   ) !7
    loc_buffTNE(:,:,:,:) = V(:, iend-gc:iend,    jend-gc:jend, kend-gc:kend   ) !8

    stag = 0
    rtag = 0

    ! excange along non-domain boundary
    ! corner-sides
    ! z-directional rectangles
    ndata = NUMB_VAR*gr_ngc*gr_ngc*gr_nz
    ! send to TW & receive from BE
    send = bl_corner_sideBC(1)
    recv = bl_corner_sideBC(4)
    call MPI_Sendrecv(loc_buffTW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BE & receive from TW
    send = bl_corner_sideBC(4)
    recv = bl_corner_sideBC(1)
    call MPI_Sendrecv(loc_buffBE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TE & receive from BW
    send = bl_corner_sideBC(2)
    recv = bl_corner_sideBC(3)
    call MPI_Sendrecv(loc_buffTE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BW & receive from TE
    send = bl_corner_sideBC(3)
    recv = bl_corner_sideBC(2)
    call MPI_Sendrecv(loc_buffBW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)

    ! y-directional rectangles
    ndata = NUMB_VAR*gr_ngc*gr_ngc*gr_ny
    ! send to DW & receive from UE
    send = bl_corner_sideBC(5)
    recv = bl_corner_sideBC(8)
    call MPI_Sendrecv(loc_buffDW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffUE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to UE & receive from DW
    send = bl_corner_sideBC(8)
    recv = bl_corner_sideBC(5)
    call MPI_Sendrecv(loc_buffUE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffDW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to DE & receive from UW
    send = bl_corner_sideBC(6)
    recv = bl_corner_sideBC(7)
    call MPI_Sendrecv(loc_buffDE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffUW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to UW & receive from DE
    send = bl_corner_sideBC(7)
    recv = bl_corner_sideBC(6)
    call MPI_Sendrecv(loc_buffUW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffDE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)

    ! x-directional rectangles
    ndata = NUMB_VAR*gr_ngc*gr_ngc*gr_nx
    ! send to BS & receive from TN
    send = bl_corner_sideBC(9)
    recv = bl_corner_sideBC(12)
    call MPI_Sendrecv(loc_buffBS, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTN, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TN & receive from BS
    send = bl_corner_sideBC(12)
    recv = bl_corner_sideBC(9)
    call MPI_Sendrecv(loc_buffTN, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBS, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BN & receive from TS
    send = bl_corner_sideBC(10)
    recv = bl_corner_sideBC(11)
    call MPI_Sendrecv(loc_buffBN, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTS, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TS & receive from BN
    send = bl_corner_sideBC(11)
    recv = bl_corner_sideBC(10)
    call MPI_Sendrecv(loc_buffTS, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBN, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    !!! done: corner-side exchanges

    ! exchange corner-corners
    ndata = NUMB_VAR*gr_ngc*gr_ngc*gr_ngc
    ! send to BSW & receive from TNE
    send = bl_corner_cornerBC(1)
    recv = bl_corner_cornerBC(8)
    call MPI_Sendrecv(loc_buffBSW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTNE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TNE & receive from BSW
    send = bl_corner_cornerBC(8)
    recv = bl_corner_cornerBC(1)
    call MPI_Sendrecv(loc_buffTNE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBSW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BSE & receive from TNW
    send = bl_corner_cornerBC(2)
    recv = bl_corner_cornerBC(7)
    call MPI_Sendrecv(loc_buffBSE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTNW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TNW & receive from BSE
    send = bl_corner_cornerBC(7)
    recv = bl_corner_cornerBC(2)
    call MPI_Sendrecv(loc_buffTNW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBSE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BNW & receive from TSE
    send = bl_corner_cornerBC(3)
    recv = bl_corner_cornerBC(6)
    call MPI_Sendrecv(loc_buffBNW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTSE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TSE & receive from BNW
    send = bl_corner_cornerBC(6)
    recv = bl_corner_cornerBC(3)
    call MPI_Sendrecv(loc_buffTSE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBNW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to BNE & receive from TSW
    send = bl_corner_cornerBC(4)
    recv = bl_corner_cornerBC(5)
    call MPI_Sendrecv(loc_buffBNE, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffTSW, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)
    ! send to TSW & receive from BNE
    send = bl_corner_cornerBC(5)
    recv = bl_corner_cornerBC(4)
    call MPI_Sendrecv(loc_buffTSW, ndata, MPI_DOUBLE_PRECISION, send, stag, &
                      rcv_buffBNE, ndata, MPI_DOUBLE_PRECISION, recv, rtag, &
                      MPI_COMM_WORLD, stat, ierr)

    ! now apply domain BCs
    ! rectangles first
    if (bl_i == 1) then          ! left boundary
      loc_buffTW(:,:,:,:) = V(:, ibeg:ibeg+gc, jmax-gc:jmax,    kbeg:kend ) !1
      loc_buffBW(:,:,:,:) = V(:, ibeg:ibeg+gc,      j0:j0+gc,   kbeg:kend ) !3
      loc_buffDW(:,:,:,:) = V(:, ibeg:ibeg+gc,    jbeg:jend,      k0:k0+gc) !5
      loc_buffUW(:,:,:,:) = V(:, ibeg:ibeg+gc,    jbeg:jend, kmax-gc:kmax ) !7
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffTW(:,i,:,:) = loc_buffTW(:,1,:,:)
          rcv_buffBW(:,i,:,:) = loc_buffBW(:,1,:,:)
          rcv_buffDW(:,i,:,:) = loc_buffDW(:,1,:,:)
          rcv_buffUW(:,i,:,:) = loc_buffUW(:,1,:,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffTW(:,i,:,:) = loc_buffTW(:,gr_ngc-i+1,:,:)
          rcv_buffBW(:,i,:,:) = loc_buffBW(:,gr_ngc-i+1,:,:)
          rcv_buffDW(:,i,:,:) = loc_buffDW(:,gr_ngc-i+1,:,:)
          rcv_buffUW(:,i,:,:) = loc_buffUW(:,gr_ngc-i+1,:,:)
        end do
        rcv_buffTW(VELX_VAR,:,:,:) = -rcv_buffTW(VELX_VAR,:,:,:)
        rcv_buffBW(VELX_VAR,:,:,:) = -rcv_buffBW(VELX_VAR,:,:,:)
        rcv_buffDW(VELX_VAR,:,:,:) = -rcv_buffDW(VELX_VAR,:,:,:)
        rcv_buffUW(VELX_VAR,:,:,:) = -rcv_buffUW(VELX_VAR,:,:,:)
      end if  ! sim_bcTypeX
    end if    ! bl_i

    if (bl_i == bl_iProcs) then  ! right boundary
      loc_buffTE(:,:,:,:) = V(:, iend-gc:iend, jmax-gc:jmax,    kbeg:kend ) !2
      loc_buffBE(:,:,:,:) = V(:, iend-gc:iend,      j0:j0+gc,   kbeg:kend ) !4
      loc_buffDE(:,:,:,:) = V(:, iend-gc:iend,    jbeg:jend,      k0:k0+gc) !6
      loc_buffUE(:,:,:,:) = V(:, iend-gc:iend,    jbeg:jend, kmax-gc:kmax ) !8
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffTE(:,i,:,:) = loc_buffTE(:,gr_ngc,:,:)
          rcv_buffBE(:,i,:,:) = loc_buffBE(:,gr_ngc,:,:)
          rcv_buffDE(:,i,:,:) = loc_buffDE(:,gr_ngc,:,:)
          rcv_buffUE(:,i,:,:) = loc_buffUE(:,gr_ngc,:,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffTE(:,i,:,:) = loc_buffTE(:,gr_ngc-i+1,:,:)
          rcv_buffBE(:,i,:,:) = loc_buffBE(:,gr_ngc-i+1,:,:)
          rcv_buffDE(:,i,:,:) = loc_buffDE(:,gr_ngc-i+1,:,:)
          rcv_buffUE(:,i,:,:) = loc_buffUE(:,gr_ngc-i+1,:,:)
        end do
        rcv_buffTE(VELX_VAR,:,:,:) = -rcv_buffTE(VELX_VAR,:,:,:)
        rcv_buffBE(VELX_VAR,:,:,:) = -rcv_buffBE(VELX_VAR,:,:,:)
        rcv_buffDE(VELX_VAR,:,:,:) = -rcv_buffDE(VELX_VAR,:,:,:)
        rcv_buffUE(VELX_VAR,:,:,:) = -rcv_buffUE(VELX_VAR,:,:,:)
      end if  ! sim_bcTypeX
    end if    ! bl_i

    if (bl_j == 1) then          ! bottom boundary
      loc_buffBW(:,:,:,:) = V(:,      i0:i0+gc, jbeg:jbeg+gc,    kbeg:kend ) !3
      loc_buffBE(:,:,:,:) = V(:, imax-gc:imax,  jbeg:jbeg+gc,    kbeg:kend ) !4
      loc_buffBS(:,:,:,:) = V(:,    ibeg:iend,  jbeg:jbeg+gc,      k0:k0+gc) !9
      loc_buffBN(:,:,:,:) = V(:,    ibeg:iend,  jbeg:jbeg+gc, kmax-gc:kmax ) !10
      if (sim_bcTypeY == "outflow") then
        do j = 1, gr_ngc
          rcv_buffBW(:,:,j,:) = loc_buffBW(:,:,1,:)
          rcv_buffBE(:,:,j,:) = loc_buffBE(:,:,1,:)
          rcv_buffBS(:,:,j,:) = loc_buffBS(:,:,1,:)
          rcv_buffBN(:,:,j,:) = loc_buffBN(:,:,1,:)
        end do
      else if (sim_bcTypeY == "reflect") then
        do j = 1, gr_ngc
          rcv_buffBW(:,:,j,:) = loc_buffBW(:,:,gr_ngc-j+1,:)
          rcv_buffBE(:,:,j,:) = loc_buffBE(:,:,gr_ngc-j+1,:)
          rcv_buffBS(:,:,j,:) = loc_buffBS(:,:,gr_ngc-j+1,:)
          rcv_buffBN(:,:,j,:) = loc_buffBN(:,:,gr_ngc-j+1,:)
        end do
        rcv_buffBW(VELY_VAR,:,:,:) = -rcv_buffBW(VELY_VAR,:,:,:)
        rcv_buffBE(VELY_VAR,:,:,:) = -rcv_buffBE(VELY_VAR,:,:,:)
        rcv_buffBS(VELY_VAR,:,:,:) = -rcv_buffBS(VELY_VAR,:,:,:)
        rcv_buffBN(VELY_VAR,:,:,:) = -rcv_buffBN(VELY_VAR,:,:,:)
      end if  ! sim_bcTypeY
    end if    ! bl_j

    if (bl_j == bl_jProcs) then  ! top boundary
      loc_buffTW(:,:,:,:) = V(:,      i0:i0+gc, jend-gc:jend,    kbeg:kend ) !1
      loc_buffTE(:,:,:,:) = V(:, imax-gc:imax,  jend-gc:jend,    kbeg:kend ) !2
      loc_buffTS(:,:,:,:) = V(:,    ibeg:iend,  jend-gc:jend,      k0:k0+gc) !11
      loc_buffTN(:,:,:,:) = V(:,    ibeg:iend,  jend-gc:jend, kmax-gc:kmax ) !12
      if (sim_bcTypeY == "outflow") then
        do j = 1, gr_ngc
          rcv_buffTW(:,:,j,:) = loc_buffTW(:,:,gr_ngc,:)
          rcv_buffTE(:,:,j,:) = loc_buffTE(:,:,gr_ngc,:)
          rcv_buffTS(:,:,j,:) = loc_buffTS(:,:,gr_ngc,:)
          rcv_buffTN(:,:,j,:) = loc_buffTN(:,:,gr_ngc,:)
        end do
      else if (sim_bcTypeY == "reflect") then
        do j = 1, gr_ngc
          rcv_buffTW(:,:,j,:) = loc_buffTW(:,:,gr_ngc-j+1,:)
          rcv_buffTE(:,:,j,:) = loc_buffTE(:,:,gr_ngc-j+1,:)
          rcv_buffTS(:,:,j,:) = loc_buffTS(:,:,gr_ngc-j+1,:)
          rcv_buffTN(:,:,j,:) = loc_buffTN(:,:,gr_ngc-j+1,:)
        end do
        rcv_buffTW(VELY_VAR,:,:,:) = -rcv_buffTW(VELY_VAR,:,:,:)
        rcv_buffTE(VELY_VAR,:,:,:) = -rcv_buffTE(VELY_VAR,:,:,:)
        rcv_buffTS(VELY_VAR,:,:,:) = -rcv_buffTS(VELY_VAR,:,:,:)
        rcv_buffTN(VELY_VAR,:,:,:) = -rcv_buffTN(VELY_VAR,:,:,:)
      end if  ! sim_bcTypeY
    end if    ! bl_j

    if (bl_k == 1) then  ! down boundary
      loc_buffDW(:,:,:,:) = V(:,      i0:i0+gc,    jbeg:jend,  kbeg:kbeg+gc) !5
      loc_buffDE(:,:,:,:) = V(:, imax-gc:imax,     jbeg:jend,  kbeg:kbeg+gc) !6
      loc_buffBS(:,:,:,:) = V(:,    ibeg:iend,       j0:j0+gc, kbeg:kbeg+gc) !9
      loc_buffTS(:,:,:,:) = V(:,    ibeg:iend,  jmax-gc:jmax,  kbeg:kbeg+gc) !11
      if (sim_bcTypeZ == "outflow") then
        do k = 1, gr_ngc
          rcv_buffDW(:,:,:,k) = loc_buffDW(:,:,:,gr_ngc)
          rcv_buffDE(:,:,:,k) = loc_buffDE(:,:,:,gr_ngc)
          rcv_buffBS(:,:,:,k) = loc_buffBS(:,:,:,gr_ngc)
          rcv_buffTS(:,:,:,k) = loc_buffTS(:,:,:,gr_ngc)
        end do
      else if (sim_bcTypeZ == "reflect") then
        do k = 1, gr_ngc
          rcv_buffDW(:,:,:,k) = loc_buffDW(:,:,:,gr_ngc-k+1)
          rcv_buffDE(:,:,:,k) = loc_buffDE(:,:,:,gr_ngc-k+1)
          rcv_buffBS(:,:,:,k) = loc_buffBS(:,:,:,gr_ngc-k+1)
          rcv_buffTS(:,:,:,k) = loc_buffTS(:,:,:,gr_ngc-k+1)
        end do
        rcv_buffDW(VELZ_VAR,:,:,:) = -rcv_buffDW(VELZ_VAR,:,:,:)
        rcv_buffDE(VELZ_VAR,:,:,:) = -rcv_buffDE(VELZ_VAR,:,:,:)
        rcv_buffBS(VELZ_VAR,:,:,:) = -rcv_buffBS(VELZ_VAR,:,:,:)
        rcv_buffTS(VELZ_VAR,:,:,:) = -rcv_buffTS(VELZ_VAR,:,:,:)
      end if  ! sim_bcTypeZ
    end if    ! bl_k

    if (bl_k == bl_kProcs) then  ! up boundary
      loc_buffUW(:,:,:,:) = V(:,      i0:i0+gc,    jbeg:jend,  kend-gc:kend) !7
      loc_buffUE(:,:,:,:) = V(:, imax-gc:imax,     jbeg:jend,  kend-gc:kend) !8
      loc_buffBN(:,:,:,:) = V(:,    ibeg:iend,       j0:j0+gc, kend-gc:kend) !10
      loc_buffTN(:,:,:,:) = V(:,    ibeg:iend,  jmax-gc:jmax,  kend-gc:kend) !12
      if (sim_bcTypeZ == "outflow") then
        do k = 1, gr_ngc
          rcv_buffUW(:,:,:,k) = loc_buffUW(:,:,:,gr_ngc)
          rcv_buffUE(:,:,:,k) = loc_buffUE(:,:,:,gr_ngc)
          rcv_buffBN(:,:,:,k) = loc_buffBN(:,:,:,gr_ngc)
          rcv_buffTN(:,:,:,k) = loc_buffTN(:,:,:,gr_ngc)
        end do
      else if (sim_bcTypeZ == "reflect") then
        do k = 1, gr_ngc
          rcv_buffUW(:,:,:,k) = loc_buffUW(:,:,:,gr_ngc-k+1)
          rcv_buffUE(:,:,:,k) = loc_buffUE(:,:,:,gr_ngc-k+1)
          rcv_buffBN(:,:,:,k) = loc_buffBN(:,:,:,gr_ngc-k+1)
          rcv_buffTN(:,:,:,k) = loc_buffTN(:,:,:,gr_ngc-k+1)
        end do
        rcv_buffUW(VELZ_VAR,:,:,:) = -rcv_buffUW(VELZ_VAR,:,:,:)
        rcv_buffUE(VELZ_VAR,:,:,:) = -rcv_buffUE(VELZ_VAR,:,:,:)
        rcv_buffBN(VELZ_VAR,:,:,:) = -rcv_buffBN(VELZ_VAR,:,:,:)
        rcv_buffTN(VELZ_VAR,:,:,:) = -rcv_buffTN(VELZ_VAR,:,:,:)
      end if  ! sim_bcTypeZ
    end if    ! bl_k
    !!! done: domain BCs for rectangles

    ! now lets do the cubes for domain BCs
    if (bl_i == 1) then  ! left boundary
      loc_buffBSW(:,:,:,:) = rcv_buffBS(:,1:gr_ngc,:,:)
      loc_buffBNW(:,:,:,:) = rcv_buffBN(:,1:gr_ngc,:,:)
      loc_buffTSW(:,:,:,:) = rcv_buffTS(:,1:gr_ngc,:,:)
      loc_buffTNW(:,:,:,:) = rcv_buffTN(:,1:gr_ngc,:,:)
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffBSW(:,i,:,:) = loc_buffBSW(:,1,:,:)
          rcv_buffBNW(:,i,:,:) = loc_buffBNW(:,1,:,:)
          rcv_buffTSW(:,i,:,:) = loc_buffTSW(:,1,:,:)
          rcv_buffTNW(:,i,:,:) = loc_buffTNW(:,1,:,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffBSW(:,i,:,:) = loc_buffBSW(:,gr_ngc-i+1,:,:)
          rcv_buffBNW(:,i,:,:) = loc_buffBNW(:,gr_ngc-i+1,:,:)
          rcv_buffTSW(:,i,:,:) = loc_buffTSW(:,gr_ngc-i+1,:,:)
          rcv_buffTNW(:,i,:,:) = loc_buffTNW(:,gr_ngc-i+1,:,:)
        end do
        rcv_buffBSW(VELX_VAR,:,:,:) = -rcv_buffBSW(VELX_VAR,:,:,:)
        rcv_buffBNW(VELX_VAR,:,:,:) = -rcv_buffBNW(VELX_VAR,:,:,:)
        rcv_buffTSW(VELX_VAR,:,:,:) = -rcv_buffTSW(VELX_VAR,:,:,:)
        rcv_buffTNW(VELX_VAR,:,:,:) = -rcv_buffTNW(VELX_VAR,:,:,:)
      end if  ! sim_bcTypeX
    end if    ! bl_i

    if (bl_i == bl_iProcs) then  ! right boundary
      loc_buffBSE(:,:,:,:) = rcv_buffBS(:,gr_nx-gc:gr_nx,:,:)
      loc_buffBNE(:,:,:,:) = rcv_buffBN(:,gr_nx-gc:gr_nx,:,:)
      loc_buffTSE(:,:,:,:) = rcv_buffTS(:,gr_nx-gc:gr_nx,:,:)
      loc_buffTNE(:,:,:,:) = rcv_buffTN(:,gr_nx-gc:gr_nx,:,:)
      if (sim_bcTypeX == "outflow") then
        do i = 1, gr_ngc
          rcv_buffBSE(:,i,:,:) = loc_buffBSE(:,gr_ngc,:,:)
          rcv_buffBNE(:,i,:,:) = loc_buffBNE(:,gr_ngc,:,:)
          rcv_buffTSE(:,i,:,:) = loc_buffTSE(:,gr_ngc,:,:)
          rcv_buffTNE(:,i,:,:) = loc_buffTNE(:,gr_ngc,:,:)
        end do
      else if (sim_bcTypeX == "reflect") then
        do i = 1, gr_ngc
          rcv_buffBSE(:,i,:,:) = loc_buffBSE(:,gr_ngc-i+1,:,:)
          rcv_buffBNE(:,i,:,:) = loc_buffBNE(:,gr_ngc-i+1,:,:)
          rcv_buffTSE(:,i,:,:) = loc_buffTSE(:,gr_ngc-i+1,:,:)
          rcv_buffTNE(:,i,:,:) = loc_buffTNE(:,gr_ngc-i+1,:,:)
        end do
        rcv_buffBSE(VELX_VAR,:,:,:) = -rcv_buffBSE(VELX_VAR,:,:,:)
        rcv_buffBNE(VELX_VAR,:,:,:) = -rcv_buffBNE(VELX_VAR,:,:,:)
        rcv_buffTSE(VELX_VAR,:,:,:) = -rcv_buffTSE(VELX_VAR,:,:,:)
        rcv_buffTNE(VELX_VAR,:,:,:) = -rcv_buffTNE(VELX_VAR,:,:,:)
      end if  ! sim_bcTypeX
    end if    ! bl_i

    if (bl_j == 1) then   ! bottom boundary
      loc_buffBSW(:,:,:,:) = rcv_buffDW(:,:,1:gr_ngc,:)
      loc_buffBSE(:,:,:,:) = rcv_buffDE(:,:,1:gr_ngc,:)
      loc_buffBNW(:,:,:,:) = rcv_buffUW(:,:,1:gr_ngc,:)
      loc_buffBNE(:,:,:,:) = rcv_buffUE(:,:,1:gr_ngc,:)
      if (sim_bcTypeY == "outflow") then
        do j = 1, gr_ngc
          rcv_buffBSW(:,:,j,:) = loc_buffBSW(:,:,1,:)
          rcv_buffBSE(:,:,j,:) = loc_buffBSE(:,:,1,:)
          rcv_buffBNW(:,:,j,:) = loc_buffBNW(:,:,1,:)
          rcv_buffBNE(:,:,j,:) = loc_buffBNE(:,:,1,:)
        end do
      else if (sim_bcTypeY == "reflect") then
        do j = 1, gr_ngc
          rcv_buffBSW(:,:,j,:) = loc_buffBSW(:,:,gr_ngc-j+1,:)
          rcv_buffBSE(:,:,j,:) = loc_buffBSE(:,:,gr_ngc-j+1,:)
          rcv_buffBNW(:,:,j,:) = loc_buffBNW(:,:,gr_ngc-j+1,:)
          rcv_buffBNE(:,:,j,:) = loc_buffBNE(:,:,gr_ngc-j+1,:)
        end do
        rcv_buffBSW(VELY_VAR,:,:,:) = -rcv_buffBSW(VELY_VAR,:,:,:)
        rcv_buffBSE(VELY_VAR,:,:,:) = -rcv_buffBSE(VELY_VAR,:,:,:)
        rcv_buffBNW(VELY_VAR,:,:,:) = -rcv_buffBNW(VELY_VAR,:,:,:)
        rcv_buffBNE(VELY_VAR,:,:,:) = -rcv_buffBNE(VELY_VAR,:,:,:)
      end if  ! sim_bcTypeY
    end if    ! bl_j

    if (bl_j == bl_jProcs) then   ! top boundary
      loc_buffTSW(:,:,:,:) = rcv_buffDW(:,:,gr_ny-gc:gr_ny,:)
      loc_buffTSE(:,:,:,:) = rcv_buffDE(:,:,gr_ny-gc:gr_ny,:)
      loc_buffTNW(:,:,:,:) = rcv_buffUW(:,:,gr_ny-gc:gr_ny,:)
      loc_buffTNE(:,:,:,:) = rcv_buffUE(:,:,gr_ny-gc:gr_ny,:)
      if (sim_bcTypeY == "outflow") then
        do j = 1, gr_ngc
          rcv_buffTSW(:,:,j,:) = loc_buffTSW(:,:,gr_ngc,:)
          rcv_buffTSE(:,:,j,:) = loc_buffTSE(:,:,gr_ngc,:)
          rcv_buffTNW(:,:,j,:) = loc_buffTNW(:,:,gr_ngc,:)
          rcv_buffTNE(:,:,j,:) = loc_buffTNE(:,:,gr_ngc,:)
        end do
      else if (sim_bcTypeY == "reflect") then
        do j = 1, gr_ngc
          rcv_buffTSW(:,:,j,:) = loc_buffTSW(:,:,gr_ngc-j+1,:)
          rcv_buffTSE(:,:,j,:) = loc_buffTSE(:,:,gr_ngc-j+1,:)
          rcv_buffTNW(:,:,j,:) = loc_buffTNW(:,:,gr_ngc-j+1,:)
          rcv_buffTNE(:,:,j,:) = loc_buffTNE(:,:,gr_ngc-j+1,:)
        end do
        rcv_buffTSW(VELY_VAR,:,:,:) = -rcv_buffTSW(VELY_VAR,:,:,:)
        rcv_buffTSE(VELY_VAR,:,:,:) = -rcv_buffTSE(VELY_VAR,:,:,:)
        rcv_buffTNW(VELY_VAR,:,:,:) = -rcv_buffTNW(VELY_VAR,:,:,:)
        rcv_buffTNE(VELY_VAR,:,:,:) = -rcv_buffTNE(VELY_VAR,:,:,:)
      end if  ! sim_bcTypeY
    end if    ! bl_j

    if (bl_k == 1) then   ! down boundary
      loc_buffBSW(:,:,:,:) = rcv_buffBW(:,:,:,1:gr_ngc)
      loc_buffBSE(:,:,:,:) = rcv_buffBE(:,:,:,1:gr_ngc)
      loc_buffTSW(:,:,:,:) = rcv_buffTW(:,:,:,1:gr_ngc)
      loc_buffTSE(:,:,:,:) = rcv_buffTE(:,:,:,1:gr_ngc)
      if (sim_bcTypeZ == "outflow") then
        do k = 1, gr_ngc
          rcv_buffBSW(:,:,:,k) = loc_buffBSW(:,:,:,1)
          rcv_buffBSE(:,:,:,k) = loc_buffBSE(:,:,:,1)
          rcv_buffTSW(:,:,:,k) = loc_buffTSW(:,:,:,1)
          rcv_buffTSE(:,:,:,k) = loc_buffTSE(:,:,:,1)
        end do
      else if (sim_bcTypeZ == "reflect") then
        do k = 1, gr_ngc
          rcv_buffBSW(:,:,:,k) = loc_buffBSW(:,:,:,gr_ngc-k+1)
          rcv_buffBSE(:,:,:,k) = loc_buffBSE(:,:,:,gr_ngc-k+1)
          rcv_buffTSW(:,:,:,k) = loc_buffTSW(:,:,:,gr_ngc-k+1)
          rcv_buffTSE(:,:,:,k) = loc_buffTSE(:,:,:,gr_ngc-k+1)
        end do
        rcv_buffBSW(VELZ_VAR,:,:,:) = -rcv_buffBSW(VELZ_VAR,:,:,:)
        rcv_buffBSE(VELZ_VAR,:,:,:) = -rcv_buffBSE(VELZ_VAR,:,:,:)
        rcv_buffTSW(VELZ_VAR,:,:,:) = -rcv_buffTSW(VELZ_VAR,:,:,:)
        rcv_buffTSE(VELZ_VAR,:,:,:) = -rcv_buffTSE(VELZ_VAR,:,:,:)
      end if  ! sim_bcTypeZ
    end if    ! bl_k

    if (bl_k == bl_kProcs) then   ! up boundary
      loc_buffBNW(:,:,:,:) = rcv_buffBW(:,:,:,gr_nz-gc:gr_nz)
      loc_buffBNE(:,:,:,:) = rcv_buffBE(:,:,:,gr_nz-gc:gr_nz)
      loc_buffTNW(:,:,:,:) = rcv_buffTW(:,:,:,gr_nz-gc:gr_nz)
      loc_buffTNE(:,:,:,:) = rcv_buffTE(:,:,:,gr_nz-gc:gr_nz)
      if (sim_bcTypeZ == "outflow") then
        do k = 1, gr_ngc
          rcv_buffBNW(:,:,:,k) = loc_buffBNW(:,:,:,gr_ngc)
          rcv_buffBNE(:,:,:,k) = loc_buffBNE(:,:,:,gr_ngc)
          rcv_buffTNW(:,:,:,k) = loc_buffTNW(:,:,:,gr_ngc)
          rcv_buffTNE(:,:,:,k) = loc_buffTNE(:,:,:,gr_ngc)
        end do
      else if (sim_bcTypeZ == "reflect") then
        do k = 1, gr_ngc
          rcv_buffBNW(:,:,:,k) = loc_buffBNW(:,:,:,gr_ngc-k+1)
          rcv_buffBNE(:,:,:,k) = loc_buffBNE(:,:,:,gr_ngc-k+1)
          rcv_buffTNW(:,:,:,k) = loc_buffTNW(:,:,:,gr_ngc-k+1)
          rcv_buffTNE(:,:,:,k) = loc_buffTNE(:,:,:,gr_ngc-k+1)
        end do
        rcv_buffBNW(VELZ_VAR,:,:,:) = -rcv_buffBNW(VELZ_VAR,:,:,:)
        rcv_buffBNE(VELZ_VAR,:,:,:) = -rcv_buffBNE(VELZ_VAR,:,:,:)
        rcv_buffTNW(VELZ_VAR,:,:,:) = -rcv_buffTNW(VELZ_VAR,:,:,:)
        rcv_buffTNE(VELZ_VAR,:,:,:) = -rcv_buffTNE(VELZ_VAR,:,:,:)
      end if  ! sim_bcTypeZ
    end if    ! bl_k
    !!! done: domain BCs for cubes
    !!! done: all domain BCs


    ! filling halo
    ! sides
    V(:,      i0:i0+gc,      j0:j0+gc, kbeg:kend) = rcv_buffBW   !3
    V(:, imax-gc:imax,       j0:j0+gc, kbeg:kend) = rcv_buffBE   !4
    V(:,      i0:i0+gc, jmax-gc:jmax,  kbeg:kend) = rcv_buffTW   !1
    V(:, imax-gc:imax,  jmax-gc:jmax,  kbeg:kend) = rcv_buffTE   !2

    V(:,      i0:i0+gc, jbeg:jend,      k0:k0+gc) = rcv_buffDW   !5
    V(:, imax-gc:imax,  jbeg:jend,      k0:k0+gc) = rcv_buffDE   !6
    V(:,      i0:i0+gc, jbeg:jend, kmax-gc:kmax ) = rcv_buffUW   !7
    V(:, imax-gc:imax,  jbeg:jend, kmax-gc:kmax ) = rcv_buffUE   !8

    V(:, ibeg:iend,      j0:j0+gc,      k0:k0+gc) = rcv_buffBS   !9
    V(:, ibeg:iend,      j0:j0+gc, kmax-gc:kmax ) = rcv_buffBN   !10
    V(:, ibeg:iend, jmax-gc:jmax,       k0:k0+gc) = rcv_buffTS   !11
    V(:, ibeg:iend, jmax-gc:jmax,  kmax-gc:kmax ) = rcv_buffTN   !12

    ! corners
    V(:,      i0:i0+gc,     j0:j0+gc,      k0:k0+gc) = rcv_buffBSW   !1
    V(:, imax-gc:imax,      j0:j0+gc,      k0:k0+gc) = rcv_buffBSE   !2
    V(:,      i0:i0+gc,     j0:j0+gc, kmax-gc:kmax ) = rcv_buffBNW   !3
    V(:, imax-gc:imax,      j0:j0+gc, kmax-gc:kmax ) = rcv_buffBNE   !4

    V(:,      i0:i0+gc, jmax-gc:jmax,      k0:k0+gc) = rcv_buffTSW   !5
    V(:, imax-gc:imax,  jmax-gc:jmax,      k0:k0+gc) = rcv_buffTSE   !6
    V(:,      i0:i0+gc, jmax-gc:jmax, kmax-gc:kmax ) = rcv_buffTNW   !7
    V(:, imax-gc:imax,  jmax-gc:jmax, kmax-gc:kmax ) = rcv_buffTNE   !8

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
  end subroutine bc_corner


end module bc
