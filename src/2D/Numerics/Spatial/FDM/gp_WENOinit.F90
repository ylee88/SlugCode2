subroutine gp_WENOinit()

#include "definition.h"

  use linalg
  use gp_data
  use GP, only: qp, SE, SE_intg, SE_intgVec
  use grid_data, only: gr_dx, gr_dy
  use sim_data,  only: sim_gpRadii,  &
                       sim_gpEll,    &
                       sim_gpEldel,  &
                       sim_gpKernel, &
                       sim_gpSigma,  &
                       sim_gpSigdel

  implicit none

  real(qp), dimension(2*sim_gpRadii+1, 2*sim_gpRadii+1) :: C, L
  real(qp), dimension(2, 2*sim_gpRadii+1) :: T, Z
  real(qp), dimension(2*sim_gpRadii+1) :: stencil

  real(qp), dimension(sim_gpRadii+1, sim_gpRadii+1) :: Ck, Lk
  real(qp), dimension(2, sim_gpRadii+1, sim_gpRadii+1) :: Zk
  real(qp), dimension(2, sim_gpRadii+1) :: Tk, linW

  real(qp), dimension(sim_gpRadii+1, sim_gpRadii+1, sim_gpRadii+1) :: Zm

  ! used for eigensystem. need to be double precisions
  real, dimension(sim_gpRadii+1, sim_gpRadii+1) :: Cm
  real, dimension(sim_gpRadii+1) :: eigvals
  real, dimension(66*sim_gpRadii+66) :: WORK
  real, dimension(sim_gpRadii+1, sim_gpRadii+1) :: Pvecs

  real(qp), dimension(2*sim_gpRadii+2, sim_gpRadii+1) :: Zmat
  real(qp), dimension(2*sim_gpRadii+2) :: Zvec

  real(qp), dimension(sim_gpRadii+1) :: ul, Pk
  real(qp), dimension(2*sim_gpRadii+1) :: un

  real(qp) :: eldel, sigdel

  integer :: LR, i, j, N, M, ROW, COL, R, LDA, LWORK, INFO


  ! initialize vars
  R = sim_gpRadii
  N = 2*R+1
  C = 0._qp
  T = 0._qp
  L = 0._qp

  ! take care of \el
  if (sim_gpEll == 0. .and. sim_gpEldel == 0.) then
    call abort_slug("both \el and \el/\Delta are zero. exiting")
  elseif (sim_gpEldel == 0.) then
    sim_gpEldel = sim_gpEll/gr_dx
  elseif (sim_gpEll == 0.) then
    sim_gpEll = sim_gpEldel*gr_dx
  end if
  eldel = sim_gpEldel

  ! take care of \sigma
  if (sim_gpSigma == 0. .and. sim_gpSigdel == 0.) then
    call abort_slug("both \sigma and \sigma/\Delta are zero. exiting")
  elseif (sim_gpSigdel == 0.) then
    sim_gpSigdel = sim_gpSigma/gr_dx
  elseif (sim_gpSigma == 0.) then
    sim_gpSigma = sim_gpSigdel*gr_dx
  end if
  sigdel = sim_gpSigdel

  ! select kernel
  if (sim_gpKernel == 'SE') then
    gp_kernel      => SE
    gp_intgKernel  => SE_intg
    gp_predVec     => SE_intgVec
  else
    call abort_slug('Unrecognized sim_gpKernel')
  end if

  ! first we work on the global stencil (2*R+1 points)
  ! build stencils
  do i = 1, N
    stencil(i) = REAL(i-R-1, qp)
  end do

  !calculate covariance and prediction vector
  do i = 1, N
    do j = 1, N
      C(i,j) = gp_intgKernel(stencil(i), stencil(j), eldel)
    end do
    T(1, i) = gp_predVec(-0.5_qp, stencil(i), eldel)
    T(2, i) = gp_predVec( 0.5_qp, stencil(i), eldel)
  end do

  !solve for weights
  call chol(C, N, L)
  ! call solve_Axb(C, gp4_v, u, L, N)
  call solve_CZT(C, Z, T, L, N)


  !now lets work on the ENO stencils
  N = R+1
  Zm = 0._qp
  do m = 1, N
    !loop over the stencils
    !build k-th stencil
    do i = 1, N
      stencil(i) = REAL(i-1-R+m-1, qp)
    end do
    do i = 1, N
      do j = 1, N
        Ck(i,j) = gp_intgKernel(stencil(i), stencil(j), eldel)
      end do
      Tk(1, i) = gp_predVec(-0.5_qp, stencil(i), eldel)
      Tk(2, i) = gp_predVec( 0.5_qp, stencil(i), eldel)
    end do
    !solve for weights
    call chol(Ck, N, Lk)
    ! call solve_Axb(Ck, gp4_vk(:,m), uk, Lk, N)
    call solve_CZT(Ck, Zk(:,:,m), Tk, Lk, N)
    ! zm-vectors for eigen system
    do i = 1, N
      do j = 1, N
        Pk(j) = gp_predVec(stencil(j), stencil(i), eldel)
      end do
      call solve_Axb(Ck, Zm(:,i,m), Pk(:), Lk, N)    ! under the eqn 39
      ! Zm(:,:,m) are the same for all m's; checked via numpy
    end do

  end do



  !lastly we solve for linear GP-WENO weights
  ROW = 2*R+1
  COL = R+1
  ul = 1._qp
  un = 1._qp
  linW = 0._qp
  do LR = 1, 2
    Zmat = 0._qp
    do m = 1, COL
      Zmat(m:m+R,m) = Zk(LR,:,m)
      Zmat(2*R+2,m) = dot_product(Zk(LR,:,m), ul)
    end do

    Zvec(1:ROW) = Z(LR,:)
    Zvec(2*R+2) = dot_product(Z(LR,:), un)

    call LSTSQ(ROW+1, COL, Zmat, linW(LR,:), Zvec)      ! ean 34
  end do


  ! now we calculate P vectors for gp_betas
  N = R+1
  LDA = N
  LWORK = 66*N
  ! quad <- double precision
  ! sigdel = sim_gpSigdel
  ! sigdel = eldel
  do i = 1, N
     do j = 1, N
        ! quad  <-  double precision
        Cm(i,j) = gp_kernel(REAL(i, qp), REAL(j, qp), sigdel)    ! calc. K. See eqn 42 and 43
     end do
  end do
  call DSYEV('V', 'L', N, Cm, LDA, eigvals, WORK, LWORK, INFO)   ! compute eigvals/eigvecs of C |-> C, W
  if (INFO < 0) then
     call abort_slug("error in dsyev, info != 1")
  end if

  do i = 1, N
     do j = 1, N
        Pvecs(j,i) = dot_product(Cm(:,i), Zm(j,:,1))/sqrt(eigvals(i))
     end do
  end do


  !now we truncate to double precision
  gp_Zk    = Zk
  gp_linW  = linW
  gp_Pvecs = Pvecs

end subroutine gp_WENOinit

