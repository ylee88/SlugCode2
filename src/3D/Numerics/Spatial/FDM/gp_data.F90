module gp_data

#include "definition.h"

  use gp_interface, only: kernel

  implicit none

  ! GP data
  real, allocatable, dimension(:,:), save :: gp_Pvecs

  ! these will be truncated from quadruple precision
  real, allocatable, dimension(:,:), save :: gp_linW
  real, allocatable, dimension(:,:,:), save :: gp_Zk

  procedure(kernel), pointer, save :: gp_kernel
  procedure(kernel), pointer, save :: gp_intgKernel
  procedure(kernel), pointer, save :: gp_predVec

  ! real(KIND=16), allocatable, dimension(:) :: gp4_v
  ! real(KIND=16), allocatable, dimension(:,:) :: gp4_w, gp4_z , gp4_vk
  ! real(KIND=16), allocatable, dimension(:,:,:) :: gp4_zk, gp4_Zvecs

  ! for MultiD
  ! integer, save :: gp_Npts
  ! real,    save :: gp_Xdel, gp_Ydel
  !
  ! real(KIND=8), allocatable, dimension(:) :: gpM_v
  ! real(KIND=8), allocatable, dimension(:,:) :: gpM_w, gpM_z
  !
  ! real(KIND=16), allocatable, dimension(:) :: gpM4_v
  ! real(KIND=16), allocatable, dimension(:,:) :: gpM4_w, gpM4_z


end module gp_data
