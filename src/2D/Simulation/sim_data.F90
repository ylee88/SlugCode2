module sim_data

#include "definition.h"

  implicit none

  !! numerics
  real, save :: sim_cfl, sim_tmax, sim_outputIntervalTime, sim_dt, sim_WENeps
  integer, save :: sim_order, sim_nStep, sim_Torder, sim_mval
  character(len=MAX_STRING_LENGTH), save :: sim_name, sim_riemann, sim_WENO
  logical, save :: sim_charLimiting, sim_RK, sim_fixDt, sim_nlim, sim_reconMultiD

  !! ICs
  real, save :: sim_gamma
  real, save :: sim_smallPres
  character(len=MAX_STRING_LENGTH), save :: sim_icType

  !! BCs
  character(len=MAX_STRING_LENGTH), save :: sim_bcTypeX, sim_bcTypeY
  integer                         , save :: sim_xBC, sim_yBC

  !! IO
  integer, save :: sim_ioNfreq
  real,    save :: sim_ioTfreq
  logical, save :: sim_hdf5, sim_pIO

  !! GP
  ! character(len=MAX_STRING_LENGTH), save :: sim_quad, sim_gp_kernel
  ! real, save :: sim_sigdel, sim_sigma, sim_matern_nu, sim_RQ_alpha

end module sim_data
