module sim_data

#include "definition.h"

  implicit none

  !! numerics
  real, save :: sim_cfl, sim_tmax, sim_outputIntervalTime, sim_dt, sim_WENeps
  integer, save :: sim_order, sim_nStep, sim_Torder, sim_mval
  character(len=MAX_STRING_LENGTH), save :: sim_name, sim_riemann, sim_WENO
  character(len=30), save :: sim_spatialMethod, sim_temporalMethod
  logical, save :: sim_charLimiting, sim_RK, sim_fixDt, sim_nlim

  !! ICs
  real, save :: sim_gamma
  real, save :: sim_smallPres
  character(len=MAX_STRING_LENGTH), save :: sim_icType

  !! BCs
  character(len=MAX_STRING_LENGTH), save :: sim_bcTypeX
  integer                         , save :: sim_xBC
  logical                         , save :: sim_cornerBC

  !! IO
  integer, save :: sim_ioNfreq
  real,    save :: sim_ioTfreq
  logical, save :: sim_hdf5, sim_pIO

  !! GP
  logical, save :: sim_gpWENO
  character(len=MAX_STRING_LENGTH), save :: sim_gpKernel
  integer, save :: sim_gpRadii
  real,    save :: sim_gpEll, sim_gpEldel, sim_gpSigma, sim_gpSigdel

  !! this is only for DMR BC
  real :: sim_t

end module sim_data
