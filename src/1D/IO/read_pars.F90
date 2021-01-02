subroutine read_pars(filename)

#include "definition.h"

  use grid_data
  use sim_data
  use block_data
  use read_initFile

  implicit none

  character(len=*), intent(IN) :: filename

  !read sim_data
  call read_initFileChar (filename,'sim_spatialMethod', sim_spatialMethod)
  call read_initFileChar (filename,'sim_temporalMethod', sim_temporalMethod)

  call read_initFileInt (filename,'sim_nstep',    sim_nStep)

  call read_initFileReal(filename,'sim_dt',       sim_dt)
  call read_initFileReal(filename,'sim_cfl',      sim_cfl)
  call read_initFileReal(filename,'sim_tmax',     sim_tmax)
  call read_initFileReal(filename,'sim_WENeps',   sim_WENeps)
  call read_initFileReal(filename,'sim_outputIntervalTime',sim_outputIntervalTime)

  call read_initFileChar(filename,'sim_name',    sim_name)
  call read_initFileChar(filename,'sim_riemann', sim_riemann)
  call read_initFileChar(filename,'sim_WENO',    sim_WENO)

  call read_initFileBool(filename,'sim_charLimiting', sim_charLimiting)
  call read_initFileBool(filename,'sim_fixDt', sim_fixDt)
  call read_initFileBool(filename,'sim_nlim',  sim_nlim)
  call read_initFileBool(filename,'sim_hdf5',  sim_hdf5)
  call read_initFileBool(filename,'sim_pIO',   sim_pIO)

  call read_initFileBool(filename,'sim_cornerBC', sim_cornerBC)

  call read_initFileChar(filename,'sim_icType',   sim_icType)
  call read_initFileReal(filename,'sim_gamma',    sim_gamma)
  call read_initFileReal(filename,'sim_smallPres', sim_smallPres)

  call read_initFileChar(filename,'sim_bcTypeX', sim_bcTypeX)

  call read_initFileReal(filename,'sim_ioTfreq',  sim_ioTfreq)
  call read_initFileInt (filename,'sim_ioNfreq',  sim_ioNfreq)
  call read_initFileInt (filename,'sim_mval',     sim_mval)


  !read grid data
  allocate(gr_xend(NDIM)); allocate(gr_xbeg(NDIM))
  call read_initFileInt (filename,'gr_nx',   gr_nx)
  call read_initFileInt (filename,'gr_ngc',  gr_ngc)

  call read_initFileReal(filename,'gr_xbeg', gr_xbeg(XDIM))
  call read_initFileReal(filename,'gr_xend', gr_xend(XDIM))

  !read block data
  call read_initFileInt (filename,'bl_iProcs',   bl_iProcs)


  ! !read in GP pars
  call read_initFileBool(filename,'sim_gpWENO',   sim_gpWENO)
  call read_initFileChar(filename,'sim_gpKernel', sim_gpKernel)
  call read_initFileInt (filename,'sim_gpRadii',  sim_gpRadii)
  call read_initFileReal(filename,'sim_gpEll',     sim_gpEll)
  call read_initFileReal(filename,'sim_gpEldel',  sim_gpEldel)
  call read_initFileReal(filename,'sim_gpSigdel',  sim_gpSigdel)

end subroutine read_pars
