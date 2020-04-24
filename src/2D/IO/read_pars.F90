subroutine read_pars(filename)

#include "definition.h"

  use grid_data
  use sim_data
  use block_data
  use read_initFile

  implicit none

  character(len=*), intent(IN) :: filename

  !read sim_data
  call read_initFileInt (filename,'sim_order',    sim_order)
  call read_initFileInt (filename,'sim_Torder',   sim_Torder)
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
  call read_initFileBool(filename,'sim_RK',    sim_RK)
  call read_initFileBool(filename,'sim_fixDt', sim_fixDt)
  call read_initFileBool(filename,'sim_nlim',  sim_nlim)
  call read_initFileBool(filename,'sim_hdf5',  sim_hdf5)
  call read_initFileBool(filename,'sim_pIO',   sim_pIO)

  call read_initFileBool(filename,'sim_cornerBC', sim_cornerBC)

  call read_initFileChar(filename,'sim_icType',   sim_icType)
  call read_initFileReal(filename,'sim_gamma',    sim_gamma)
  call read_initFileReal(filename,'sim_smallPres', sim_smallPres)

  call read_initFileChar(filename,'sim_bcTypeX', sim_bcTypeX)
  call read_initFileChar(filename,'sim_bcTypeY', sim_bcTypeY)

  call read_initFileReal(filename,'sim_ioTfreq',  sim_ioTfreq)
  call read_initFileInt (filename,'sim_ioNfreq',  sim_ioNfreq)
  call read_initFileInt (filename,'sim_mval',     sim_mval)


  !read grid data
  allocate(gr_xend(NDIM)); allocate(gr_xbeg(NDIM))
  call read_initFileInt (filename,'gr_nx',   gr_nx)
  call read_initFileInt (filename,'gr_ny',   gr_ny)
  call read_initFileInt (filename,'gr_ngc',  gr_ngc)

  call read_initFileReal(filename,'gr_xbeg', gr_xbeg(XDIM))
  call read_initFileReal(filename,'gr_xend', gr_xend(XDIM))
  call read_initFileReal(filename,'gr_ybeg', gr_xbeg(YDIM))
  call read_initFileReal(filename,'gr_yend', gr_xend(YDIM))

  !read block data
  call read_initFileInt (filename,'bl_iProcs',   bl_iProcs)
  call read_initFileInt (filename,'bl_jProcs',   bl_jProcs)


  ! !read in GP pars
  ! gpM_radius = 0.0
  ! call read_initFileReal(filename,'gr_radius',   gr_radius)
  ! call read_initFileChar(filename,'gp_quad'  ,  gp_quad  )
  ! call read_initFileChar(filename,'gp_kernel',  gp_kernel)
  ! call read_initFileReal(filename,'gp_ell'    ,  gp_el    )
  ! call read_initFileReal(filename,'gp_eldel' ,  gp_eldel )
  ! call read_initFileInt (filename,'gp_radius',  gp_radius)
  ! call read_initFileReal(filename,'gpM_radius'    ,  gpM_radius)

end subroutine read_pars
