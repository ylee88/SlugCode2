module grid_data
  implicit none

  integer, allocatable, dimension(:), save :: gr_i0, gr_ibeg, gr_iend, gr_imax
  real,    allocatable, dimension(:), save :: gr_xCoord
  real,    allocatable, dimension(:), save :: gr_xbeg, gr_xend

  integer, save                            :: gr_nx, gr_ngc, gr_glb_nx
  real,    save                            :: gr_dx

  !! dummies
  real,    save                            :: gr_dy, gr_dz

  real, allocatable, dimension(:,:), save   :: gr_maxSpeed  ! maximum speed
  real, allocatable, dimension(:,:), save   :: gr_U ! conservative vars
  real, allocatable, dimension(:,:), save   :: gr_V ! primitive vars

  real, allocatable, dimension(:,:,:), save :: gr_vL   ! left Riemann states
  real, allocatable, dimension(:,:,:), save :: gr_vR   ! right Riemann states

  real, allocatable, dimension(:,:,:), save :: gr_fL   ! left interface flux
  real, allocatable, dimension(:,:,:), save :: gr_fR   ! right interface flux
  real, allocatable, dimension(:,:,:), save :: gr_flux ! fluxes

end module grid_data
