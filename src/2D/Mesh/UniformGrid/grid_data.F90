module grid_data
  implicit none

  integer, allocatable, dimension(:), save :: gr_i0, gr_ibeg, gr_iend, gr_imax
  real,    allocatable, dimension(:), save :: gr_xCoord, gr_yCoord
  real,    allocatable, dimension(:), save :: gr_xbeg, gr_xend

  integer, save                            :: gr_nx , gr_ny, gr_ngc, gr_glb_nx, gr_glb_ny
  real,    save                            :: gr_dx, gr_dy

  real, allocatable, dimension(:,:),   save   :: gr_maxSpeed  ! maximum speed
  real, allocatable, dimension(:,:,:), save   :: gr_U ! conservative vars
  real, allocatable, dimension(:,:,:), save   :: gr_V ! primitive vars

  real, allocatable, dimension(:,:,:,:), save :: gr_vL   ! left Riemann states
  real, allocatable, dimension(:,:,:,:), save :: gr_vR   ! right Riemann states

  real, allocatable, dimension(:,:,:,:), save :: gr_fL   ! left interface flux
  real, allocatable, dimension(:,:,:,:), save :: gr_fR   ! right interface flux
  real, allocatable, dimension(:,:,:,:), save :: gr_flux ! fluxes

  !GP vars
  ! real, save                               :: gr_radius
  ! integer, save                            :: gr_gp_stencilPts, gr_Tcells
  ! real, allocatable, dimension(:  ), save  :: gr_GPv
  ! real, allocatable, dimension(:,:), save  :: gr_GPZ
  ! real, allocatable, dimension(:, :), save :: gr_GP_stencil

end module grid_data
