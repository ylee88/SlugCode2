subroutine grid_finalize()

  use grid_data, only: gr_xCoord, &
                       gr_i0,     &
                       gr_imax,   &
                       gr_ibeg,   &
                       gr_iend,   &
                       gr_U,      &
                       gr_V,      &
                       gr_vL,     &
                       gr_vR,     &
                       gr_fL,     &
                       gr_fR,     &
                       gr_flux,   &
                       gr_maxSpeed

  implicit none

  if (allocated(gr_xCoord) .eqv. .true.) deallocate(gr_xCoord)
  if (allocated(gr_i0)     .eqv. .true.) deallocate(gr_i0)
  if (allocated(gr_imax)   .eqv. .true.) deallocate(gr_imax)
  if (allocated(gr_ibeg)   .eqv. .true.) deallocate(gr_ibeg)
  if (allocated(gr_iend)   .eqv. .true.) deallocate(gr_iend)

  if (allocated(gr_U)    .eqv. .true.) deallocate(gr_U)
  if (allocated(gr_V)    .eqv. .true.) deallocate(gr_V)

  if (allocated(gr_vL)   .eqv. .true.) deallocate(gr_vL)
  if (allocated(gr_vR)   .eqv. .true.) deallocate(gr_vR)

  if (allocated(gr_fL)   .eqv. .true.) deallocate(gr_fL)
  if (allocated(gr_fR)   .eqv. .true.) deallocate(gr_fR)
  if (allocated(gr_flux) .eqv. .true.) deallocate(gr_flux)

  if (allocated(gr_maxSpeed) .eqv. .true.) deallocate(gr_maxSpeed)


  return
end subroutine grid_finalize
