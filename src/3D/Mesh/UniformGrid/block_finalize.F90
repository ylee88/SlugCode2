subroutine block_finalize()

  use block_data, only: bl_grid,            &
                        bl_grid_ext,        &
                        bl_BC,              &
                        bl_corner_cornerBC, &
                        bl_corner_sideBC

  implicit none

  integer :: ierr

  if (allocated(bl_grid) .eqv. .true.) deallocate(bl_grid)
  if (allocated(bl_grid_ext) .eqv. .true.) deallocate(bl_grid_ext)
  if (allocated(bl_BC) .eqv. .true.) deallocate(bl_BC)
  if (allocated(bl_corner_cornerBC) .eqv. .true.) deallocate(bl_corner_cornerBC)
  if (allocated(bl_corner_sideBC) .eqv. .true.) deallocate(bl_corner_sideBC)

  call MPI_Finalize(ierr)

  return
end subroutine block_finalize
