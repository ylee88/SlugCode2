subroutine abort_slug(msg)

#include "definition.h"

  use mpi

  implicit none

  character(len=MAX_STRING_LENGTH), intent(IN) :: msg
  integer :: ierr

  write(*,*)'ERROR: ', msg
  call MPI_Abort(MPI_COMM_WORLD, ierr, ierr)
  stop

end subroutine abort_slug
