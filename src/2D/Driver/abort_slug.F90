subroutine abort_slug(msg)

#include "definition.h"

  use mpi

  implicit none

  character(len=*), intent(IN) :: msg
  integer :: errcode, ierr

  write(*,*)'ERROR: ', msg
  call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
  stop

end subroutine abort_slug
