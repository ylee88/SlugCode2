module num_interface

#include "definition.h"

  use grid_data, only: gr_imax

  implicit none

  abstract interface
    subroutine temporal(dt)
      implicit none

      real, intent(IN) :: dt

    end subroutine temporal
  end interface

  abstract interface
    subroutine spatial(dt, radius, stencilV, reconL, reconR, dir)
      implicit none

      real, intent(IN) :: dt
      integer, intent(IN) :: radius, dir
      real, dimension(NUMB_VAR, radius*2+1), intent(IN) :: stencilV
      real, dimension(NSYS_VAR) :: reconL, reconR
    end subroutine spatial
  end interface


end module num_interface
