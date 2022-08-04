module num_interface

#include "definition.h"

  implicit none

  abstract interface
    subroutine temporal(dt)
      implicit none

      real, intent(IN) :: dt

    end subroutine temporal
  end interface

  abstract interface
    subroutine spatial(dt, radius, stnclL, stnclR, reconL, reconR)
      implicit none

      real, intent(IN) :: dt
      integer, intent(IN) :: radius
      real, dimension(:, :), intent(IN) :: stnclL, stnclR
      real, dimension(NSYS_VAR) :: reconL, reconR
    end subroutine spatial
  end interface


end module num_interface
