module gp_interface

#include "definition.h"

  implicit none

  abstract interface
    function kernel(x, y, eldel) result(f)
      implicit none

      real(kind=16), intent(IN) :: x, y, eldel
      real(kind=16) :: f
    end function kernel
  end interface

end module gp_interface
