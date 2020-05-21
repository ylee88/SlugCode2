module GP

#include "definition.h"

  use gp_data

  real :: rt2 = SQRT(2.)
  integer, parameter :: qp=kind(0.q0)           ! quad precision

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Kernel Functions !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!! Squared Exponential!!!!!!!!!!!!!!!!!!
  function SE(x, y, eldel) result(f)

    implicit none

    real(qp), intent(IN) :: x, y, eldel
    real(qp) :: f, r

    r = abs(x-y)
    f = EXP( -0.5*(r/eldel)**2 )

    return
  end function SE

  function SE_intg(x1, x2, eldel) result(f)

    !exact quadrature, only good for SE kernel
    real(qp), intent(IN) :: x1, x2, eldel

    real(qp) :: f, yxp, yxn, yxm

    yxp = (x1 - x2 + 1.)/(rt2*eldel)
    yxn = (x1      -x2)/(rt2*eldel)
    yxm = (x1 - x2 -1.)/(rt2*eldel)


    f = SQRT(PI)*(eldel)**2 *( yxp*ERF(yxp) + yxm*ERF(yxm) &
      - 2.*( yxn*ERF(yxn) + 1./SQRT(PI) *EXP(-yxn**2) ) &
      + 1./SQRT(PI) * ( EXP(-yxp**2) + exp(-yxm**2) ) )
    return
  end function SE_intg

  function SE_intgVec(x, t, eldel) result(f)
    implicit none
    real(qp), intent(IN) :: x, t, eldel
    real(qp) :: f

    f = eldel*SQRT(.5*PI)*(ERF( (x + .5 - t)/(rt2*eldel)) - ERF( (x - .5 - t)/(rt2*eldel) ))

  end function SE_intgVec




end module GP
