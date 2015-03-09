#ifndef TEST_SYSTEM_04_H
#define TEST_SYSTEM_04_H

#include "ExplicitOde.h"

// fortex2.html
// initial condition y(*) = 0.0
// exact solution: y1(t) = e^t + e^-t - 2, y2(t) = e^t - e^-t - t
class TestSystem04 : public ExplicitOde {
public:
  TestSystem04( )
  : ExplicitOde( 2 )
  {}

  virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
  {
    ydot(0) = y(1) + t;
    ydot(1) = y(0) + 1;

    return 0;
  }

};

#endif // TEST_SYSTEM_04_H


