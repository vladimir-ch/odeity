#ifndef TEST_SYSTEM_03_H
#define TEST_SYSTEM_03_H

#include "ExplicitOde.h"

// Behind and Beyond the MATLAB ODE Suite
// van der pool oscillator
// initial conditions: y(0) = 0.0, y(1) = 0.25
class TestSystem03 : public ExplicitOde {
public:
  TestSystem03( )
  : ExplicitOde( 2 )
  {}

  virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
  {
    ydot(0) = y(1);
    ydot(1) = y(1)*(1.0-y(0)*y(0)) - y(0);

    return 0;
  }

};

#endif // TEST_SYSTEM_03_H

