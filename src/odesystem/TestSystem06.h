#ifndef TEST_SYSTEM_06_H
#define TEST_SYSTEM_06_H

#include "ExplicitOde.h"

// Matlab ODE Suite - chm6ode
// initial condition y(*) = 0.0
class TestSystem06 : public ExplicitOde {
public:
  TestSystem06( )
  : ExplicitOde( 4 )
  {}

  virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
  {
    realtype K = std::pow(M_E, 20.7 - 1500.0/y(0) );
    ydot(0) = 1.3 * (y(2)-y(0)) + 10400.0 * K * y(1);
    ydot(1) = 1880.0 * (y(3) - y(1)*(1+K));
    ydot(2) = 1752.0 - 269.0 * y(2) + 267.0 * y(0);
    ydot(3) = 0.1 + 320.0 * y(1) - 321.0 * y(3);

    return 0;
  }

};

#endif // TEST_SYSTEM_06_H



