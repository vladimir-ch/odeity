#ifndef TEST_SYSTEM_05_H
#define TEST_SYSTEM_05_H

#include "ExplicitOde.h"

// Stiff DETEST A3
// initial condition y(*) = 1.0
class TestSystem05 : public ExplicitOde {
public:
  TestSystem05( )
  : ExplicitOde( 4 )
  {}

  virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
  {
    ydot(0) = -1.0e4 * y(0) + 1.0e2 * y(1) - 1.0e1 * y(2) + y(3);
    ydot(1) = -1.0e3 * y(1) + 1.0e1 * y(2) - 1.0e1 * y(3);
    ydot(2) = -y(2) + 1.0e1 * y(3);
    ydot(3) = -1.0e-1 * y(3);

    return 0;
  }

};

#endif // TEST_SYSTEM_04_H


