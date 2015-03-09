#ifndef TEST_SYSTEM_02_H
#define TEST_SYSTEM_02_H

#include "ExplicitOde.h"

// non-stiff DETEST B1
class TestSystem02 : public ExplicitOde {
public:
  TestSystem02( )
  : ExplicitOde( 2 )
  {}

  virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
  {
    ydot(0) = 2.0 * ( y(0) - y(0) * y(1) );
    ydot(1) = -( y(1) - y(0) * y(1) );

    return 0;
  }

};

#endif // TEST_SYSTEM_02_H
