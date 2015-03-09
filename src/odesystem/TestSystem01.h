#ifndef TEST_SYSTEM_01_H
#define TEST_SYSTEM_01_H

#include "ExplicitOde.h"

// fortex1.html
// initial condition y(0) = 0.0
class TestSystem01 : public ExplicitOde
{
    public:
        TestSystem01() : ExplicitOde( 1 ) {}

        virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
        {
            ydot(0) = -2.0 * t * y(0) + 4.0 * t;

            return 0;
        }

};

#endif // TEST_SYSTEM_01_H
