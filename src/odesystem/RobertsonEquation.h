#ifndef ROBERTSON_EQUATION_H
#define ROBERTSON_EQUATION_H

#include "ExplicitOde.h"

class RobertsonEquation : public ExplicitOde
{
    public:
        RobertsonEquation() : ExplicitOde( 3 ) {}

        virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
        {
            ydot(0) = -0.04e0 * y(0) + 1.0e4 * y(1) * y(2);
            ydot(1) = -y(0) - y(2);
            ydot(2) = 3.0e7 * y(1) * y(1);

            return 0;
        }

};

#endif // ROBERTSON_EQUATION_H

