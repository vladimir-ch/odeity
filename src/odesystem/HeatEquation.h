#ifndef HEAT_EQUATION_H
#define HEAT_EQUATION_H

#include "ExplicitOde.h"

#include <cmath>

class HeatEquation : public ExplicitOde
{
    public:

        HeatEquation( const int numberOfEquations, const realtype diffusionCoefficient )
            : ExplicitOde( numberOfEquations ),
            diffusionCoefficient_( diffusionCoefficient )
        {}

        virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
        {
            realtype u, ul, ur;
            int neq = y.size();
            realtype spatialStepsize = 1.0 / ( neq - 1 );

            for ( int i = 0; i < neq; i++)
            {
                u = y(i);
                ul = ( i == 0 ) ? y(1) : y(i-1);
                ur = ( i == neq-1 ) ? y(i-1) : y(i+1);

                ydot(i) = diffusionCoefficient_ * ( ul - 2*u + ur ) / spatialStepsize /spatialStepsize
                    - (1.0/diffusionCoefficient_) * u * ( u-1.0)*(u-0.5)
                    - std::fabs(ur - ul) / 2 / spatialStepsize * diffusionCoefficient_;
            }

            return 0;
        }

        realtype getDiffusionCoefficient() const
        {
            return diffusionCoefficient_;
        }

    private:
        realtype diffusionCoefficient_;
};

#endif // HEAT_EQUATION_H
