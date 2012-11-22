/* Copyright (c) 2010 Vladimir Chalupecky
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


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
