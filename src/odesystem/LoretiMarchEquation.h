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

#ifndef LORETI_MARCH_EQUATION_H
#define LORETI_MARCH_EQUATION_H

#include "MolOdeSystem.h"
#include "../utils/Vector.h"

class ParameterHandler;

class LoretiMarchEquation : public MolOdeSystem<2>
{
    public:

        LoretiMarchEquation();

        int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );

        std::string name() const;
        std::string componentName( int i ) const { return "phase_field"; }
        realtype interfaceWidth() const;
        realtype potentialCoefficient() const;

        void printInfo() const;
        static void declareParameters( ParameterHandler& prm );
        void getParameters( ParameterHandler& prm );

        void attachGrid( const RectangularGrid<2>& grid );

    private:

        realtype xi_, xiInv_, xiSqr_, twoOverXi_, twoOverXiPow3_;
        realtype a_;
        Vector<realtype> w_;

        realtype PsiDer( const realtype u);
        realtype PsiDerDer( const realtype u);

};


inline
std::string
LoretiMarchEquation::name() const
{
  return "Loreti-March equation for motion V=H-2H_{ss}-H^3";
}


inline
realtype
LoretiMarchEquation::interfaceWidth() const
{
  return xi_;
}


inline
realtype
LoretiMarchEquation::potentialCoefficient() const
{
    return a_;
}


#endif // LORETI_MARCH_EQUATION_H

