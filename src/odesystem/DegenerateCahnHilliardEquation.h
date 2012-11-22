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


#ifndef DEG_CAHN_HILLIARD_EQUATION_H
#define DEG_CAHN_HILLIARD_EQUATION_H

#include "MolOdeSystem.h"
#include "../utils/Vector.h"

class ParameterHandler;

class DegenerateCahnHilliardEquation : public MolOdeSystem<2>
{
    public:

        DegenerateCahnHilliardEquation();

        int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );

        std::string name() const;
        std::string componentName( int i ) const { return "phase_field"; }

        realtype interfaceWidth() const;
        realtype alpha() const;
        realtype mobilityCoefficient() const;

        void printInfo() const;
        static void declareParameters( ParameterHandler& prm );
        void getParameters( ParameterHandler& prm );
        void attachGrid( const RectangularGrid<2>& grid );

    private:

        realtype xi_, xiPow2Inv_, xiPow2_;
        realtype xiPowAlpha_;
        realtype alpha_, beta_;
        realtype eps_;
        Vector<realtype> w_;
        Vector<realtype> mAvgX_;
        Vector<realtype> mAvgY_;

        realtype f0( const realtype u );
        realtype M( const realtype u );
};



inline
std::string
DegenerateCahnHilliardEquation::name() const
{
    return "Cahn-Hilliard equation with degenerate mobility";
}



inline
realtype
DegenerateCahnHilliardEquation::interfaceWidth() const
{
    return xi_;
}



inline
realtype
DegenerateCahnHilliardEquation::alpha() const
{
    return alpha_;
}



inline
realtype
DegenerateCahnHilliardEquation::mobilityCoefficient() const
{
    return beta_;
}



#endif // DEG_CAHN_HILLIARD_EQUATION_H


