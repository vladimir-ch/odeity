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

#ifndef CAHN_HILLIARD_EQUATION_H
#define CAHN_HILLIARD_EQUATION_H

#include "MolOdeSystem.h"
#include "../utils/Vector.h"

class ParameterHandler;

class CahnHilliardEquation : public MolOdeSystem<2>
{
  public:

    CahnHilliardEquation();

    int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );
    int jacobian( const Vector<realtype>& v, Vector<realtype>& Jv,
        const realtype t, const Vector<realtype>& y, const Vector<realtype>& fy,
        const Vector<realtype>& tmp );

    std::string name() const;
    std::string componentName( int i ) const { return "phase_field"; }
    realtype interfaceWidth() const;
    realtype potentialCoefficient() const;

    void printInfo() const;
    static void declareParameters( ParameterHandler& prm );
    void getParameters( ParameterHandler& prm );

    void attachGrid( const RectangularGrid<2>& grid );

  private:

    realtype xi_, xiInv_, xiSqr_;
    realtype a_;
    Vector<realtype> w_;

    realtype f0( const realtype u );
    realtype f0deriv( const realtype u);

};


inline
std::string
CahnHilliardEquation::name() const
{
  return "Cahn-Hilliard equation with constant mobility";
}


inline
realtype
CahnHilliardEquation::interfaceWidth() const
{
  return xi_;
}


inline
realtype
CahnHilliardEquation::potentialCoefficient() const
{
    return a_;
}


#endif // CAHN_HILLIARD_EQUATION_H

