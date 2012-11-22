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

#ifndef ALLEN_CAHN_EQUATION_H
#define ALLEN_CAHN_EQUATION_H

#include "MolOdeSystem.h"

class ParameterHandler;

class AllenCahnEquation : public MolOdeSystem<2>
{
    public:

        AllenCahnEquation();

        virtual int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );

        realtype interfaceWidth() const;
        realtype forcingTerm() const;

        std::string name() const;
        void printInfo() const;
        static void declareParameters( ParameterHandler& prm );
        void getParameters( ParameterHandler& prm );

        std::string componentName( int i ) const { return std::string("phase_field"); }

    private:

        realtype xi_, xiSqrInv_;
        realtype F_;
};


inline
std::string
AllenCahnEquation::name() const
{
    return "Allen-Cahn equation";
}


inline
realtype
AllenCahnEquation::interfaceWidth() const
{
    return xi_;
}



inline
realtype
AllenCahnEquation::forcingTerm() const
{
    return F_;
}



#endif // ALLEN_CAHN_EQUATION_H

