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

#ifndef FUNCTION_H
#define FUNCTION_H

#include "FunctionTime.h"
#include "../geometry/Point.h"
#include "../utils/ParameterHandler.h"


template<int dim>
class Function : public FunctionTime
{
    public:

        static const unsigned int dimension = dim;
        static const double smallNumber = 1.0e-2;

        const unsigned int numberOfComponents;

        Function( const unsigned int n_components = 1, const realtype initial_time = 0.0 )
            :
                FunctionTime( initial_time ),
                numberOfComponents( n_components )
        {
            Assert( n_components > 0, ExcZero() );
        }

        virtual ~Function()
        {}

        Function & operator = ( const Function & f );
        virtual realtype value( const Point<dim>& p, const unsigned int component = 0 );
        realtype operator () ( const Point<dim>& p, const unsigned int component = 0 )
        {
            return value( p, component );
        }

        virtual void getParameters( ParameterHandler & prm ) {}
        virtual void printInfo() const {}

        DeclException2 (ExcNumberOfComponents,
                int, int,
                << "You can only assign function objects with the same "
                << "number of components. However, here the left operand "
                << "has " << arg1 << " components, and the right operand "
                << arg2 << " components.");
};

#endif // FUNCTION_H

