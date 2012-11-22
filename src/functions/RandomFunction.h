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

#ifndef RANDOM_FUNCTION_H
#define RANDOM_FUNCTION_H

#include "Function.h"


template<int dim>
class RandomFunction : public Function<dim>
{
    public:

        RandomFunction( const realtype mean = 0.0, const realtype sigma = 1.0, const unsigned int n_components = 1 );
        RandomFunction( const RandomFunction & other )
            :
                mean_( other.mean_ ), sigma_( other.sigma_ ), valid_( false )
        {}

        realtype value ( const Point<dim>& p, const unsigned int component = 0 );

        /**
         * Returns random double between 0 and 1.
         * Based on Wichmann & Hill Applied Stats AS183
         * 1982 Vol 31 188-190.
         */
        realtype uniform();

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

        void printInfo() const;

    private:
        realtype mean_, sigma_;
        realtype r1_, r2_, cachedRho_;
        bool valid_;

        int x,y,z; /// variables for uniform random number generator
};

#include "RandomFunction.impl.h"

template<int dim>
Function<dim> *
createRandomFunction()
{
    return new RandomFunction<dim>();
}

#endif // RANDOM_FUNCTION_H

