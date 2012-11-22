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

#ifndef TWO_SCALE_SINE_FUNCTION_H
#define TWO_SCALE_SINE_FUNCTION_H

#include "Function.h"


class TwoScaleSineFunction : public Function<2>
{
    public:
        TwoScaleSineFunction(
                const realtype a = 0.1,
                const realtype b = 0.3,
                const realtype c = 0.5,
                const int m = 1,
                const int n = 16,
                const bool smooth = true,
                const realtype smoothXi = 0.0625,
                const unsigned int n_components = 1 );

        realtype value( const Point<2>& p, const unsigned int component = 0 );
        void printInfo() const;

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

    private:
        realtype a_;
        realtype b_;
        realtype c_;
        int m_;
        int n_;
        bool smooth_;
        realtype smoothXi_;
};

#include "TwoScaleSineFunction.impl.h"

Function<2> *
createTwoScaleSineFunction()
{
    return new TwoScaleSineFunction();
}

#endif // TWO_SCALE_SINE_FUNCTION_H


