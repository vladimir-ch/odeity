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

#ifndef RECTANGLE_FUNCTION_H
#define RECTANGLE_FUNCTION_H

#include "Function.h"


class RectangleFunction : public Function<2>
{
    public:
        RectangleFunction(
                const Point<2>& lowerLeftCorner = Point<2>( 0.2, 0.2 ),
                const Point<2>& upperRightCorner = Point<2>( 0.8, 0.8 ),
                const unsigned int n_components = 1 );

        realtype value( const Point<2>& p, const unsigned int component = 0 );
        void printInfo() const;

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

    private:
        Point<2> lowerLeftCorner_;
        Point<2> upperRightCorner_;
};

#include "RectangleFunction.impl.h"

Function<2> *
createRectangleFunction()
{
    return new RectangleFunction();
}

#endif // RECTANGLE_FUNCTION_H


