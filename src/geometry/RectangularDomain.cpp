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

#include "RectangularDomain.h"

template<>
void
RectangularDomain<1>::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Rectangular domain");
        minExtents_[xDim] = prm.get_double( "x-min" );
        maxExtents_[xDim] = prm.get_double( "x-max" );
    prm.leave_subsection();
}


template<>
void
RectangularDomain<2>::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Rectangular domain");
        minExtents_[xDim] = prm.get_double( "x-min" );
        maxExtents_[xDim] = prm.get_double( "x-max" );
        minExtents_[yDim] = prm.get_double( "y-min" );
        maxExtents_[yDim] = prm.get_double( "y-max" );
    prm.leave_subsection();
}


template<>
void
RectangularDomain<3>::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Rectangular domain");
        minExtents_[xDim] = prm.get_double( "x-min" );
        maxExtents_[xDim] = prm.get_double( "x-max" );
        minExtents_[yDim] = prm.get_double( "y-min" );
        maxExtents_[yDim] = prm.get_double( "y-max" );
        minExtents_[zDim] = prm.get_double( "z-min" );
        maxExtents_[zDim] = prm.get_double( "z-max" );
    prm.leave_subsection();
}

