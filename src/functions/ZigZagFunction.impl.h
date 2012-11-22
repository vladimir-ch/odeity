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

#include "../utils/LogStream.h"

ZigZagFunction::ZigZagFunction( const unsigned int n_components )
    :
        Function<2>( n_components )
{}



realtype
ZigZagFunction::value( const Point<2>& p, const unsigned int component )
{
    const double K = 1.0 + sqrt(5.0)/2.0;
    double d0x;
    double y = p(1);
    double ax = fabs(p(0)/0.15);

    if ( ax < 1.0 )
        d0x = -K + (1+K)*ax;
    else if ( 1.0 <= ax && ax < 2.0)
        d0x = 2-ax;
    else
        d0x = 0.0;

    d0x = 1.0 + 0.3*d0x;
    
    return tanh((y-d0x)/0.015625);
}


void
ZigZagFunction::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Zig-zag function from Baensch, Morin, Nochetto" << endl;
    logger << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
}

