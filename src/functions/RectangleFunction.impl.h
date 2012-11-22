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


RectangleFunction::RectangleFunction(
        const Point<2>& lowerLeftCorner,
        const Point<2>& upperRightCorner,
        const unsigned int n_components )
:
    Function<2>( n_components ),
    lowerLeftCorner_( lowerLeftCorner ),
    upperRightCorner_( upperRightCorner )
{}


    
realtype
RectangleFunction::value( const Point<2>& p, const unsigned int component )
{
    realtype x = p(0);
    realtype y = p(1);

    if ( lowerLeftCorner_(0) <= x && x <= upperRightCorner_(0) && lowerLeftCorner_(1) <= y && y <= upperRightCorner_(1) )
        return 1.0; //-smallNumber;
    else
        return -1.0; //+smallNumber;

    // return std::min(
            // tanh((x-lowerLeftCorner_(0))/0.03125)+
            // tanh((upperRightCorner_(0)-x)/0.03125),
            // tanh((y-lowerLeftCorner_(1))/0.03125)+
            // tanh((upperRightCorner_(1)-y)/0.03125))-1.0;
}



void
RectangleFunction::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Rectangle function" << endl;
    logger << "~~~~~~~~~~~~~~~~~~" << endl;
    logger << "*Lower left corner*:: [" << lowerLeftCorner_( 0 ) << "," << lowerLeftCorner_( 1 ) << "]" << endl;
    logger << "*Upper right corner*:: [" << upperRightCorner_( 0 ) << "," << upperRightCorner_( 1 ) << "]" << endl;
}


void
RectangleFunction::declareParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Rectangle" );
        prm.declare_entry( "lower left x", "0.2", Patterns::Double() );
        prm.declare_entry( "lower left y", "0.2", Patterns::Double() );
        prm.declare_entry( "upper right x", "0.8", Patterns::Double() );
        prm.declare_entry( "upper right y", "0.8", Patterns::Double() );
    prm.leave_subsection();
}



void RectangleFunction::getParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Rectangle" );
        lowerLeftCorner_(0) = prm.get_double( "lower left x" );
        lowerLeftCorner_(1) = prm.get_double( "lower left y" );
        upperRightCorner_(0) = prm.get_double( "upper right x" );
        upperRightCorner_(1) = prm.get_double( "upper right y" );
    prm.leave_subsection();
}

