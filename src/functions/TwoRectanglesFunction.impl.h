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


TwoRectanglesFunction::TwoRectanglesFunction(
        const Point<2>& lowerLeftCorner1,
        const Point<2>& upperRightCorner1,
        const Point<2>& lowerLeftCorner2,
        const Point<2>& upperRightCorner2,
        const unsigned int n_components )
:
    Function<2>( n_components ),
    lowerLeftCorner1_( lowerLeftCorner1 ),
    upperRightCorner1_( upperRightCorner1 ),
    lowerLeftCorner2_( lowerLeftCorner2 ),
    upperRightCorner2_( upperRightCorner2 )
{}


    
realtype
TwoRectanglesFunction::value( const Point<2>& p, const unsigned int component )
{
    realtype x = p(0);
    realtype y = p(1);

    // if ( (lowerLeftCorner1_(0) <= x && x <= upperRightCorner1_(0) && lowerLeftCorner1_(1) <= y && y <= upperRightCorner1_(1)) ||
         // (lowerLeftCorner2_(0) <= x && x <= upperRightCorner2_(0) && lowerLeftCorner2_(1) <= y && y <= upperRightCorner2_(1))
       // )
        // return 1.0; //-smallNumber;
    // else
        // return -1.0; //+smallNumber;

    return
        std::max(
        std::min(
            tanh((x-lowerLeftCorner1_(0))/0.015625)+
            tanh((upperRightCorner1_(0)-x)/0.015625),
            tanh((y-lowerLeftCorner1_(1))/0.015625)+
            tanh((upperRightCorner1_(1)-y)/0.015625))-1,
        std::min(
            tanh((x-lowerLeftCorner2_(0))/0.015625)+
            tanh((upperRightCorner2_(0)-x)/0.015625),
            tanh((y-lowerLeftCorner2_(1))/0.015625)+
            tanh((upperRightCorner2_(1)-y)/0.015625))-1);
}



void
TwoRectanglesFunction::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Two rectangles function" << endl;
    logger << "~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    logger << "*Lower left corner 1*:: [" << lowerLeftCorner1_( 0 ) << "," << lowerLeftCorner1_( 1 ) << "]" << endl;
    logger << "*Upper right corner 1*:: [" << upperRightCorner1_( 0 ) << "," << upperRightCorner1_( 1 ) << "]" << endl;
    logger << "*Lower left corner 2*:: [" << lowerLeftCorner2_( 0 ) << "," << lowerLeftCorner2_( 1 ) << "]" << endl;
    logger << "*Upper right corner 2*:: [" << upperRightCorner2_( 0 ) << "," << upperRightCorner2_( 1 ) << "]" << endl;
}


void
TwoRectanglesFunction::declareParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Two rectangles" );
        prm.declare_entry( "center x 1", "0.3", Patterns::Double() );
        prm.declare_entry( "center y 1", "0.5", Patterns::Double() );
        prm.declare_entry( "width 1", "0.2", Patterns::Double() );
        prm.declare_entry( "height 1", "0.8", Patterns::Double() );
        prm.declare_entry( "center x 2", "0.7", Patterns::Double() );
        prm.declare_entry( "center y 2", "0.5", Patterns::Double() );
        prm.declare_entry( "width 2", "0.2", Patterns::Double() );
        prm.declare_entry( "height 2", "0.8", Patterns::Double() );
    prm.leave_subsection();
}



void TwoRectanglesFunction::getParameters( ParameterHandler & prm )
{
    Point<2> c1, c2;
    double width1, height1;
    double width2, height2;

    prm.enter_subsection( "Two rectangles" );
        c1(0) = prm.get_double( "center x 1" );
        c1(1) = prm.get_double( "center y 1" );
        width1 = prm.get_double( "width 1" );
        height1 = prm.get_double( "height 1" );
        c2(0) = prm.get_double( "center x 2" );
        c2(1) = prm.get_double( "center y 2" );
        width2 = prm.get_double( "width 2" );
        height2 = prm.get_double( "height 2" );
    prm.leave_subsection();

    lowerLeftCorner1_(0) = c1(0) - width1/2.0;
    lowerLeftCorner1_(1) = c1(1) - height1/2.0;
    upperRightCorner1_(0) = c1(0) + width1/2.0;
    upperRightCorner1_(1) = c1(1) + height1/2.0;
    lowerLeftCorner2_(0) = c2(0) - width2/2.0;
    lowerLeftCorner2_(1) = c2(1) - height2/2.0;
    upperRightCorner2_(0) = c2(0) + width2/2.0;
    upperRightCorner2_(1) = c2(1) + height2/2.0;
}

