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


TwoScaleSineFunction::TwoScaleSineFunction(
        const realtype a,
        const realtype b,
        const realtype c,
        const int m, 
        const int n, 
        const bool smooth,
        const realtype smoothXi,
        const unsigned int n_components )
:
    a_(a),
    b_(b),
    c_(c),
    m_(m),
    n_(n),
    smooth_( smooth ),
    smoothXi_( smoothXi )
{}


    
realtype
TwoScaleSineFunction::value( const Point<2>& p, const unsigned int component )
{
    realtype x = p(0);
    realtype y = a_ + b_*std::sin( m_*M_PI*x ) + c_*std::sin( n_*M_PI*x );

    if ( smooth_ )
        return (1.0-smallNumber)*std::tanh( (y-p(1)) / smoothXi_ );
    else
        return ( y-p(1) > 0.0 ) ? -1.0 : 1.0;
}



void
TwoScaleSineFunction::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Two-scale sine function ( a + b*sin(m*PI*x) + c*sin(n*PI*x) )" << endl;
    logger << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    logger << "*a*:: " << a_ << endl;
    logger << "*b*:: " << b_ << endl;
    logger << "*c*:: " << c_ << endl;
    logger << "*m*:: " << m_ << endl;
    logger << "*n*:: " << n_ << endl;
    if ( smooth_ )
    {
        logger << "*Smoothed*:: true " << endl;
        logger << "*Interface width*:: " << smoothXi_ << endl;
    }
    else
        logger << "*Smoothed*:: false " << endl;
}


void
TwoScaleSineFunction::declareParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Two-scale sine" );
        prm.declare_entry( "a", "0.5", Patterns::Double() );
        prm.declare_entry( "b", "0.1", Patterns::Double() );
        prm.declare_entry( "c", "0.3", Patterns::Double() );
        prm.declare_entry( "m", "1", Patterns::Integer() );
        prm.declare_entry( "n", "16", Patterns::Integer() );
        prm.declare_entry( "smooth", "yes", Patterns::Bool() );
        prm.declare_entry( "interface width", "0.0625", Patterns::Double() );
    prm.leave_subsection();
}



void TwoScaleSineFunction::getParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Two-scale sine" );
        a_ = prm.get_double( "a" );
        b_ = prm.get_double( "b" );
        c_ = prm.get_double( "c" );
        m_ = prm.get_double( "m" );
        n_ = prm.get_double( "n" );
        smooth_ = prm.get_bool( "smooth" );
        smoothXi_ = prm.get_double( "interface width" );
    prm.leave_subsection();
}

