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


WavyCircleFunction::WavyCircleFunction( const Point<2>& center,
        const realtype radius,
        const realtype waveSize,
        const int numWaves, 
        const bool smooth,
        const realtype smoothXi,
        const unsigned int n_components )
:
    Function<2>( n_components ),
    center_( center ),
    radius_( radius ),
    waveSize_( waveSize ),
    numWaves_( numWaves ),
    smooth_( smooth ),
    smoothXi_( smoothXi )
{}


    
realtype
WavyCircleFunction::value( const Point<2>& p, const unsigned int component )
{
    Point<2> tp = p - center_;
    realtype norm = tp.norm();
    realtype phi = std::atan2( tp(1), tp(0) );
    realtype rho = radius_ + waveSize_ * std::cos( numWaves_ * phi );

    if ( smooth_ )
        return (1.0-smallNumber)*std::tanh( (rho-norm) / smoothXi_ );
    else
        return ( norm - rho > 0.0 ) ? -1.0+smallNumber : 1.0-smallNumber;
}



void
WavyCircleFunction::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Wavy circle function" << endl;
    logger << "~~~~~~~~~~~~~~~~~~~~" << endl;
    logger << "*Center*:: [" << center_( 0 ) << "," << center_( 1 ) << "]" << endl;
    logger << "*Radius*:: " << radius_ << endl;
    logger << "*Peak size*:: " << waveSize_ << endl;
    logger << "*Number of peaks*:: " << numWaves_ << endl;
    if ( smooth_ )
    {
        logger << "*Smoothed*:: true " << endl;
        logger << "*Interface width*:: " << smoothXi_ << endl;
    }
    else
        logger << "*Smoothed*:: false " << endl;
}


void
WavyCircleFunction::declareParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Wavy circle" );
        prm.declare_entry( "center x", "0.5", Patterns::Double() );
        prm.declare_entry( "center y", "0.5", Patterns::Double() );
        prm.declare_entry( "radius", "0.3", Patterns::Double() );
        prm.declare_entry( "peak size", "0.1", Patterns::Double() );
        prm.declare_entry( "number of peaks", "8", Patterns::Integer() );
        prm.declare_entry( "smooth", "yes", Patterns::Bool() );
        prm.declare_entry( "interface width", "0.0625", Patterns::Double() );
    prm.leave_subsection();
}



void WavyCircleFunction::getParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Wavy circle" );
        center_(0) = prm.get_double( "center x" );
        center_(1) = prm.get_double( "center y" );
        waveSize_ = prm.get_double( "peak size" );
        numWaves_ = prm.get_integer( "number of peaks" );
        radius_ = prm.get_double( "radius" );
        smooth_ = prm.get_bool( "smooth" );
        smoothXi_ = prm.get_double( "interface width" );
    prm.leave_subsection();
}

