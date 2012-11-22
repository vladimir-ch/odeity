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

template<int dim>
RandomFunction<dim>::RandomFunction( const realtype mean, const realtype sigma, const unsigned int n_components )
    :
        Function<dim>( n_components ),
        mean_( mean ),
        sigma_( sigma ),
        valid_( false ),
        x( 2 ),
        y( 100 ),
        z( 20000 )
{}



template<int dim>
realtype
RandomFunction<dim>::value ( const Point<dim>& p, const unsigned int component )
{
    if ( not valid_ )
    {
        r1_ = uniform();
        r2_ = uniform();
        cachedRho_ = sqrt( - realtype(2) * log(realtype(1) - r2_));
        valid_ = true;
    } else {
        valid_ = false;
    }

    const realtype pi = realtype(3.14159265358979323846);

    return cachedRho_ *( valid_ ? cos(realtype(2)*pi*r1_) : sin(realtype(2)*pi*r1_) )*sigma_ + mean_;
}




template<int dim>
realtype
RandomFunction<dim>::uniform()
{
    realtype ran;

    x = 171 * ( x % 177 ) - 2 * ( x / 177 );
    if ( x < 0 )
        x += 30269;
    y = 172 * ( y % 176 ) - 35 * ( y / 176 );
    if ( y < 0 )
        y += 30307;
    z = 170 * ( z % 178 ) - 63 * ( z / 178 );
    if ( z < 0 )
        z += 30323;
    ran = realtype( x / 30269.0 ) + realtype( y / 30307.0 )
        + realtype( z / 30323.0 ) ;
    ran -= std::floor( ran );

    return ran ;
}



template<int dim>
void
RandomFunction<dim>::declareParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Random" );
        prm.declare_entry( "mean value", "0.0", Patterns::Double() );
        prm.declare_entry( "deviation", "1.0", Patterns::Double() );
    prm.leave_subsection();
}



template<int dim>
void
RandomFunction<dim>::getParameters( ParameterHandler & prm )
{
    prm.enter_subsection( "Random" );
        mean_ = prm.get_double( "mean value" );
        sigma_ = prm.get_double( "deviation" );
    prm.leave_subsection();
}


template<int N>
void
RandomFunction<N>::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Random function" << endl;
    logger << "~~~~~~~~~~~~~~~" << endl;
    logger << "*Mean value*:: " << mean_ << endl;
    logger << "*Standard deviation*:: " << sigma_ << endl;
}
