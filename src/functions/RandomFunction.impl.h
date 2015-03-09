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
