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

