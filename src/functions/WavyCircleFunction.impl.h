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

