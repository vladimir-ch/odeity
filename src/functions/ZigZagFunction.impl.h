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

