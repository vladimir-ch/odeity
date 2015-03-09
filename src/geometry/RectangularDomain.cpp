#include "RectangularDomain.h"

template<>
void
RectangularDomain<1>::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Rectangular domain");
        minExtents_[xDim] = prm.get_double( "x-min" );
        maxExtents_[xDim] = prm.get_double( "x-max" );
    prm.leave_subsection();
}


template<>
void
RectangularDomain<2>::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Rectangular domain");
        minExtents_[xDim] = prm.get_double( "x-min" );
        maxExtents_[xDim] = prm.get_double( "x-max" );
        minExtents_[yDim] = prm.get_double( "y-min" );
        maxExtents_[yDim] = prm.get_double( "y-max" );
    prm.leave_subsection();
}


template<>
void
RectangularDomain<3>::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Rectangular domain");
        minExtents_[xDim] = prm.get_double( "x-min" );
        maxExtents_[xDim] = prm.get_double( "x-max" );
        minExtents_[yDim] = prm.get_double( "y-min" );
        maxExtents_[yDim] = prm.get_double( "y-max" );
        minExtents_[zDim] = prm.get_double( "z-min" );
        maxExtents_[zDim] = prm.get_double( "z-max" );
    prm.leave_subsection();
}

