#ifndef RECTANGULAR_DOMAIN_H
#define RECTANGULAR_DOMAIN_H

#include "../utils/Exceptions.h"
#include "../utils/ParameterHandler.h"
#include "../utils/LogStream.h"

#include <sundials/sundials_types.h>

const unsigned int xDim = 0;
const unsigned int yDim = 1;
const unsigned int zDim = 2;



template<int dim>
class RectangularDomain
{
    public:

        static const unsigned int dimension = dim;

        RectangularDomain();

        void setExtents( const unsigned int d, const realtype minExtent, const realtype maxExtent );
        realtype minExtent( const unsigned int d ) const;
        realtype maxExtent( const unsigned int d ) const;

        realtype size( const unsigned int d ) const;

        virtual void printInfo() const;

        static void declareParameters( ParameterHandler& prm );
        virtual void getParameters( ParameterHandler& prm );

        DeclException1( ExcIncorrectDimension,
                int,
                << "Incorrect dimension " << arg1 <<". Correct dimensions are 1, 2 or 3." );

    private:

        realtype minExtents_[dim];
        realtype maxExtents_[dim];

};



template<int dim>
RectangularDomain<dim>::RectangularDomain()
{
    Assert( 1 <= dim && dim <= 3, ExcIncorrectDimension( dim ) );

    for( unsigned int i = 0; i < dim; i++ )
    {
        minExtents_[i] = 0.0;
        maxExtents_[i] = 1.0;
    }
}



template<int dim>
inline
void
RectangularDomain<dim>::setExtents( const unsigned int d, const realtype minExtent, const realtype maxExtent )
{
    Assert( d < dim, ExcIndexRange( index, 0, dim ) );

    minExtents_[d] = minExtent;
    maxExtents_[d] = maxExtent;
}




template<int dim>
inline
realtype
RectangularDomain<dim>::minExtent( const unsigned int d ) const
{
    Assert( d < dim, ExcIndexRange( d, 0, dim ) );
    return minExtents_[d];
}



template<int dim>
inline
realtype
RectangularDomain<dim>::maxExtent( const unsigned int d ) const
{
    Assert( d < dim, ExcIndexRange( d, 0, dim ) );
    return maxExtents_[d];
}



template<int dim>
inline
realtype
RectangularDomain<dim>::size( const unsigned int d ) const
{
    Assert( d < dim, ExcIndexRange( d, 0, dim ) );
    return maxExtents_[d] - minExtents_[d];
}



template<int dim>
void
RectangularDomain<dim>::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Rectangular domain" << endl;
    logger << "------------------" << endl;
    logger << "*Dimension*:: " << dim << endl;
    logger << "*Extents*:: ";
    for( unsigned int i = 0; i < dim; i++ )
    {
        logger << "[ " << minExtents_[i] << ", " << maxExtents_[i] << " ]";
        if( i != dim-1 )
            logger << " x ";
    }
    logger << endl;
}



template<int dim>
void
RectangularDomain<dim>::declareParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Rectangular domain");
        prm.declare_entry( "x-min", "0.0", Patterns::Double(), "Lower x extent of the domain" );
        prm.declare_entry( "x-max", "1.0", Patterns::Double(), "Upper x extent of the domain" );
        prm.declare_entry( "y-min", "0.0", Patterns::Double(), "Lower y extent of the domain" );
        prm.declare_entry( "y-max", "1.0", Patterns::Double(), "Upper y extent of the domain" );
        prm.declare_entry( "z-min", "0.0", Patterns::Double(), "Lower z extent of the domain" );
        prm.declare_entry( "z-max", "1.0", Patterns::Double(), "Upper z extent of the domain" );
    prm.leave_subsection();
}




#endif // RECTANGULAR_DOMAIN_H

