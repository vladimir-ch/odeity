#ifndef FUNCTION_H
#define FUNCTION_H

#include "FunctionTime.h"
#include "../geometry/Point.h"
#include "../utils/ParameterHandler.h"


template<int dim>
class Function : public FunctionTime
{
    public:

        static const unsigned int dimension = dim;
        static const double smallNumber = 1.0e-2;

        const unsigned int numberOfComponents;

        Function( const unsigned int n_components = 1, const realtype initial_time = 0.0 )
            :
                FunctionTime( initial_time ),
                numberOfComponents( n_components )
        {
            Assert( n_components > 0, ExcZero() );
        }

        virtual ~Function()
        {}

        Function & operator = ( const Function & f );
        virtual realtype value( const Point<dim>& p, const unsigned int component = 0 );
        realtype operator () ( const Point<dim>& p, const unsigned int component = 0 )
        {
            return value( p, component );
        }

        virtual void getParameters( ParameterHandler & prm ) {}
        virtual void printInfo() const {}

        DeclException2 (ExcNumberOfComponents,
                int, int,
                << "You can only assign function objects with the same "
                << "number of components. However, here the left operand "
                << "has " << arg1 << " components, and the right operand "
                << arg2 << " components.");
};

#endif // FUNCTION_H

