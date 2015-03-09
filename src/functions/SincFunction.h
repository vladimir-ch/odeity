#ifndef SINC_FUNCTION_H
#define SINC_FUNCTION_H


#include "Function.h"

template<int dim>
class SincFunction : public Function<dim>
{
  public:

    SincFunction( const unsigned int n_components = 1 );
    ~SincFunction() {}

    realtype value ( const Point<dim>& p, const unsigned int component = 0 );

    void printInfo() const;
};



template<int dim>
Function<dim> *
createSincFunction()
{
    return new SincFunction<dim>();
}

#include "SincFunction.impl.h"

#endif // SINC_FUNCTION_H

