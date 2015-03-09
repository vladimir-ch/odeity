#ifndef ZERO_FUNCTION_H
#define ZERO_FUNCTION_H

#include "Function.h"


template<int dim>
class ZeroFunction : public Function<dim>
{
  public:
    ZeroFunction( const unsigned int n_components = 1 );
    realtype value ( const Point<dim>& p, const unsigned int component = 0 );
    
    void printInfo() const;
};

#include "ZeroFunction.impl.h"

template<int dim>
Function<dim> *
createZeroFunction()
{
    return new ZeroFunction<dim>();
}

#endif // ZERO_FUNCTION_H
