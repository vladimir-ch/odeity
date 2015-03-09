#include "Function.h"


template<int dim>
Function<dim> &
Function<dim>::operator= ( const Function<dim>& f )
{
  Assert( numberOfComponents == f.numberOfComponents,
      ExcNumberOfComponents( numberOfComponents, f.numberOfComponents ) );

  return *this;
}



template<int dim>
realtype
Function<dim>::value( const Point<dim>& p, const unsigned int component )
{
  Assert (false, ExcPureFunctionCalled());
  return 0.0;
}

template class Function<1>;
template class Function<2>;
template class Function<3>;

