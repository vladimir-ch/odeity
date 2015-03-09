#include "../utils/LogStream.h"

template<int dim>
ZeroFunction<dim>::ZeroFunction( const unsigned int n_components )
    :
        Function<dim>( n_components )
{}



template<int dim>
realtype
ZeroFunction<dim>::value( const Point<dim>& p, const unsigned int component )
{
    return 0.0;
}


template<int N>
void
ZeroFunction<N>::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Zero function" << endl;
    logger << "~~~~~~~~~~~~~" << endl;
}

