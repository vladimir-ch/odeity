#include "../utils/LogStream.h"

template<int dim>
SincFunction<dim>::SincFunction( const unsigned int n_components )
    :
        Function<dim>( n_components )
{}



template<int dim>
realtype
SincFunction<dim>::value( const Point<dim>& p, const unsigned int component )
{
    realtype dist = p.norm();
    if ( dist == 0.0 )
        return 1.0;

    return std::sin(dist)/dist;
}


template<int N>
void
SincFunction<N>::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Sinc function" << endl;
    logger << "~~~~~~~~~~~~~" << endl;
}

