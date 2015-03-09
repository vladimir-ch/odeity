#ifndef RANDOM_FUNCTION_H
#define RANDOM_FUNCTION_H

#include "Function.h"


template<int dim>
class RandomFunction : public Function<dim>
{
    public:

        RandomFunction( const realtype mean = 0.0, const realtype sigma = 1.0, const unsigned int n_components = 1 );
        RandomFunction( const RandomFunction & other )
            :
                mean_( other.mean_ ), sigma_( other.sigma_ ), valid_( false )
        {}

        realtype value ( const Point<dim>& p, const unsigned int component = 0 );

        /**
         * Returns random double between 0 and 1.
         * Based on Wichmann & Hill Applied Stats AS183
         * 1982 Vol 31 188-190.
         */
        realtype uniform();

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

        void printInfo() const;

    private:
        realtype mean_, sigma_;
        realtype r1_, r2_, cachedRho_;
        bool valid_;

        int x,y,z; /// variables for uniform random number generator
};

#include "RandomFunction.impl.h"

template<int dim>
Function<dim> *
createRandomFunction()
{
    return new RandomFunction<dim>();
}

#endif // RANDOM_FUNCTION_H

