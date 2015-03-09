#ifndef WAVY_CIRCLE_FUNCTION_H
#define WAVY_CIRCLE_FUNCTION_H

#include "Function.h"


class WavyCircleFunction : public Function<2>
{
    public:
        WavyCircleFunction( const Point<2>& center = Point<2>( 0.5, 0.5 ),
                const realtype radius = 0.5,
                const realtype waveSize = 0.1,
                const int numWaves = 8, 
                const bool smooth = true,
                const realtype smoothXi = 0.0625,
                const unsigned int n_components = 1 );

        realtype value( const Point<2>& p, const unsigned int component = 0 );
        void printInfo() const;

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

    private:
        Point<2> center_;
        realtype radius_;
        realtype waveSize_;
        int numWaves_;
        bool smooth_;
        realtype smoothXi_;
};

#include "WavyCircleFunction.impl.h"

Function<2> *
createWavyCircleFunction()
{
    return new WavyCircleFunction();
}

#endif // WAVY_CIRCLE_FUNCTION_H

