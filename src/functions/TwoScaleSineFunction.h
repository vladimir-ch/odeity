#ifndef TWO_SCALE_SINE_FUNCTION_H
#define TWO_SCALE_SINE_FUNCTION_H

#include "Function.h"


class TwoScaleSineFunction : public Function<2>
{
    public:
        TwoScaleSineFunction(
                const realtype a = 0.1,
                const realtype b = 0.3,
                const realtype c = 0.5,
                const int m = 1,
                const int n = 16,
                const bool smooth = true,
                const realtype smoothXi = 0.0625,
                const unsigned int n_components = 1 );

        realtype value( const Point<2>& p, const unsigned int component = 0 );
        void printInfo() const;

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

    private:
        realtype a_;
        realtype b_;
        realtype c_;
        int m_;
        int n_;
        bool smooth_;
        realtype smoothXi_;
};

#include "TwoScaleSineFunction.impl.h"

Function<2> *
createTwoScaleSineFunction()
{
    return new TwoScaleSineFunction();
}

#endif // TWO_SCALE_SINE_FUNCTION_H


