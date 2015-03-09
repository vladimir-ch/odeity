#ifndef TWO_RECTANGLES_FUNCTION_H
#define TWO_RECTANGLES_FUNCTION_H

#include "Function.h"


class TwoRectanglesFunction : public Function<2>
{
    public:
        TwoRectanglesFunction(
                const Point<2>& lowerLeftCorner1 = Point<2>( 0.2, 0.2 ),
                const Point<2>& upperRightCorner1 = Point<2>( 0.4, 0.8 ),
                const Point<2>& lowerLeftCorner2 = Point<2>( 0.6, 0.2 ),
                const Point<2>& upperRightCorner2 = Point<2>( 0.8, 0.8 ),
                const unsigned int n_components = 1 );

        realtype value( const Point<2>& p, const unsigned int component = 0 );
        void printInfo() const;

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

    private:
        Point<2> lowerLeftCorner1_;
        Point<2> upperRightCorner1_;
        Point<2> lowerLeftCorner2_;
        Point<2> upperRightCorner2_;
};

#include "TwoRectanglesFunction.impl.h"

Function<2> *
createTwoRectanglesFunction()
{
    return new TwoRectanglesFunction();
}

#endif // TWO_RECTANGLES_FUNCTION_H


