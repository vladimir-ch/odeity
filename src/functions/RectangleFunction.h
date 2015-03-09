#ifndef RECTANGLE_FUNCTION_H
#define RECTANGLE_FUNCTION_H

#include "Function.h"


class RectangleFunction : public Function<2>
{
    public:
        RectangleFunction(
                const Point<2>& lowerLeftCorner = Point<2>( 0.2, 0.2 ),
                const Point<2>& upperRightCorner = Point<2>( 0.8, 0.8 ),
                const unsigned int n_components = 1 );

        realtype value( const Point<2>& p, const unsigned int component = 0 );
        void printInfo() const;

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

    private:
        Point<2> lowerLeftCorner_;
        Point<2> upperRightCorner_;
};

#include "RectangleFunction.impl.h"

Function<2> *
createRectangleFunction()
{
    return new RectangleFunction();
}

#endif // RECTANGLE_FUNCTION_H


