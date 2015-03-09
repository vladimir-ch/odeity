#ifndef ZIG_ZAG_FUNCTION_H
#define ZIG_ZAG_FUNCTION_H


#include "Function.h"

class ZigZagFunction : public Function<2>
{
  public:

    ZigZagFunction( const unsigned int n_components = 1 );
    ~ZigZagFunction() {}

    realtype value ( const Point<2>& p, const unsigned int component = 0 );

    void printInfo() const;
};



Function<2> *
createZigZagFunction()
{
    return new ZigZagFunction();
}

#include "ZigZagFunction.impl.h"

#endif // ZIG_ZAG_FUNCTION_H

