#ifndef ODE_SYSTEM_H
#define ODE_SYSTEM_H

#include <string>

class OdeSystemBase
{

    public:
        OdeSystemBase() {}

        virtual int numberOfEquations() const = 0;
        virtual std::string name() const = 0;
};

#endif // ODE_SYSTEM_H

