#ifndef RUNGE_KUTTA_23_H
#define RUNGE_KUTTA_23_H

#include "ExplicitRungeKuttaBase.h"

class RungeKutta23 : public ExplicitRungeKuttaBase {
  public:
    RungeKutta23();
};

OdeIntegratorBase * createRungeKutta23Solver();

#endif // RUNGE_KUTTA_23_H

