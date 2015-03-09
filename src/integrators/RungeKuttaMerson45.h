#ifndef RUNGE_KUTTA_MERSON_45_H
#define RUNGE_KUTTA_MERSON_45_H

#include "ExplicitRungeKuttaBase.h"

class RungeKuttaMerson45 : public ExplicitRungeKuttaBase {
  public:
    RungeKuttaMerson45();
};

OdeIntegratorBase * createRungeKuttaMerson45Solver();

#endif // RUNGE_KUTTA_MERSON_45_H
