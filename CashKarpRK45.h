#ifndef CASH_KARP_RK45_H
#define CASH_KARP_RK45_H

#include "ExplicitRungeKuttaBase.h"

class CashKarpRK45 : public ExplicitRungeKuttaBase {
  public:
    CahsKarpRK45();
};

OdeIntegratorBase * createCahsKarpRK45Solver();

#endif // CASH_KARP_RK45_H

