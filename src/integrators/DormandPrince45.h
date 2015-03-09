#ifndef DORMAND_PRINCE_45_H
#define DORMAND_PRINCE_45_H

#include "ExplicitRungeKuttaBase.h"

class DormandPrince45 : public ExplicitRungeKuttaBase {
  public:
    DormandPrince45();
};


OdeIntegratorBase * createDormandPrince45Solver();

#endif // DORMAND_PRINCE_45_H
