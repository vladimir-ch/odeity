#include "MolOdeSystem.h"
#include "../geometry/RectangularGrid.h"


template<int dim>
int
MolOdeSystem<dim>::numberOfEquations() const
{
  Assert( grid_ != 0, ExcNotInitialized() );
  return grid_->numberOfNodes();
}

template class MolOdeSystem<1>;
template class MolOdeSystem<2>;
template class MolOdeSystem<3>;

