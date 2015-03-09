#ifndef MOL_ODE_SYSTEM_H
#define MOL_ODE_SYSTEM_H

#include "ExplicitOde.h"
#include "../utils/Exceptions.h"

template<int dim> class RectangularGrid;
class ParameterHandler;

template<int dim>
class MolOdeSystem : public ExplicitOde
{
  public:
    MolOdeSystem( bool hasJacobian, int numComponents );

    int numberOfEquations() const;
    virtual void attachGrid( const RectangularGrid<dim>& grid );

    const RectangularGrid<dim>* grid() const { return grid_; }

    int numberOfComponents() const { return numComponents_; }
    virtual std::string componentName( int i ) const = 0;

    virtual void printInfo() const = 0;
    virtual void getParameters( ParameterHandler& prm ) = 0;

  protected:

    const RectangularGrid<dim>* grid_;
    int numComponents_;

};



template<int dim>
inline
MolOdeSystem<dim>::MolOdeSystem( bool hasJacobian, int numComponents )
  :
    ExplicitOde( hasJacobian ),
    grid_( 0 ),
    numComponents_( numComponents )
{}



template<int dim>
inline
void
MolOdeSystem<dim>::attachGrid( const RectangularGrid<dim>& grid )
{
  grid_ = &grid;
}

#endif // MOD_ODE_SYSTEM_H
