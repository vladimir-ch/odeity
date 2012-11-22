/* Copyright (c) 2010 Vladimir Chalupecky
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


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
