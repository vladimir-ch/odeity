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

#ifndef RECTANGULAR_GRID_H
#define RECTANGULAR_GRID_H

#include "RectangularDomain.h"
#include "Point.h"
#include "../utils/Exceptions.h"
#include "../utils/ParameterHandler.h"
#include "../utils/LogStream.h"

#include <sundials/sundials_types.h>


template<int dim>
class RectangularGrid
{
  public:

    // static const unsigned int dimension = dim; // commented due to name
    // conflict

    RectangularGrid();
    RectangularGrid( const RectangularDomain<dim>& domain );

    void attachToRectangularDomain( const RectangularDomain<dim>& domain );

    void update();

    unsigned int dimension( const unsigned int index ) const;
    void setDimension( const unsigned int xdim );
    void setDimension( const unsigned int xdim, const unsigned int ydim );
    void setDimension( const unsigned int xdim, const unsigned int ydim, const unsigned int zdim );

    realtype spatialStep( const unsigned int index ) const;
    realtype spatialStepInv( const unsigned int index ) const;
    realtype spatialStepPow2Inv( const unsigned int index ) const;
    realtype spatialStepPow4Inv( const unsigned int index ) const;

    Point<dim> operator () ( const unsigned int xindex ) const;
    Point<dim> operator () ( const unsigned int xindex, const unsigned int yindex ) const;
    Point<dim> operator () ( const unsigned int xindex, const unsigned int yindex, const unsigned int zindex ) const;

    unsigned int numberOfNodes() const;
    unsigned int nodeIndex( const unsigned int xindex ) const;
    unsigned int nodeIndex( const unsigned int xindex, const unsigned int yindex ) const;
    unsigned int nodeIndex( const unsigned int xindex, const unsigned int yindex, const unsigned int zindex ) const;

    static void declareParameters( ParameterHandler &prm );
    virtual void getParameters( ParameterHandler &prm );

    virtual void printInfo() const;

    DeclException1( ExcIncorrectDimension,
        int,
        << "Incorrect dimension " << arg1 <<". Correct dimensions are 1, 2 or 3." );

  private:

    unsigned int numNodes_[dim];
    realtype spatialStep_[dim];
    realtype spatialStepInv_[dim];
    realtype spatialStepPow2Inv_[dim];
    realtype spatialStepPow4Inv_[dim];

    bool updated_;

    const RectangularDomain<dim>* domain_;
};



template<int dim>
inline
RectangularGrid<dim>::RectangularGrid() :
  updated_( false ),
  domain_( 0 )
{
  Assert( 1 <= dim && dim <= 3, ExcIncorrectDimension( dim ) );

  for( unsigned int i = 0; i < dim; i++ )
  {
    numNodes_[i] = 0;
    spatialStep_[i] = 0.0;
  }  
}



template<int dim>
inline
RectangularGrid<dim>::RectangularGrid( const RectangularDomain<dim>& domain )
{
  Assert( 1 <= dim && dim <= 3, ExcIncorrectDimension( dim ) );

  attachToRectangularDomain( domain );
}



template<int dim>
inline
void
RectangularGrid<dim>::attachToRectangularDomain( const RectangularDomain<dim>& domain )
{
  domain_ = &domain;
  updated_ = false;
}



template<int dim>
inline
unsigned int
RectangularGrid<dim>::dimension( const unsigned int index ) const
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );
  return numNodes_[index];
}



template<>
inline
void
RectangularGrid<1>::setDimension( const unsigned int xdim )
{
  Assert( dim == 1, ExcImpossibleInDim(dim) );
  numNodes_[0] = xdim;
  updated_ = false;
}



template<>
inline
void
RectangularGrid<2>::setDimension( const unsigned int xdim, const unsigned int ydim )
{
  Assert( dim == 2, ExcImpossibleInDim(dim) );
  numNodes_[0] = xdim;
  numNodes_[1] = ydim;
  updated_ = false;
}



template<>
inline
void
RectangularGrid<3>::setDimension( const unsigned int xdim, const unsigned int ydim, const unsigned int zdim )
{
  Assert( dim == 3, ExcImpossibleInDim(dim) );
  numNodes_[0] = xdim;
  numNodes_[1] = ydim;
  numNodes_[2] = zdim;
  updated_ = false;
}



template<int dim>
inline
realtype
RectangularGrid<dim>::spatialStep( const unsigned int index ) const
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );

  return spatialStep_[index];
}



template<int dim>
inline
realtype
RectangularGrid<dim>::spatialStepInv( const unsigned int index ) const
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );

  return spatialStepInv_[index];
}



template<int dim>
inline
realtype
RectangularGrid<dim>::spatialStepPow2Inv( const unsigned int index ) const
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );

  return spatialStepPow2Inv_[index];
}



template<int dim>
inline
realtype
RectangularGrid<dim>::spatialStepPow4Inv( const unsigned int index ) const
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );

  return spatialStepPow4Inv_[index];
}



template<int dim>
inline
void
RectangularGrid<dim>::update()
{
  Assert( domain_ != 0, ExcMessage("Grid is not attached to any domain") );

  for( unsigned int i = 0; i < dim; i++ )
  {
    Assert( numNodes_[i] != 0, ExcZero() );
    spatialStep_[i] = domain_->size( i ) / ( numNodes_[i] - 1 );
    spatialStepInv_[i] = 1.0 / spatialStep_[i];
    spatialStepPow2Inv_[i] = 1.0 / (spatialStep_[i]*spatialStep_[i]);
    spatialStepPow4Inv_[i] = 1.0 / (spatialStep_[i]*spatialStep_[i]*spatialStep_[i]*spatialStep_[i]);
  }

  updated_ = true;
}



template<int dim>
inline
Point<dim>
RectangularGrid<dim>::operator () ( const unsigned int xindex ) const
{
  Assert( dim == 1, ExcImpossibleInDim(dim) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );
  Assert( xindex < numNodes_[0], ExcIndexRange( xindex, 0, numNodes_[0] ) );

  return Point<dim>( domain_->minExtent(0) + xindex*spatialStep_[0] );
}



template<int dim>
inline
Point<dim>
RectangularGrid<dim>::operator () ( const unsigned int xindex, const unsigned int yindex ) const
{
  Assert( dim == 2, ExcImpossibleInDim(dim) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );
  Assert( xindex < numNodes_[0], ExcIndexRange( xindex, 0, numNodes_[0] ) );
  Assert( yindex < numNodes_[1], ExcIndexRange( yindex, 0, numNodes_[1] ) );

  return Point<dim>( domain_->minExtent(0) + xindex*spatialStep_[0],
                     domain_->minExtent(1) + yindex*spatialStep_[1] );
}



template<int dim>
inline
Point<dim>
RectangularGrid<dim>::operator () ( const unsigned int xindex, const unsigned int yindex, const unsigned int zindex ) const
{
  Assert( dim == 3, ExcImpossibleInDim(dim) );
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );
  Assert( xindex < numNodes_[0], ExcIndexRange( xindex, 0, numNodes_[0] ) );
  Assert( yindex < numNodes_[1], ExcIndexRange( yindex, 0, numNodes_[1] ) );
  Assert( zindex < numNodes_[2], ExcIndexRange( zindex, 0, numNodes_[2] ) );

  return Point<dim>( domain_->minExtent(0) + xindex*spatialStep_[0],
                     domain_->minExtent(1) + yindex*spatialStep_[1],
                     domain_->minExtent(2) + zindex*spatialStep_[2] );
}



template<int dim>
inline
unsigned int
RectangularGrid<dim>::numberOfNodes() const
{
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );

  unsigned int num = 1;
  for( int i = 0; i < dim; i++ )
    num *= numNodes_[i];
  return num;
}



template<int dim>
inline
unsigned int
RectangularGrid<dim>::nodeIndex( const unsigned int xindex ) const
{
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );
  Assert( dim == 1, ExcImpossibleInDim(dim) );
  Assert( xindex < numNodes_[0], ExcIndexRange( xindex, 0, numNodes_[0] ) );

  return xindex;
}



template<int dim>
inline
unsigned int
RectangularGrid<dim>::nodeIndex( const unsigned int xindex, const unsigned int yindex ) const
{
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );
  Assert( dim == 2, ExcImpossibleInDim(dim) );
  Assert( xindex < numNodes_[0], ExcIndexRange( xindex, 0, numNodes_[0] ) );
  Assert( yindex < numNodes_[1], ExcIndexRange( yindex, 0, numNodes_[1] ) );

  return xindex*numNodes_[1] + yindex;
}



template<int dim>
inline
unsigned int
RectangularGrid<dim>::nodeIndex( const unsigned int xindex, const unsigned int yindex, const unsigned int zindex ) const
{
  Assert( updated_ == true, ExcMessage("Grid is not up to date") );
  Assert( dim == 3, ExcImpossibleInDim(dim) );
  Assert( xindex < numNodes_[0], ExcIndexRange( xindex, 0, numNodes_[0] ) );
  Assert( yindex < numNodes_[1], ExcIndexRange( yindex, 0, numNodes_[1] ) );
  Assert( zindex < numNodes_[2], ExcIndexRange( zindex, 0, numNodes_[2] ) );

  return (xindex*numNodes_[1] + yindex)*numNodes_[2] + zindex;
}



template<int dim>
void
RectangularGrid<dim>::declareParameters( ParameterHandler &prm )
{
    prm.enter_subsection("Rectangular grid");

    if ( dim == 1 || dim == 2 || dim == 3 )
        prm.declare_entry( "x-dimension", "100", Patterns::Integer(), "Number of grid nodes in x" );
    if ( dim == 2 || dim == 3 )
        prm.declare_entry( "y-dimension", "100", Patterns::Integer(), "Number of grid nodes in y" );
    if ( dim == 3 )
        prm.declare_entry( "z-dimension", "100", Patterns::Integer(), "Number of grid nodes in z" );

    prm.leave_subsection();
}



template<int dim>
inline
void
RectangularGrid<dim>::getParameters( ParameterHandler &prm )
{
  prm.enter_subsection("Rectangular grid");

  unsigned int xdim;
  unsigned int ydim;
  unsigned int zdim;

  if ( dim == 1 || dim == 2 || dim == 3 )
  {
    xdim = prm.get_integer( "x-dimension");
    if ( dim == 1 )
      setDimension( xdim );
  }
  if ( dim == 2 || dim == 3 )
  {
    ydim = prm.get_integer( "y-dimension");
    if ( dim == 2 )
      setDimension( xdim, ydim );
  }
  if ( dim == 3 )
  {
    zdim = prm.get_integer( "z-dimension");
    setDimension( xdim, ydim, zdim );
  }

  prm.leave_subsection();
}



template<int dim>
inline
void
RectangularGrid<dim>::printInfo() const
{
  using namespace std;

  logger << endl;
  logger << "Rectangular grid" << endl;
  logger << "----------------" << endl;

  if ( updated_ )
  {
    logger << "*Dimension*:: " << dim << endl;
    logger << "*Extents*:: ";
    logger << "[ ";
    for( unsigned int i = 0; i < dim ; i++ )
    {
      logger << numNodes_[i];
      if ( i != dim - 1 )
        logger << ", ";
    }
    logger << " ]" << endl;
    logger << "*Spatial steps*:: ";
    logger << "[ ";
    for( unsigned int i = 0; i < dim ; i++ )
    {
      logger << spatialStep_[i];
      if ( i != dim - 1 )
        logger << ", ";
    }
    logger << " ]" << endl;
  }
  else
  {
    logger << "Grid is not updated!" << endl;
  }
}


#endif // RECTANGULAR_GRID_H

