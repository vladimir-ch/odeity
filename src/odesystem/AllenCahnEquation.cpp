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

#include "AllenCahnEquation.h"
#include "../geometry/RectangularGrid.h"
#include "../utils/ParameterHandler.h"
#include "../utils/Vector.h"

#include <cmath>


AllenCahnEquation::AllenCahnEquation()
    :
        MolOdeSystem<2>( false, 1 ),
        xi_( 1.0 ),
        xiSqrInv_( 1.0 ),
        F_( 0.0 )
{}



int AllenCahnEquation::rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot )
{
    realtype u, ul, ur, uu, ud;
    unsigned int xdim = grid_->dimension( xDim );
    unsigned int ydim = grid_->dimension( yDim );
    realtype hxInv2 = 2*grid_->spatialStepInv( xDim );
    realtype hyInv2 = 2*grid_->spatialStepInv( yDim );
    realtype hxPow2Inv = grid_->spatialStepPow2Inv( xDim );
    realtype hyPow2Inv = grid_->spatialStepPow2Inv( yDim );
    int index;
    realtype gradNorm;
    realtype dx, dy;

    for ( unsigned int x = 0; x < xdim; x++ )
        for ( unsigned int y = 0; y < ydim; y++ )
        {
            index = grid_->nodeIndex(x,y);

            u = Y(index);
            ul = ( x == 0 ) ? Y(grid_->nodeIndex(1,y)) : Y(grid_->nodeIndex(x-1,y));
            ur = ( x == xdim-1 ) ? Y(grid_->nodeIndex(xdim-2,y)) : Y(grid_->nodeIndex(x+1,y));
            ud = ( y == 0 ) ? Y(grid_->nodeIndex(x,1)) : Y(grid_->nodeIndex(x,y-1));
            uu = ( y == ydim-1 ) ? Y(grid_->nodeIndex(x,ydim-2)) : Y(grid_->nodeIndex(x,y+1));

            dx = hxInv2*(ur-ul);
            dy = hyInv2*(uu-ud);
            gradNorm = std::sqrt(dx*dx + dy*dy);

            Ydot(index) = hxPow2Inv*( ul - 2*u + ur ) + hyPow2Inv*( uu - 2*u + ud )
                + xiSqrInv_*u*(1.0 - u)*(u + 1.0)
                + F_*gradNorm;
        }

    return 0;
}



void AllenCahnEquation::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Allen-Cahn equation" << endl;
    logger << "-------------------" << endl;
    logger << "*Interface width*:: " << xi_ << endl;
    logger << "*Forcing term*:: " << F_ << endl;
}



    void
AllenCahnEquation::declareParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Allen-Cahn equation");

    prm.declare_entry( "F", "0.0", Patterns::Double(), "Constant forcing term" );
    prm.declare_entry( "xi", "-1", Patterns::Double(), "Interface width" );
    prm.declare_entry( "use analytical Jacobian", "false", Patterns::Bool() );

    prm.leave_subsection();
}



void AllenCahnEquation::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Allen-Cahn equation");

    F_  = prm.get_double( "F" );
    xi_ = prm.get_double( "xi" );
    if ( xi_ < 0.0 )
        xi_ = -xi_ * std::max( grid_->spatialStep( xDim ), grid_->spatialStep( yDim ) );
    xiSqrInv_ = 1.0/(xi_*xi_);

    setUseJacobian( prm.get_bool("use analytical Jacobian") );

    prm.leave_subsection();
}


