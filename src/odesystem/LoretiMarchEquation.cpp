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

#include "LoretiMarchEquation.h"
#include "../geometry/RectangularGrid.h"
#include "../utils/ParameterHandler.h"
#include "../utils/Vector.h"
#include "../utils/LogStream.h"

#include <cmath>


LoretiMarchEquation::LoretiMarchEquation()
    :
        MolOdeSystem<2>( false, 1 ),
        xi_( 1.0 ),
        xiInv_( 1.0 ),
        xiSqr_( 1.0 ),
        twoOverXi_( 2.0 ),
        twoOverXiPow3_( 2.0 ),
        a_( 1.0 )
{}



int
LoretiMarchEquation::rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot )
{
    realtype u, ul, ur, uu, ud, w, wl, wr, wd, wu;
    unsigned int xdim = grid_->dimension( xDim );
    unsigned int ydim = grid_->dimension( yDim );
    realtype hxPow2Inv = grid_->spatialStepPow2Inv( xDim );
    realtype hyPow2Inv = grid_->spatialStepPow2Inv( yDim );
    unsigned int index, x, y;

    for ( x = 0; x < xdim; x++ )
        for ( y = 0; y < ydim; y++ )
        {
            index = grid_->nodeIndex(x,y);

            u = Y(index);
            ul = ( x == 0 ) ? Y(grid_->nodeIndex(1,y)) : Y(grid_->nodeIndex(x-1,y));
            ur = ( x == xdim-1 ) ? Y(grid_->nodeIndex(xdim-2,y)) : Y(grid_->nodeIndex(x+1,y));
            ud = ( y == 0 ) ? Y(grid_->nodeIndex(x,1)) : Y(grid_->nodeIndex(x,y-1));
            uu = ( y == ydim-1 ) ? Y(grid_->nodeIndex(x,ydim-2)) : Y(grid_->nodeIndex(x,y+1));

            w_(index) = - xiSqr_*(hxPow2Inv*( ul - 2*u + ur ) + hyPow2Inv*( uu - 2*u + ud )) + PsiDer(u);
        }

    for ( x = 0; x < xdim; x++ )
        for ( y = 0; y < ydim; y++ )
        {
            index = grid_->nodeIndex(x,y);

            w = w_(index);
            u = Y(index);
            wl = ( x == 0 ) ? w_(grid_->nodeIndex(1,y)) : w_(grid_->nodeIndex(x-1,y));
            wr = ( x == xdim-1 ) ? w_(grid_->nodeIndex(xdim-2,y)) : w_(grid_->nodeIndex(x+1,y));
            wd = ( y == 0 ) ? w_(grid_->nodeIndex(x,1)) : w_(grid_->nodeIndex(x,y-1));
            wu = ( y == ydim-1 ) ? w_(grid_->nodeIndex(x,ydim-2)) : w_(grid_->nodeIndex(x,y+1));

            Ydot(index) = twoOverXi_*(hxPow2Inv*( wl - 2*w + wr ) + hyPow2Inv*( wu - 2*w + wd ))
                - xiInv_*w - twoOverXiPow3_ * PsiDerDer(u)*w;
        }

    return 0;
}



inline
realtype LoretiMarchEquation::PsiDer( const realtype u )
{
    return a_*(u*u - 1.0)*u;
}


inline
realtype LoretiMarchEquation::PsiDerDer( const realtype u )
{
    return a_*(3.0*u*u - 1.0);
}



void
LoretiMarchEquation::printInfo( ) const
{
    using namespace std;

    logger << endl;
    logger << "Loreti-March equation" << endl;
    logger << "---------------------" << endl;
    logger << "*Interface width*:: " << xi_ << endl;
    logger << "*a*:: " << a_ << endl;
    logger << "*Jacobian implemented*:: false" << endl;
}



void
LoretiMarchEquation::attachGrid( const RectangularGrid<2>& grid )
{
    MolOdeSystem<2>::attachGrid( grid );

    w_.reinit( grid_->numberOfNodes() );
}



void
LoretiMarchEquation::declareParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Loreti-March equation");
    prm.declare_entry( "xi", "-1", Patterns::Double(), "Interface width" );
    prm.declare_entry( "a", "4", Patterns::Double(), "Coefficient in front of the double-well potential" );
    prm.leave_subsection();
}



    void
LoretiMarchEquation::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Loreti-March equation");

    // TODO: assert that the grid has been assigned

    xi_ = prm.get_double( "xi" );
    if ( xi_ < 0.0 )
        xi_ = -xi_ * std::max( grid_->spatialStep( xDim ), grid_->spatialStep( yDim ) );
    xiSqr_ = xi_*xi_;
    xiInv_ = 1.0/xi_;
    twoOverXi_ = 2.0/xi_;
    twoOverXiPow3_ = 2.0/(xi_*xi_*xi_);

    a_ = prm.get_double( "a" );

    prm.leave_subsection();
}


