#include "CahnHilliardEquation.h"
#include "../geometry/RectangularGrid.h"
#include "../utils/ParameterHandler.h"
#include "../utils/Vector.h"
#include "../utils/LogStream.h"

#include <cmath>
#include <omp.h>


CahnHilliardEquation::CahnHilliardEquation()
    :
        MolOdeSystem<2>( true, 1 ),
        xi_( 1.0 ),
        xiInv_( 1.0 ),
        xiSqr_( 1.0 ),
        a_( 1.0 )
{}



    int
CahnHilliardEquation::rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot )
{
    realtype u, ul, ur, uu, ud;
    unsigned int xdim = grid_->dimension( xDim );
    unsigned int ydim = grid_->dimension( yDim );
    realtype hxPow2Inv = grid_->spatialStepPow2Inv( xDim );
    realtype hyPow2Inv = grid_->spatialStepPow2Inv( yDim );
    unsigned int index, x, y;

    for ( x = 0; x < xdim; x++ )
    {
        for ( y = 0; y < ydim; y++ )
        {
            index = grid_->nodeIndex(x,y);

            u = Y(index);
            ul = ( x == 0 ) ? Y(grid_->nodeIndex(1,y)) : Y(grid_->nodeIndex(x-1,y));
            ur = ( x == xdim-1 ) ? Y(grid_->nodeIndex(xdim-2,y)) : Y(grid_->nodeIndex(x+1,y));
            ud = ( y == 0 ) ? Y(grid_->nodeIndex(x,1)) : Y(grid_->nodeIndex(x,y-1));
            uu = ( y == ydim-1 ) ? Y(grid_->nodeIndex(x,ydim-2)) : Y(grid_->nodeIndex(x,y+1));

            w_(index) = - xiSqr_*(hxPow2Inv*( ul - 2*u + ur ) + hyPow2Inv*( uu - 2*u + ud )) + f0(u);
        }
    }

    for ( x = 0; x < xdim; x++ )
    {
        for ( y = 0; y < ydim; y++ )
        {
            index = grid_->nodeIndex(x,y);

            u = w_(index);
            ul = ( x == 0 ) ? w_(grid_->nodeIndex(1,y)) : w_(grid_->nodeIndex(x-1,y));
            ur = ( x == xdim-1 ) ? w_(grid_->nodeIndex(xdim-2,y)) : w_(grid_->nodeIndex(x+1,y));
            ud = ( y == 0 ) ? w_(grid_->nodeIndex(x,1)) : w_(grid_->nodeIndex(x,y-1));
            uu = ( y == ydim-1 ) ? w_(grid_->nodeIndex(x,ydim-2)) : w_(grid_->nodeIndex(x,y+1));

            Ydot(index) = xiInv_*(hxPow2Inv*( ul - 2*u + ur ) + hyPow2Inv*( uu - 2*u + ud ));
        }
    }

    return 0;
}



    inline
realtype CahnHilliardEquation::f0( const realtype u )
{
    // double-obstacle potential ?
    // if ( u >= -1.0 || u <= 1.0 )
    // return -u;
    // else if ( u < -1.0 )
    // return -1.0e10;
    // else return 1.0e10;

    return a_*(u*u - 1.0)*u;
}



int CahnHilliardEquation::jacobian( const Vector<realtype>& v, Vector<realtype>& Jv,
        const realtype t, const Vector<realtype>& Y, const Vector<realtype>& fy,
        const Vector<realtype>& tmp )
{
    unsigned int xdim = grid_->dimension( xDim );
    unsigned int ydim = grid_->dimension( yDim );
    realtype hxPow2Inv = grid_->spatialStepPow2Inv( xDim );
    realtype hyPow2Inv = grid_->spatialStepPow2Inv( yDim );
    realtype hxPow4Inv = grid_->spatialStepPow4Inv( xDim );
    realtype hyPow4Inv = grid_->spatialStepPow4Inv( yDim );
    realtype hxPow2hyPow2Inv = hxPow2Inv*hyPow2Inv;
    unsigned int x, y;

    for ( x = 2; x < xdim-2; x++ )
        for ( y = 2; y < ydim-2; y++ )
        {
            Jv(grid_->nodeIndex(x,y)) =
                -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))

                -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
                +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
                *v(grid_->nodeIndex(x  ,y-1))
                -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))

                -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
                +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
                *v(grid_->nodeIndex(x-1,y  ))
                -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
                *v(grid_->nodeIndex(x  ,y  ))                        
                +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
                *v(grid_->nodeIndex(x+1,y  ))
                -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))

                -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
                +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
                *v(grid_->nodeIndex(x  ,y+1))
                -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))

                -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));
        }

    // lower left corner
    x = 0; y = 0;
    Jv(grid_->nodeIndex(x,y)) =
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)      *v(grid_->nodeIndex(x+2,y  ))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -4.0*(2.0*xi_*hxPow2hyPow2Inv)
        *v(grid_->nodeIndex(x+1,y+1))
        -2.0*(xi_*hyPow4Inv)      *v(grid_->nodeIndex(x  ,y+2));

    x = 1; y = 0;
    Jv(grid_->nodeIndex(x,y)) =
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)      *v(grid_->nodeIndex(x+2,y  ))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)
        *v(grid_->nodeIndex(x+1,y+1))
        -2.0*(xi_*hyPow4Inv)      *v(grid_->nodeIndex(x  ,y+2));

    x = 0; y = 1;
    Jv(grid_->nodeIndex(x,y)) =
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))

        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));

    x = 1; y = 1;
    Jv(grid_->nodeIndex(x,y)) =
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));

    // upper left corner
    x = 0; y = ydim-1;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)      *v(grid_->nodeIndex(x  ,y-2))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -4.0*(2.0*xi_*hxPow2hyPow2Inv)
        *v(grid_->nodeIndex(x+1,y-1))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)      *v(grid_->nodeIndex(x+2,y  ));

    x = 1; y = ydim-1;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)      *v(grid_->nodeIndex(x+2,y  ));

    x = 0; y = ydim-2;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)      *v(grid_->nodeIndex(x  ,y-2))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)      *v(grid_->nodeIndex(x+2,y  ))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1));

    x = 1; y = ydim-2;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)      *v(grid_->nodeIndex(x  ,y-2))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(xi_*hxPow4Inv)      *v(grid_->nodeIndex(x+2,y  ))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1));


    // lower right corner
    x = xdim-1; y = 0;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        -4.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));

    x = xdim-1; y = 1;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));

    x = xdim-2; y = 0;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));

    x = xdim-2; y = 1;
    Jv(grid_->nodeIndex(x,y)) =
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));

    // upper right corner

    x = xdim-1; y=ydim-1;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
        -4.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ));

    x = xdim-1; y=ydim-2;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1));

    x = xdim-2; y=ydim-1;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ));

    x = xdim-2; y=ydim-2;
    Jv(grid_->nodeIndex(x,y)) =
        -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
        *v(grid_->nodeIndex(x  ,y-1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
        -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
        *v(grid_->nodeIndex(x-1,y  ))
        -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
        *v(grid_->nodeIndex(x  ,y  ))                        
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
        *v(grid_->nodeIndex(x+1,y  ))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
        +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
        *v(grid_->nodeIndex(x  ,y+1))
        -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1));

    // lower boundary
    y = 0;
    for ( x = 2; x < xdim-2; x++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
            *v(grid_->nodeIndex(x-1,y  ))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
            *v(grid_->nodeIndex(x+1,y  ))
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
            +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
            *v(grid_->nodeIndex(x  ,y+1))
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
            -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));
    }

    y = 1;
    for ( x = 2; x < xdim-2; x++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
            *v(grid_->nodeIndex(x  ,y-1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
            *v(grid_->nodeIndex(x-1,y  ))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
            *v(grid_->nodeIndex(x+1,y  ))
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
            *v(grid_->nodeIndex(x  ,y+1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
            -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));
    }

    // left boundary
    x = 0;
    for ( y = 2; y < ydim-2; y++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
            *v(grid_->nodeIndex(x  ,y-1))
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
            *v(grid_->nodeIndex(x+1,y  ))
            -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
            *v(grid_->nodeIndex(x  ,y+1))
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));
    }

    x = 1;
    for ( y = 2; y < ydim-2; y++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
            *v(grid_->nodeIndex(x  ,y-1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
            *v(grid_->nodeIndex(x-1,y  ))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
            *v(grid_->nodeIndex(x+1,y  ))
            -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
            *v(grid_->nodeIndex(x  ,y+1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));
    }

    // right boundary
    x = xdim-1;
    for ( y = 2; y < ydim-2; y++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
            *v(grid_->nodeIndex(x  ,y-1))
            -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
            +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
            *v(grid_->nodeIndex(x-1,y  ))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
            *v(grid_->nodeIndex(x  ,y+1))
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));
    }

    x = xdim-2;
    for ( y = 2; y < ydim-2; y++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
            *v(grid_->nodeIndex(x  ,y-1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
            -2.0*(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
            *v(grid_->nodeIndex(x-1,y  ))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
            *v(grid_->nodeIndex(x+1,y  ))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
            *v(grid_->nodeIndex(x  ,y+1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1))
            -(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y+2));
    }

    // top boundary
    y = ydim-1;
    for ( x = 2; x < xdim-2; x++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
            +2.0*(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
            *v(grid_->nodeIndex(x  ,y-1))
            -2.0*(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
            *v(grid_->nodeIndex(x-1,y  ))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
            *v(grid_->nodeIndex(x+1,y  ))
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ));
    }

    y = ydim-2;
    for ( x = 2; x < xdim-2; x++ )
    {
        Jv(grid_->nodeIndex(x,y)) =
            -2.0*(xi_*hyPow4Inv)          *v(grid_->nodeIndex(x  ,y-2))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y-1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y-1))))
            *v(grid_->nodeIndex(x  ,y-1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y-1))
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x-2,y  ))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x-1,y  ))))
            *v(grid_->nodeIndex(x-1,y  ))
            -(6.0*xi_*hxPow4Inv+6.0*xi_*hyPow4Inv+8.0*xi_*hxPow2hyPow2Inv+(2.0/xi_)*f0deriv(Y(grid_->nodeIndex(x,y  ))*(hxPow2Inv+hyPow2Inv)))
            *v(grid_->nodeIndex(x  ,y  ))                        
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hxPow4Inv + (hxPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x+1,y  ))))
            *v(grid_->nodeIndex(x+1,y  ))
            -(xi_*hxPow4Inv)          *v(grid_->nodeIndex(x+2,y  ))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x-1,y+1))
            +(2.0*xi_*hxPow2hyPow2Inv + 4*xi_*hyPow4Inv + (hyPow2Inv/xi_)*f0deriv(Y(grid_->nodeIndex(x  ,y+1))))
            *v(grid_->nodeIndex(x  ,y+1))
            -(2.0*xi_*hxPow2hyPow2Inv)*v(grid_->nodeIndex(x+1,y+1));
    }

    return 0;
}


    inline
realtype CahnHilliardEquation::f0deriv( const realtype u )
{
    //  return 3.0*u*u - 3.0*u + 0.5;
    return a_*(12.0*u*u - 4.0);
}



void
CahnHilliardEquation::printInfo( ) const
{
    using namespace std;

    logger << endl;
    logger << "Cahn-Hilliard equation" << endl;
    logger << "----------------------" << endl;
    logger << "*Interface width*:: " << xi_ << endl;
    logger << "*a*:: " << a_ << endl;
    logger << "*Jacobian implemented*:: ";
    if ( hasJacobian() )
        logger << "true" << endl;
    else
        logger << "false" << endl;
    logger << "*Using Jacobian*:: ";
    if ( useJacobian() )
        logger << "true" << endl;
    else
        logger << "false" << endl;
}



    void
CahnHilliardEquation::attachGrid( const RectangularGrid<2>& grid )
{
    MolOdeSystem<2>::attachGrid( grid );

    w_.reinit( grid_->numberOfNodes() );
}



    void
CahnHilliardEquation::declareParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Cahn-Hilliard equation");
    prm.declare_entry( "xi", "-1", Patterns::Double(), "Interface width" );
    prm.declare_entry( "a", "1", Patterns::Double(), "Coefficient in front of the double-well potential" );
    prm.declare_entry( "use analytical Jacobian", "false", Patterns::Bool() );
    prm.leave_subsection();
}



    void
CahnHilliardEquation::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Cahn-Hilliard equation");

    // TODO: assert that the grid has been assigned

    xi_ = prm.get_double( "xi" );
    if ( xi_ < 0.0 )
        xi_ = -xi_ * std::max( grid_->spatialStep( xDim ), grid_->spatialStep( yDim ) );
    xiSqr_ = xi_*xi_;
    xiInv_ = 1.0/xi_;

    a_ = prm.get_double( "a" );

    setUseJacobian( prm.get_bool("use analytical Jacobian") );

    prm.leave_subsection();
}


