#include "DegenerateCahnHilliardEquation.h"
#include "../geometry/RectangularGrid.h"
#include "../utils/ParameterHandler.h"
#include "../utils/Vector.h"
#include "../utils/LogStream.h"

#include <cmath>
#include <limits>


DegenerateCahnHilliardEquation::DegenerateCahnHilliardEquation()
    :
        MolOdeSystem<2>( false, 1 ),
        xi_( 1.0 ),
        xiPow2Inv_( 1.0 ),
        xiPow2_( 1.0 ),
        xiPowAlpha_( 1.0 ),
        eps_( 10*std::numeric_limits<realtype>::epsilon() )
{}



int
DegenerateCahnHilliardEquation::rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot )
{
    realtype u, ul, ur, uu, ud;
    realtype ml, mr, mu, md;
    unsigned int xdim = grid_->dimension( xDim );
    unsigned int ydim = grid_->dimension( yDim );
    realtype hxPow2Inv = grid_->spatialStepPow2Inv( xDim );
    realtype hyPow2Inv = grid_->spatialStepPow2Inv( yDim );
    unsigned int index, x, y;

    // compute chemical potential
    for ( x = 0; x < xdim; x++ )
        for ( y = 0; y < ydim; y++ )
        {
            index = grid_->nodeIndex(x,y);

            u = Y(index);
            ul = ( x == 0 ) ? Y(grid_->nodeIndex(1,y)) : Y(grid_->nodeIndex(x-1,y));
            ur = ( x == xdim-1 ) ? Y(grid_->nodeIndex(xdim-2,y)) : Y(grid_->nodeIndex(x+1,y));
            ud = ( y == 0 ) ? Y(grid_->nodeIndex(x,1)) : Y(grid_->nodeIndex(x,y-1));
            uu = ( y == ydim-1 ) ? Y(grid_->nodeIndex(x,ydim-2)) : Y(grid_->nodeIndex(x,y+1));

            w_(index) = - xiPow2_*(hxPow2Inv*( ul - 2*u + ur ) + hyPow2Inv*( uu - 2*u + ud )) + f0(u);
        }

    // compute mobility M

    for ( x = 0; x < xdim-1; x++ )
        for ( y = 0; y < ydim; y++ )
        {
            mAvgX_(grid_->nodeIndex(x,y)) = 0.5*( M(Y(grid_->nodeIndex(x,y))) + M(Y(grid_->nodeIndex(x+1,y))) );
            // mAvgX_(grid_->nodeIndex(x,y)) = M(0.5*( Y(grid_->nodeIndex(x,y)) + Y(grid_->nodeIndex(x+1,y)) ));
        }

    for ( x = 0; x < xdim; x++ )
        for ( y = 0; y < ydim-1; y++ )
        {
            mAvgY_(grid_->nodeIndex(x,y)) = 0.5*( M(Y(grid_->nodeIndex(x,y))) + M(Y(grid_->nodeIndex(x,y+1))) );
            // mAvgY_(grid_->nodeIndex(x,y)) = M(0.5*( Y(grid_->nodeIndex(x,y)) + Y(grid_->nodeIndex(x,y+1)) ));
        }

    // compute divergence of M times grad of chempot

    for ( x = 0; x < xdim; x++ )
        for ( y = 0; y < ydim; y++ )
        {
            index = grid_->nodeIndex(x,y);

            u = w_(index);

            if ( x != 0 ) {
                ul = w_(grid_->nodeIndex(x-1,y));
                ml = mAvgX_(grid_->nodeIndex(x-1,y));
            } else {
                ul = w_(grid_->nodeIndex(1,y));
                ml = mAvgX_(grid_->nodeIndex(0,y));
            }
            if ( x != xdim-1 ) {
                ur = w_(grid_->nodeIndex(x+1,y));
                mr = mAvgX_(grid_->nodeIndex(x,y));
            } else {
                ur = w_(grid_->nodeIndex(xdim-2,y));
                mr = mAvgX_(grid_->nodeIndex(x-1,y));
            }
            if ( y != 0 ) {
                ud = w_(grid_->nodeIndex(x,y-1));
                md = mAvgY_(grid_->nodeIndex(x,y-1));
            } else {
                ud = w_(grid_->nodeIndex(x,1));
                md = mAvgY_(grid_->nodeIndex(x,0));
            }
            if ( y != xdim-1 ) {
                uu = w_(grid_->nodeIndex(x,y+1));
                mu = mAvgY_(grid_->nodeIndex(x,y));
            } else {
                uu = w_(grid_->nodeIndex(x,ydim-2));
                mu = mAvgY_(grid_->nodeIndex(x,y-1));
            }

            Ydot(index) = xiPow2Inv_*(
                    hxPow2Inv*( mr*(ur-u) - ml*(u-ul) ) +
                    hyPow2Inv*( mu*(uu-u) - md*(u-ud) ) );
        }

    return 0;
}



inline
realtype DegenerateCahnHilliardEquation::f0( const realtype u )
{
    // return xiPowAlpha_*(log(1.0 + u) - log(1.0 - u)) - u;

    static const double eps_inv = 1.0/eps_;
    static const double lneps = log(eps_);

    if ( u >= 1.0 - eps_ )
        return xiPowAlpha_*(1.0+log(1.0+u) - eps_inv*(1.0-u) - lneps) - u;
    else if ( u <= -1.0 + eps_ )
        return xiPowAlpha_*(-1.0-log(1.0-u) + eps_inv*(1.0+u) + lneps) - u;
    else
        return xiPowAlpha_*(log(1.0 + u) - log(1.0 - u)) - u;
}



inline
realtype DegenerateCahnHilliardEquation::M( const realtype u )
{
    return std::max(beta_*(1.0 - u*u),0.0);
    // if ( u > -1.0 && u < 1.0 )
        // return beta_*(1.0 - u*u);
    // else return 0.0;
}



void
DegenerateCahnHilliardEquation::printInfo( ) const
{
    using namespace std;

    logger << endl;
    logger << "Cahn-Hilliard equation with degenerate mobility" << endl;
    logger << "-----------------------------------------------" << endl;
    logger << "*Mobility parameter (beta)*:: " << beta_ << endl;
    logger << "*Interface width (xi)*:: " << xi_ << endl;
    logger << "*Exponent in log. free energy (alpha)*:: " << alpha_ << endl;
    logger << "*xi^alpha*:: " << 2.0*xiPowAlpha_ << endl;
    logger << "*Regularization parameter in log. free energy*:: " << eps_ << endl;
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
DegenerateCahnHilliardEquation::attachGrid( const RectangularGrid<2>& grid )
{
    MolOdeSystem<2>::attachGrid( grid );

    w_.reinit( grid_->numberOfNodes() );
    mAvgY_.reinit( grid_->numberOfNodes() );
    mAvgX_.reinit( grid_->numberOfNodes() );
}



void
DegenerateCahnHilliardEquation::declareParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Degenerate Cahn-Hilliard equation");

    prm.declare_entry( "xi", "-1.0", Patterns::Double(), "Interface width" );
    prm.declare_entry( "beta", "1.0", Patterns::Double(), "Mobility parameter" );
    prm.declare_entry( "alpha", "1.0", Patterns::Double(), "If positive, it is exponent in log. free energy, if negative, then -alpha is the quench temperature" );
    prm.declare_entry( "eps", "1.0e-3", Patterns::Double(), "Regularization parameter in log. free energy" );

    prm.leave_subsection();
}



void
DegenerateCahnHilliardEquation::getParameters( ParameterHandler& prm )
{
    prm.enter_subsection("Degenerate Cahn-Hilliard equation");

    xi_ = prm.get_double( "xi" );
    if ( xi_ < 0.0 )
    {
        realtype h = std::max( grid_->spatialStep( xDim ), grid_->spatialStep( yDim ) );
        xi_ = -xi_ * h; // * log(h);
    }
    xiPow2_ = xi_*xi_;
    xiPow2Inv_ = 1.0/xiPow2_;

    beta_ = prm.get_double("beta");
    alpha_ = prm.get_double("alpha");

    if ( alpha_ < 0.0 )
        xiPowAlpha_ = -0.5*alpha_;
    else
        xiPowAlpha_ = 0.5*std::pow( xi_, alpha_ );

    eps_ = prm.get_double( "eps" );
    if ( eps_ < 0.0 )
        eps_ = 10*std::numeric_limits<realtype>::epsilon();

    prm.leave_subsection();
}



