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

#include "CVode.h"
#include "../odesystem/ExplicitOde.h"
#include "../utils/Vector.h"
#include "../utils/LogStream.h"

#include <cvode/cvode.h>                  /* main integrator header file */
#include <cvode/cvode_spgmr.h>            /* prototypes & constants for CVSPGMR solver */
#include <cvode/cvode_spbcgs.h>           /* prototypes & constants for CVSPBCG solver */
#include <cvode/cvode_sptfqmr.h>          /* prototypes & constants for CVSTFQMR solver */
#include <cvode/cvode_bandpre.h>


extern "C"
int ExplicitOdeRhs( double t, N_Vector nvY, N_Vector nvYdot, void * f_data )
{
    Vector<realtype> y( nvY );
    Vector<realtype> ydot( nvYdot );
    ExplicitOde * explicitOde = ( ExplicitOde * ) f_data;

    return explicitOde->rhs( t, y , ydot );
}



extern "C"
int ExplicitOdeJacobian( N_Vector nvV, N_Vector nvJv,
        realtype t, N_Vector nvY, N_Vector nvFy,
        void *jac_data, N_Vector nvTmp )
{
    Vector<realtype> v( nvV );
    Vector<realtype> Jv( nvJv );
    Vector<realtype> y( nvY );
    Vector<realtype> fy( nvFy );
    Vector<realtype> tmp( nvTmp );
    ExplicitOde * explicitOde = ( ExplicitOde * ) jac_data;

    return explicitOde->jacobian( v, Jv, t, y, fy, tmp );
}



OdeIntegratorBase * createCVodeGMRESSolver()
{
    return new CVodeGMRES();
}



OdeIntegratorBase * createCVodeBiCGSolver()
{
    return new CVodeBiCG();
}



OdeIntegratorBase * createCVodeTFQMRSolver()
{
    return new CVodeTFQMR();
}



// ------------------------------------------------------------------ CVodeBase

CVodeBase::CVodeBase()
    :
        initialized_( false ),
        maxNumSteps_( 10000 )
{
    cvodeMem_ = CVodeCreate( CV_BDF, CV_NEWTON );
    AssertThrow( cvodeMem_ != 0, ExcCVodeCreateError( "CVodeCreate failed" ) );
}



CVodeBase::~CVodeBase()
{
    if ( initialized_ )
        CVodeFree( &cvodeMem_ );
}



void
CVodeBase::integrateTo( const realtype tOut )
{
    if ( saveHistory() )
    {
        int flag = CVodeSetStopTime( cvodeMem_, tOut );
        AssertThrow( flag == CV_SUCCESS, ExcCVodeError( flag ) );

        while( currentTime_ != tOut )
        {
            flag = CVode( cvodeMem_, tOut, nvCurrentState_, &currentTime_, CV_ONE_STEP );
            AssertThrow( flag == CV_SUCCESS || flag == CV_TSTOP_RETURN, ExcCVodeError( flag ) );

            updateHistory();
        }
    }
    else
    {
        int flag = CVode( cvodeMem_, tOut, nvCurrentState_, &currentTime_, CV_NORMAL );
        AssertThrow( flag == CV_SUCCESS, ExcCVodeError( flag ) );
    }
}



void
CVodeBase::setRelativeTolerance( realtype relTol )
{
    OdeIntegratorBase::setRelativeTolerance( relTol );
    int flag;
    flag = CVodeSStolerances( cvodeMem_, relTol_, absTol_ );
    Assert( flag == CV_SUCCESS, ExcCVodeSetTolerancesError( flag ) );
}



void
CVodeBase::setAbsoluteTolerance( realtype absTol )
{
    OdeIntegratorBase::setAbsoluteTolerance( absTol );
    int flag;
    flag = CVodeSStolerances( cvodeMem_, relTol_, absTol_ );
    Assert( flag == CV_SUCCESS, ExcCVodeSetTolerancesError( flag ) );
}



void
CVodeBase::assignExplicitOde( 
        ExplicitOde& odeProblem,
        const realtype initialTime,
        const Vector<realtype> &initialState )
{
    OdeIntegratorBase::assignExplicitOde( odeProblem, initialTime, initialState );

    nvCurrentState_ = N_VMake_Serial( currentState_.size(), currentState_.data() );

    if ( initialized_ )
        CVodeFree( &cvodeMem_ );

    int flag = CVodeInit( cvodeMem_, ExplicitOdeRhs, currentTime_, nvCurrentState_ );
    AssertThrow( flag == CV_SUCCESS, ExcCVodeInitError( flag ) );
    initialized_ = true;

    setRelativeTolerance( relTol_ );
    setAbsoluteTolerance( absTol_ );

    flag = CVodeSetUserData( cvodeMem_, (void*) odeProblem_ );
    Assert( flag == CV_SUCCESS, ExcCVodeSetUserDataError( flag ) );

    flag = CVodeSetMaxNumSteps( cvodeMem_, maxNumSteps_ );
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVodeSetMaxNumSteps failed"));
    
    flag = CVodeSetMaxErrTestFails( cvodeMem_, 20 );
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVodeSetMaxErrTestFails failed"));
    
    flag = CVodeSetMaxConvFails( cvodeMem_, 50 );
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVodeSetMaxConvFails failed"));
}



void
CVodeBase::declareParameters( ParameterHandler & prm )
{
    OdeIntegratorBase::declareParameters( prm );

    prm.enter_subsection( "ODE integrator" );
    prm.enter_subsection( "CVode" );
        prm.declare_entry( "maximum number of steps", "100000", Patterns::Integer() );
    prm.leave_subsection();
    prm.leave_subsection();
}



void
CVodeBase::getParameters( ParameterHandler & prm )
{
    OdeIntegratorBase::getParameters( prm );

    prm.enter_subsection( "ODE integrator" );
    prm.enter_subsection( "CVode" );
        maxNumSteps_ = prm.get_integer( "maximum number of steps" );
    prm.leave_subsection();
    prm.leave_subsection();

}



void
CVodeBase::printInfo() const
{
    OdeIntegratorBase::printInfo();

    logger << "*Maximum number of steps to t_out*:: " << maxNumSteps_ << std::endl;
}



// ----------------------------------------------------------------- CVodeSpils
CVodeSpils::CVodeSpils()
    :
        precond_( false ),
        krylovSubspaceDim_( 5 )
{}



CVodeSpils::~CVodeSpils()
{}


void
CVodeSpils::assignExplicitOde(
        ExplicitOde& odeProblem,
        const realtype initialTime,
        const Vector<realtype> &initialState )
{
    CVodeBase::assignExplicitOde( odeProblem, initialTime, initialState );
    stats_.reset();
}



IntegratorStatsBase&
CVodeSpils::stats()
{
    stats_.update( cvodeMem_ );
    return stats_;
}



void
CVodeSpils::updateHistory()
{
    realtype hlast;
    int flag = CVodeGetLastStep( cvodeMem_, &hlast);

    int qlast;
    flag = CVodeGetLastOrder( cvodeMem_, &qlast);

    stats_.updateHistory( currentTime_, hlast, qlast );
}


void
CVodeSpils::declareParameters( ParameterHandler & prm )
{
    CVodeBase::declareParameters( prm );

    prm.enter_subsection( "ODE integrator" );
    prm.enter_subsection( "CVode" );
    prm.enter_subsection( "CVodeSpils" );
        prm.declare_entry( "precondition", "false", Patterns::Bool() );
        prm.declare_entry( "Krylov subspace dimension", "30", Patterns::Integer() );
    prm.leave_subsection();
    prm.leave_subsection();
    prm.leave_subsection();
}



void
CVodeSpils::getParameters( ParameterHandler & prm )
{
    CVodeBase::getParameters( prm );

    prm.enter_subsection( "ODE integrator" );
    prm.enter_subsection( "CVode" );
    prm.enter_subsection( "CVodeSpils" );
        precond_ = prm.get_bool( "precondition" );
        krylovSubspaceDim_ = prm.get_integer( "Krylov subspace dimension" );
    prm.leave_subsection();
    prm.leave_subsection();
    prm.leave_subsection();
}



void
CVodeSpils::printInfo() const
{
    CVodeBase::printInfo();

    logger << "*Maximum dimension of Krylov subspace*:: " << krylovSubspaceDim_ << std::endl;
    logger << "*Using preconditioning BANDPRE*:: " << precond_ << std::endl;
}



// ----------------------------------------------------------------- CVodeGMRES

CVodeGMRES::CVodeGMRES()
{
    solverName_ = "CVODE with GMRES";
}



void
CVodeGMRES::assignExplicitOde( 
        ExplicitOde& odeProblem,
        const realtype initialTime,
        const Vector<realtype> &initialState )
{
    CVodeSpils::assignExplicitOde( odeProblem, initialTime, initialState );

    if ( not precond_ )
    {
        int flag = CVSpgmr( cvodeMem_, PREC_NONE, krylovSubspaceDim_ );
        AssertThrow( flag == CV_SUCCESS, ExcCVSpgmrError( flag ) );
    }
    else
    {
        int flag = CVSpgmr( cvodeMem_, PREC_RIGHT, krylovSubspaceDim_ );
        AssertThrow( flag == CV_SUCCESS, ExcCVSpgmrError( flag ) );

        flag = CVBandPrecInit( cvodeMem_, odeProblem_->numberOfEquations(), 1, 1 );
        AssertThrow( flag == CVSPILS_SUCCESS, ExcMessage("Failed to allocate BAND preconditioner") );
    }

    if ( odeProblem_->useJacobian() )
    {
        int flag;
        flag = CVSpilsSetJacTimesVecFn( cvodeMem_, ExplicitOdeJacobian );
        Assert( flag == CV_SUCCESS, ExcMessage("Call to SpilsSetJactimesVecFn failed"));
    }
}



// ------------------------------------------------------------------ CVodeBiCG

CVodeBiCG::CVodeBiCG()
{
    solverName_ = "CVODE with Bi-CGStab";
}



void
CVodeBiCG::assignExplicitOde( 
        ExplicitOde& odeProblem,
        const realtype initialTime,
        const Vector<realtype> &initialState )
{
    CVodeSpils::assignExplicitOde( odeProblem, initialTime, initialState );

    if ( not precond_ )
    {
        int flag = CVSpbcg( cvodeMem_, PREC_NONE, krylovSubspaceDim_ );
        AssertThrow( flag == CV_SUCCESS, ExcCVSpbcgError( flag ) );
    }
    else
    {
        int flag = CVSpbcg( cvodeMem_, PREC_RIGHT, krylovSubspaceDim_ );
        AssertThrow( flag == CV_SUCCESS, ExcCVSpbcgError( flag ) );

        flag = CVBandPrecInit( cvodeMem_, odeProblem_->numberOfEquations(), 1, 1 );
        AssertThrow( flag == CVSPILS_SUCCESS, ExcMessage("Failed to allocate BAND preconditioner") );
    }

    if ( odeProblem_->useJacobian() )
    {
        int flag;
        flag = CVSpilsSetJacTimesVecFn( cvodeMem_, ExplicitOdeJacobian );
        Assert( flag == CV_SUCCESS, ExcMessage("Call to SpilsSetJactimesVecFn failed"));
    }

}



// ----------------------------------------------------------------- CVodeTFQMR

CVodeTFQMR::CVodeTFQMR()
{
    solverName_ = "CVODE with TFQMR";
}



void
CVodeTFQMR::assignExplicitOde( 
        ExplicitOde& odeProblem,
        const realtype initialTime,
        const Vector<realtype> &initialState )
{
    CVodeSpils::assignExplicitOde( odeProblem, initialTime, initialState );

    if ( not precond_ )
    {
        int flag;
        flag = CVSptfqmr( cvodeMem_, PREC_NONE, krylovSubspaceDim_ );
        AssertThrow( flag == CV_SUCCESS, ExcCVSptfqmrError( flag ) );
    }
    else
    {
        int flag;
        flag = CVSptfqmr( cvodeMem_, PREC_RIGHT, krylovSubspaceDim_ );
        AssertThrow( flag == CV_SUCCESS, ExcCVSptfqmrError( flag ) );

        flag = CVBandPrecInit( cvodeMem_, odeProblem_->numberOfEquations(), 1, 1 );
        AssertThrow( flag == CVSPILS_SUCCESS, ExcMessage("Failed to allocate BAND preconditioner") );
    }

    if ( odeProblem_->useJacobian() )
    {
        int flag;
        flag = CVSpilsSetJacTimesVecFn( cvodeMem_, ExplicitOdeJacobian );
        Assert( flag == CV_SUCCESS, ExcMessage("Call to SpilsSetJactimesVecFn failed"));
    }

}
