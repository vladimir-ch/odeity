#include "ExplicitRungeKuttaBase.h"
#include "../odesystem/ExplicitOde.h"
#include "../utils/LogStream.h"

#include <cmath>

ExplicitRungeKuttaBase::ExplicitRungeKuttaBase(
        const int  orderError,
        const int  stages,
        const bool fsal )
:
    RungeKuttaBase( orderError ),
    stages_( stages ),
    fsal_( fsal ),
    a_( stages, stages ),
    b_( stages ),
    c_( stages ),
    e_( stages ),
    k_( stages )
{}



void
ExplicitRungeKuttaBase::printInfo() const
{
    using namespace std;

    RungeKuttaBase::printInfo();

    logger << "*Number of stages*:: " << stages_ << endl;
    if ( fsal_ )
        logger << "*FSAL*:: true" << endl;
    else
        logger << "*FSAL*:: false" << endl;
    logger << "*Small number*:: " << tooSmall_ 
        << " (quantity used to determine when a step size is too small for the precision available)" << endl << endl;
}



void
ExplicitRungeKuttaBase::assignExplicitOde(
        ExplicitOde& odeProblem,
        const realtype initialTime,
        const Vector<realtype>& initialState )
{
    RungeKuttaBase::assignExplicitOde( odeProblem, initialTime, initialState );

    stats_.reset();

    int neq = odeProblem.numberOfEquations();

    for ( int i = 0; i < stages_; i++ )
        k_[i].reinit( neq );

    computeTooSmall();

    odeProblem_->rhs( currentTime_, currentState_, k_[0] ); 
    stats_.incRhsEvaluations();
}



IntegratorStatsBase&
ExplicitRungeKuttaBase::stats()
{
    return (IntegratorStatsBase&) stats_;
}



void
ExplicitRungeKuttaBase::computeTooSmall()
{
    realtype cDiff = 1.0;
    realtype diff;
    for ( int i = 0; i < stages_ - 1; i++)
        for ( int j = i+1; j < stages_; j++)
        {
            diff = std::fabs( c_(i) - c_(j) );
            if ( diff > 0 )
                cDiff = std::min( cDiff, diff );
        }
    tooSmall_ = 10 * epsilon_ / cDiff;
}



void
ExplicitRungeKuttaBase::estimateInitialStepsize( )
{
    minStepsize_ = std::max( tiny_, tooSmall_ * ( endTime_ - currentTime_ ) );

    stepsize_ = endTime_ - currentTime_;

    for ( unsigned int i = 0; i < currentState_.size() ; i++ )
    {
        realtype tol = relTol_ * std::fabs( currentState_(i) ) + absTol_;
        realtype ypk = std::fabs( k_[0](i) );
        if ( ypk * std::pow( stepsize_, 5.0 ) > tol )
            stepsize_ = std::pow( tol/ypk, 0.2 );
    }

    stepsize_ = std::max( minStepsize_, std::min( stepsize_, endTime_ - currentTime_) );
}



void
ExplicitRungeKuttaBase::startStep()
{
    last_ = false;

    realtype dt = endTime_ - currentTime_;
    if ( dt < stepsize_ ) // if in the current stepsize we could reach the endTime, set stepsize appropriately
    {
        stepsize_ = dt;
        last_ = true;
    }
    else if ( dt < 2.0*stepsize_ ) // if in two current stepsizes we could reach the endTime, halve the stepsize
    {
        stepsize_ *= 0.5;
    }
}



void
ExplicitRungeKuttaBase::calculateSolutionPoint()
{
    for ( int s = 1; s < stages_; s++ )
    {
        newState_ = currentState_;
        for ( int i = 0; i < s; i++ )
            if ( a_(s,i) != 0.0 )
                newState_.add( stepsize_ * a_(s,i), k_[i] );

        odeProblem_->rhs( currentTime_ + stepsize_*c_(s), newState_, k_[s] );
        stats_.incRhsEvaluations();
    }

    newTime_ = currentTime_ + stepsize_;

    if ( not fsal_ )
    {
        newState_ = currentState_;
        for ( int s = 0; s < stages_; s++ )
            if ( b_(s) != 0.0 )
                newState_.add( stepsize_ * b_(s), k_[s] );
    }
}



void
ExplicitRungeKuttaBase::estimateError()
{
    newLocalError_ = 0.0;
    for ( int s = 0; s < stages_; s++)
        if ( e_(s) != 0.0 )
            newLocalError_.add( e_(s), k_[s] );
    newLocalError_ *= stepsize_;
    newLocalErrorNorm_ = newLocalError_.rms_norm( weights_ );
}



void
ExplicitRungeKuttaBase::newStepsizeAfterAcceptedStep()
{
    realtype factor;

    if ( stats_.getAcceptedSteps() >= 1 )
        factor = std::pow( oldLocalErrorNorm_, stabilizationBeta_) / std::pow( newLocalErrorNorm_, stabilizationAlpha_);
    else
        factor = std::pow( 1.0 / newLocalErrorNorm_, 1 / (realtype)(orderError_ + 1) );

    factor = std::min( increaseFactor_, std::max( decreaseFactor_, safety_ * factor ) );

    stepsize_ *= factor;

}



void
ExplicitRungeKuttaBase::newStepsizeAfterRejectedStep()
{
    realtype factor;

    factor = std::pow( 1.0 / newLocalErrorNorm_, 1.0 / (realtype)( orderError_ + 1.0 ) );
    factor = std::min( increaseFactor_, std::max( decreaseFactor_, safety_ * factor ) );
    stepsize_ *= factor;

    stats_.incRejectedSteps();
}



void
ExplicitRungeKuttaBase::completeStep()
{
    stats_.incAcceptedSteps();
    // stats_.updateTimeStepsize( currentTime_ , stepsize_ );

    RungeKuttaBase::completeStep();

    calculateWeights( currentState_ );
}



void
ExplicitRungeKuttaBase::prepareNextStep()
{
    RungeKuttaBase::prepareNextStep();

    if ( fsal_ )
        swap( k_[0], k_[stages_-1] );
    else
    {
        odeProblem_->rhs( currentTime_, currentState_, k_[0] );
        stats_.incRhsEvaluations();
    }
}



void
ExplicitRungeKuttaBase::updateHistory()
{
    stats_.updateHistory( oldTime_, oldStepsize_ );
}

