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


#include "RungeKuttaBase.h"
#include "../odesystem/ExplicitOde.h"
#include "../utils/LogStream.h"

#include <limits>

RungeKuttaBase::RungeKuttaBase( int orderError )
    :
        orderError_( orderError ),
        newLocalErrorNorm_( 0.0 ),
        oldLocalErrorNorm_( 0.0 ),
        minStepsize_( 10*std::numeric_limits<realtype>::epsilon() ),
        stepsize_( 0.0 ),
        endTime_( 0.0 ),
        newTime_( 0.0 ),
        last_( false ),
        safety_( 0.9 ),
        decreaseFactor_( 1.0e-1 ),
        increaseFactor_( 1.0e1 )
        //    maxErrorTestFailures_( 10 ),
{
    setStabilizationFactor( 0.1 );
}



void RungeKuttaBase::printInfo() const
{
    using namespace std;

    OdeIntegratorBase::printInfo();

    logger << "*Order of error*:: " << orderError_ << endl;
    logger << "*Stabilization factor alpha*:: " << stabilizationAlpha_ << endl;
    logger << "*Stabilization factor beta*:: " << stabilizationBeta_ << endl;
    logger << "*Safety factor*:: " << safety_ << endl;
    logger << "*Decrease factor*:: " << decreaseFactor_ << endl;
    logger << "*Increase factor*:: " << increaseFactor_ << endl;
}



void RungeKuttaBase::assignExplicitOde (
        ExplicitOde& odeProblem,
        realtype initialTime,
        const Vector<realtype>& initialState )
{
    OdeIntegratorBase::assignExplicitOde( odeProblem, initialTime, initialState );

    int neq = odeProblem.numberOfEquations();

    newState_.reinit( neq );
    localError_.reinit( neq );
    newLocalError_.reinit( neq );
    weights_.reinit( neq );
    calculateWeights( currentState_ );

    setStabilizationFactor( 0.1 );
}



void RungeKuttaBase::integrateTo( realtype tOut )
{
    Assert( odeProblem_ != 0, ExcNotInitialized() );

    AssertThrow( tOut - currentTime_ > std::numeric_limits<realtype>::epsilon(),
            ExcFinalTimeRange( tOut, currentTime_ ) );

    endTime_ = tOut;

    initializeIntegration();

    estimateInitialStepsize();

    AssertThrow( endTime_ - currentTime_ >= minStepsize_,
            ExcEndTimeTooClose( currentTime_, endTime_, minStepsize_ ) );

    while ( currentTime_ != endTime_ )
    {
        performIntegrationStep();
        if ( saveHistory() )
            updateHistory();
    }
}    



void RungeKuttaBase::setStabilizationFactor( realtype stabilizationBeta )
{
    AssertThrow( 0.0 <= stabilizationBeta && stabilizationBeta <= 0.2,
            ExcNotInAdmissibleRange( stabilizationBeta, 0.0, 0.2 ) );
    stabilizationBeta_ = stabilizationBeta;
    stabilizationAlpha_ = 1.0 / (realtype) ( orderError_ + 1 ) - 7.5e-1 * stabilizationBeta;
}



void RungeKuttaBase::calculateWeights( const Vector<realtype>& y )
{
    weights_.invabsequ( relTol_, y, absTol_ );
}



void RungeKuttaBase::handleRejectedStep()
{
    newStepsizeAfterRejectedStep();

    AssertThrow( stepsize_ >= minStepsize_, ExcTimestepTooSmall( stepsize_, minStepsize_ ) );
}



void RungeKuttaBase::completeStep()
{
    oldStepsize_ = stepsize_;
    oldTime_ = currentTime_;

    if ( last_ )
        currentTime_ = endTime_;
    else
        currentTime_ = newTime_;
}



void RungeKuttaBase::prepareNextStep()
{
    swap( currentState_, newState_ );
    swap( localError_, newLocalError_ );
    oldLocalErrorNorm_ = newLocalErrorNorm_;
}



void RungeKuttaBase::performIntegrationStep()
{
    while ( true )
    {
        startStep();
        calculateSolutionPoint();
        estimateError();

        if ( not testAccuracy() ) // failed step
        {
            handleRejectedStep();
        }
        else // successfull step
        {
            completeStep();
            newStepsizeAfterAcceptedStep();
            prepareNextStep();

            return;
        }
    }
}



inline
bool RungeKuttaBase::testAccuracy()
{
    return ( newLocalErrorNorm_ <= 1.0 ); // return true if accuracy is good
}

