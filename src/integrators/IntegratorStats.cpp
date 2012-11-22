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


#include "IntegratorStats.h"
#include "../utils/LogStream.h"
#include "../utils/Exceptions.h"

#include <cvode/cvode.h>                  /* main integrator header file */
#include <cvode/cvode_spgmr.h>            /* prototypes & constants for CVSPGMR solver */
#include <cvode/cvode_spbcgs.h>           /* prototypes & constants for CVSPBCG solver */
#include <cvode/cvode_sptfqmr.h>          /* prototypes & constants for CVSTFQMR solver */
#include <cvode/cvode_band.h>
#include <cvode/cvode_bandpre.h>

#include <fstream>
#include <iomanip>

/* --------------------------------------------------------------------------*/

IntegratorStatsBase::IntegratorStatsBase()
    :
        historyFile_( 0 )
{
    reset();
}



void
IntegratorStatsBase::setHistoryFile( std::ofstream & historyFile )
{
    historyFile_ = &historyFile;
    historyFile_->width( 18 );
    historyFile_->precision( 12 );
    historyFile_->setf( std::ios::scientific, std::ios::floatfield );
}



void
IntegratorStatsBase::updateTimeStepsize( realtype time, realtype stepsize )
{
    Assert( historyFile_ != 0, ExcNotInitialized() );

    if ( stepsize > maxTimeStepsize_ )
        maxTimeStepsize_ = stepsize;
    if ( stepsize < minTimeStepsize_ )
        minTimeStepsize_ = stepsize;

    *historyFile_ << std::endl << time << " " << stepsize;
}



void
IntegratorStatsBase::reset()
{
    maxTimeStepsize_ = 0.0;
    minTimeStepsize_ = 1.0;
}



/* --------------------------------------------------------------------------*/

RungeKuttaStats::RungeKuttaStats()
{
    reset();
}



void RungeKuttaStats::reset()
{
    IntegratorStatsBase::reset();

    rhsEvaluations_ = 0;
    acceptedSteps_ = 0;
    rejectedSteps_ = 0;
}



void
RungeKuttaStats::updateHistory( realtype time, realtype stepsize )
{
    IntegratorStatsBase::updateTimeStepsize( time, stepsize );
}



void RungeKuttaStats::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Runge-Kutta statistics" << endl;
    logger << "----------------------" << endl;
    logger << endl;
    logger << "[frame=\"topbot\",grid=\"rows\"]" << endl;
    logger << "`-----------------'--------------" << endl;
    logger << "RHS evaluations   " << setw( 15 ) << rhsEvaluations_ << endl;
    logger << "Accepted steps    " << setw( 15 ) << acceptedSteps_ << endl;
    logger << "Rejected steps    " << setw( 15 ) << rejectedSteps_ << endl;
    logger << "Maximum stepsize  " << setw( 15 ) << setprecision( 8 ) << maxTimeStepsize_ << endl;
    logger << "Minimum stepsize  " << setw( 15 ) << setprecision( 8 ) << minTimeStepsize_ << endl;
    logger << "---------------------------------" << endl;
}


/* --------------------------------------------------------------------------*/

RungeKuttaChebyshevStats::RungeKuttaChebyshevStats()
{
    reset();
}



void RungeKuttaChebyshevStats::reset()
{
    RungeKuttaStats::reset();

    nspraditers_ = 0;
    maxStage_ = 0;
}



void
RungeKuttaChebyshevStats::updateHistory( realtype time, realtype stepsize, int stage )
{
    Assert( historyFile_ != 0, ExcNotInitialized() );

    if ( stage > maxStage_ )
        maxStage_ = stage;

    RungeKuttaStats::updateHistory( time, stepsize );
    *historyFile_ << " " << stage;
}



void RungeKuttaChebyshevStats::printInfo() const
{
    using namespace std;

    RungeKuttaStats::printInfo();
    logger << "Spect. rad. iters " << setw( 15 ) << nspraditers_ << endl;
    if ( maxStage_ != -1 )
        logger << "Maximum stage used" << setw( 15 ) << maxStage_ << endl;
    logger << "---------------------------------" << endl;
}



// ----------------------------------------------------------------- CVodeStats
CVodeStats::CVodeStats()
{
    reset();
}



void CVodeStats::reset()
{
    IntegratorStatsBase::reset();

    nsteps = 0;
    nfevals = 0;
    nlinsetups = 0;
    netfails = 0;

    nniters = 0;
    nncfails = 0;

    maxOrder_ = 0;
}



void
CVodeStats::updateHistory( realtype time, realtype stepsize, int order )
{
    Assert( historyFile_ != 0, ExcNotInitialized() );

    updateTimeStepsize( time, stepsize );

    if ( order > maxOrder_ )
        maxOrder_ = order;

    *historyFile_ << " " << order;
}



void
CVodeStats::update( void * cvodeMem )
{
    int dummyi;
    realtype dummyf;

    int flag = CVodeGetIntegratorStats( cvodeMem,
            &nsteps, &nfevals, &nlinsetups, &netfails,
            &dummyi, &dummyi, &dummyf, &dummyf, &dummyf, &dummyf);
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVodeGetIntegratorStats failed"));

    flag = CVodeGetNonlinSolvStats( cvodeMem, &nniters, &nncfails);
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVodeGetNonlinSolvStats failed"));
}



void
CVodeStats::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "CVode statistics" << endl;
    logger << "----------------" << endl;
    logger << endl;
    logger << "[frame=\"topbot\",grid=\"rows\"]" << endl;
    logger << "`-------------------------------'--------------" << endl;
    logger << "RHS evaluations                 " << setw( 15 ) << nfevals << endl;
    logger << "Accepted steps                  " << setw( 15 ) << nsteps << endl;
    logger << "Error test failures             " << setw( 15 ) << netfails << endl;
    logger << "Linear solver setups            " << setw( 15 ) << nlinsetups << endl;
    logger << "Nr. of nonlinear iterations     " << setw( 15 ) << nniters << endl;
    logger << "Nr. of nonlinear conv. failures " << setw( 15 ) << nncfails << endl;
    if ( maxTimeStepsize_ > 0.0 )
        logger << "Maximum stepsize                " << setw( 15 ) << setprecision( 8 ) << maxTimeStepsize_ << endl;
    if ( minTimeStepsize_ < 1.0 )
        logger << "Minimum stepsize                " << setw( 15 ) << setprecision( 8 ) << minTimeStepsize_ << endl;
    if ( maxOrder_ > 0 )
        logger << "Maximum order                   " << setw( 15 ) << setprecision( 8 ) << maxOrder_ << endl;
    logger << "-----------------------------------------------" << endl;
}


// ------------------------------------------------------------ CVodeSpilsStats
void
CVodeSpilsStats::update( void * cvodeMem )
{
    CVodeStats::update( cvodeMem );

    int flag = CVSpilsGetNumLinIters( cvodeMem, &nliters );
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVSpilsGetNumLinIters failed"));

    flag = CVSpilsGetNumConvFails( cvodeMem, &nlcfails);
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVSpilsGetNumConvFails failed"));

    flag = CVSpilsGetNumPrecEvals( cvodeMem, &npevals); 
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVSpilsGetNumPrecEvals failed"));

    flag = CVSpilsGetNumPrecSolves( cvodeMem, &npsolves); 
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVSpilsGetNumPrecSolves failed"));

    flag = CVSpilsGetNumJtimesEvals(cvodeMem, &njvevals);
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVSpilsGetNumJtimesEvals failed"));

    flag = CVSpilsGetNumRhsEvals( cvodeMem, &nfevalsLS);
    Assert( flag == CV_SUCCESS, ExcMessage("Call to CVSpilsGetNumRhsEvals failed"));
}



void
CVodeSpilsStats::update( void * cvodeMem, void * bpData )
{
    update( cvodeMem );
    if ( bpData )
    {
        int flag = CVBandPrecGetNumRhsEvals( bpData, &nfevalsBP );
        Assert( flag == CV_SUCCESS, ExcMessage("Call to CVBandPreGetNumRhsEvals failed"));
    }
}



void
CVodeSpilsStats::printInfo() const
{
    using namespace std;

    CVodeStats::printInfo();
    logger << "Linear iterations               " << setw( 15 ) << nliters << endl;
    logger << "Linear convergence failures     " << setw( 15 ) << nlcfails << endl;
    logger << "Preconditioner evaluations      " << setw( 15 ) << npevals << endl;
    logger << "Preconditioner solves           " << setw( 15 ) << npsolves << endl;
    logger << "Jacobian-vector evaluations     " << setw( 15 ) << njvevals << endl;
    logger << "RHS evals for FD Jac-vec prod.  " << setw( 15 ) << nfevalsLS << endl;
    logger << "RHS calls in BANDPRE            " << setw( 15 ) << nfevalsBP << endl;
    logger << "-----------------------------------------------" << endl;
    logger << "Total number of RHS calls       " << setw( 15 ) << nfevalsLS+nfevals+nfevalsBP << endl;
    logger << "-----------------------------------------------" << endl;
}

