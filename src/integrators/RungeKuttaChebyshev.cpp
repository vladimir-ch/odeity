#include "RungeKuttaChebyshev.h"
#include "../odesystem/ExplicitOde.h"
#include "../utils/LogStream.h"

#include <cmath>


RungeKuttaChebyshev::RungeKuttaChebyshev()
  :
    RungeKuttaBase( 2 ),
    spectralRadius_( 0.0 ),
    computeSpectralRadius_( true ),
    jacobianAtT_( false ),
    computeSRFlag_( 0 ),
    maxIterations_( 50 )
{
  solverName_ = "Stabilized recursive Runge-Kutta-Chebyshev";
  maxStage_ = (int)( std::sqrt( relTol_ / ( 10.0 * epsilon_) ) );
  maxStage_ = std::max( maxStage_, 2 );
}



IntegratorStatsBase&
RungeKuttaChebyshev::stats()
{
  return (IntegratorStatsBase&) stats_;
}



void RungeKuttaChebyshev::assignExplicitOde( ExplicitOde& odeProblem,
    const realtype initialTime,
    const Vector<realtype>& initialState )
{
  RungeKuttaBase::assignExplicitOde( odeProblem, initialTime, initialState );

  stats_.reset();

  fn_.reinit( odeProblem.numberOfEquations() );
  fev_.reinit( odeProblem.numberOfEquations() );
  temp1_.reinit( odeProblem.numberOfEquations() );
  temp2_.reinit( odeProblem.numberOfEquations() );
  yjm1_.reinit( odeProblem.numberOfEquations() );
  yjm2_.reinit( odeProblem.numberOfEquations() );
  eigenVector_.reinit( odeProblem.numberOfEquations() );

  odeProblem.rhs( currentTime_, currentState_, fn_ );
  stats_.incRhsEvaluations();

  computeSpectralRadius_ = true;
  jacobianAtT_ = false;
}


void RungeKuttaChebyshev::printInfo() const
{
  using namespace std;

  RungeKuttaBase::printInfo();

  logger << endl;
  logger << "Concrete solver information" << endl;
  logger << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  logger << "*Maximum number of stages*:: " << maxStage_ << endl;
  logger << "*Maximum number of iterations for spectral radius computation*:: " << maxIterations_ << endl;
  logger << "*Small number*:: " << small_ 
    << " (if spectral radius is less than this quantity, it does not restrict the time step size)" << endl << endl;
}



void RungeKuttaChebyshev::initializeIntegration()
{
  if ( computeSpectralRadius_ )
  {
    estimateSpectralRadius();
    jacobianAtT_ = true;
  }

  small_ = 1.0 / ( endTime_ - currentTime_ );
}



void RungeKuttaChebyshev::estimateInitialStepsize()
{
  realtype dt = endTime_ - currentTime_;

  minStepsize_ = 10.0 * epsilon_ * dt;
  stepsize_ = dt;
  if ( spectralRadius_ * stepsize_ > 1.0 )
    stepsize_ = 1.0 / spectralRadius_;

  stepsize_ = std::max( stepsize_, minStepsize_ );
  temp1_.equ( 1.0, currentState_, stepsize_, fn_ );
  odeProblem_->rhs( currentTime_ + stepsize_, temp1_, temp2_ );
  stats_.incRhsEvaluations();
  temp2_.add( -1.0, fn_ );

  realtype estimate = stepsize_ * temp2_.rms_norm( weights_ );

  if ( 0.1 * stepsize_ < dt * std::sqrt( estimate ) )
    stepsize_ = std::max( 0.1 * stepsize_ / std::sqrt( estimate ) , minStepsize_ );
  else
    stepsize_ = dt;
}



void RungeKuttaChebyshev::estimateSpectralRadius()
{
  if ( stats_.getAcceptedSteps() == 0 )
    eigenVector_ = fn_;

  realtype yNorm = currentState_.l2_norm();
  realtype evNorm = eigenVector_.l2_norm();
  realtype dyNorm;

  if ( yNorm != 0.0 && evNorm != 0.0 )
  {
    dyNorm = yNorm * sqrtEpsilon_;
    eigenVector_.sadd( dyNorm / evNorm, 1.0, currentState_ );
  }
  else if ( yNorm != 0.0)
  {
    dyNorm = yNorm * sqrtEpsilon_;
    eigenVector_.equ( 1.0, currentState_, sqrtEpsilon_, currentState_ );
  }
  else if ( evNorm != 0.0 )
  {
    dyNorm = epsilon_;
    eigenVector_.scale( dyNorm / evNorm );
  }
  else
  {
    dyNorm = epsilon_;
    eigenVector_ = dyNorm;
  }

  //--------------------------------------------
  //  Now iterate with a nonlinear power method.
  //--------------------------------------------
  realtype sigma = 0.0;
  realtype sigmal;
  realtype dfNorm = 0.0;

  int iter;
  for ( iter = 1; iter <= maxIterations_; iter++ )
  {
    odeProblem_->rhs( currentTime_, eigenVector_, fev_ );
    stats_.incRhsEvaluations();

    dfNorm = 0.0;
    fev_.add( -1.0, fn_ );
    dfNorm = fev_.l2_norm();

    sigmal = sigma;
    sigma = dfNorm/dyNorm;

    //----------------------------------------------------------
    //  sprad is a little bigger than the estimate sigma of the 
    //  spectral radius, so is more likely to be an upper bound.
    //----------------------------------------------------------
    spectralRadius_ = 1.2 * sigma;
    if ( iter >= 2 && std::fabs( sigma - sigmal) <= std::max( sigma, small_) * 0.01 )
    {
      eigenVector_.add( -1.0, currentState_ );
      ++stats_.nspraditers_;
      break;
      // return;
    }
    //--------------------------------------
    //  The next v(*) is the change in f
    //  scaled so that norm(v - yn) = dynrm.
    //--------------------------------------
    if ( dfNorm != 0.0 )
      eigenVector_.equ( 1.0, currentState_, dyNorm/dfNorm, fev_);
    else
    {
      //-------------------------------------------------------
      //  The new v(*) degenerated to yn(*)--"randomly" perturb 
      //  current approximation to the eigenvector by changing 
      //  the sign of one component.
      //-------------------------------------------------------
      int index = 1 + ( iter % odeProblem_->numberOfEquations() );
      eigenVector_(index) = currentState_(index)
        - (eigenVector_(index) - currentState_(index));
    }
  }

  AssertThrow( iter <= maxIterations_, ExcMaxIterExceededSR( maxIterations_ ) );
}



void RungeKuttaChebyshev::startStep()
{
  if ( computeSpectralRadius_ )
  {
    estimateSpectralRadius();
    jacobianAtT_ = true;
  }

  last_ = false;
  if ( 1.2 * stepsize_ > endTime_ - currentTime_ )
  {
    stepsize_ = endTime_ - currentTime_;
    last_ = true;
  }

  stages_ = 1 + (int)( std::sqrt( 1.54 * stepsize_ * spectralRadius_ + 1.0 ) );
  if ( stages_ > maxStage_ )
  {
    stages_ = maxStage_;
    stepsize_ = ( stages_*stages_ - 1 )/( 1.54 * spectralRadius_ );
    last_ = false;
  }
  // stats_.updateStage( currentTime_, stages_ );
}



void RungeKuttaChebyshev::calculateSolutionPoint()
{
  realtype temp1, temp2;
  realtype w0, w1;
  realtype ajm1, bj, bjm1, bjm2;
  realtype arg;
  realtype mus, mu, nu;
  realtype thj, thjm1, thjm2;
  realtype zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2;
  realtype zj, dzj, d2zj;

  // initialization
  w0 = 1.0 + 2.0/(13.0 * stages_ * stages_ );
  temp1 = w0 * w0 - 1.0;
  temp2 = sqrt( temp1 );
  arg = stages_ * std::log( w0 + temp2 );
  w1 = std::sinh( arg ) * temp1 / ( std::cosh( arg )*stages_*temp2 - w0*std::sinh(arg) );
  bjm1 = 1.0 / (4.0 * w0*w0 );
  bjm2 = bjm1;

  // compute first stage
  yjm2_ = currentState_;
  mus = w1 * bjm1;
  yjm1_.equ( 1.0, currentState_, stepsize_ * mus, fn_ );
  thjm2  = 0.0;
  thjm1  = mus;
  zjm1   = w0;
  zjm2   = 1.0;
  dzjm1  = 1.0;
  dzjm2  = 0.0;
  d2zjm1 = 0.0;
  d2zjm2 = 0.0;

  // compute other stages s = 2,..., stages_
  for ( int s = 2; s <= stages_; s++)
  {
    zj   =   2.0*w0*zjm1 - zjm2;
    dzj  =   2.0*w0*dzjm1 - dzjm2 + 2.0*zjm1;
    d2zj =   2.0*w0*d2zjm1 - d2zjm2 + 4.0*dzjm1;
    bj   =   d2zj/(dzj*dzj);
    ajm1 =   1.0 - zjm1*bjm1;
    mu   =   2.0*w0*bj/bjm1;
    nu   = - bj/bjm2;
    mus  =   mu*w1/w0;

    odeProblem_->rhs( currentTime_ + stepsize_ * thjm1, yjm1_, newState_ );
    stats_.incRhsEvaluations();

    newState_.sadd( stepsize_ * mus, mu, yjm1_, nu, yjm2_, 1.0 - mu - nu, currentState_ );
    newState_.add( - stepsize_ * mus * ajm1, fn_ );

    thj = mu*thjm1 + nu*thjm2 + mus*(1.0 - ajm1);

    if ( s < stages_ )
    {
      swap( yjm2_, yjm1_ );
      yjm1_ = newState_;

      thjm2  = thjm1;
      thjm1  = thj;
      bjm2   = bjm1;
      bjm1   = bj;
      zjm2   = zjm1;
      zjm1   = zj;
      dzjm2  = dzjm1;
      dzjm1  = dzj;
      d2zjm2 = d2zjm1;
      d2zjm1 = d2zj;
    }
  }

  newTime_ = currentTime_ + stepsize_;

  odeProblem_->rhs( newTime_, newState_, temp1_ );
  stats_.incRhsEvaluations();
}



void RungeKuttaChebyshev::estimateError()
{
  temp2_.maxabs( currentState_, newState_ );
  calculateWeights( temp2_ );

  newLocalError_.equ( 0.8, currentState_, -0.8, newState_);
  newLocalError_.add( 0.4*stepsize_, fn_, -0.4*stepsize_, temp1_ );

  newLocalErrorNorm_ = newLocalError_.rms_norm( weights_ );       
}



void RungeKuttaChebyshev::newStepsizeAfterRejectedStep()
{
  stepsize_ = safety_ * stepsize_ / std::pow( newLocalErrorNorm_, 1.0/3.0 );
}



void RungeKuttaChebyshev::handleRejectedStep()
{
  RungeKuttaBase::handleRejectedStep();

  stats_.incRejectedSteps();
  computeSpectralRadius_ = not jacobianAtT_; // if we have jacobian, do not compute it
}



void RungeKuttaChebyshev::newStepsizeAfterAcceptedStep()
{
  realtype temp1, temp2, factor = 10.0;

  if ( stats_.getAcceptedSteps() == 0 )
  {
    temp1 = std::pow( newLocalErrorNorm_, 1.0/3.0 );
    if ( 0.8 < factor * temp1 )
      factor = 0.8 / temp1;
  }
  else
  {
    temp1 = 0.8 * stepsize_ * std::pow( oldLocalErrorNorm_, 1.0/3.0 );
    temp2 = oldStepsize_ * std::pow( newLocalErrorNorm_, 2.0/3.0 );
    if ( temp1 < factor * temp2 )
      factor = temp1 / temp2;
  }
  stepsize_ *= std::max( 0.1, factor);
  stepsize_ = std::max( minStepsize_, stepsize_ );
}



void RungeKuttaChebyshev::completeStep()
{
  stats_.incAcceptedSteps();
  // stats_.updateTimeStepsize( currentTime_ , stepsize_ );

  RungeKuttaBase::completeStep();
}



void RungeKuttaChebyshev::prepareNextStep()
{
  RungeKuttaBase::prepareNextStep();

  jacobianAtT_ = false; // we assume that jacobian is not constant
  computeSpectralRadius_ = false;

  computeSRFlag_ = ( computeSRFlag_ + 1 ) % 25;
  if ( computeSRFlag_ == 0 )
    computeSpectralRadius_ = not jacobianAtT_;

  swap( fn_, temp1_ );
}



void
RungeKuttaChebyshev::updateHistory()
{
  stats_.updateHistory( oldTime_, oldStepsize_, stages_ );
}



OdeIntegratorBase * createRungeKuttaChebyshevSolver()
{
    return new RungeKuttaChebyshev();
}

