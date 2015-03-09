#include "OdeIntegratorBase.h"
#include "../odesystem/ExplicitOde.h"
#include "../utils/ParameterHandler.h"
#include "../utils/LogStream.h"

#include <limits>
#include <cmath>

OdeIntegratorBase::OdeIntegratorBase()
  :
    odeProblem_( 0 ),
    relTol_( 1.0e-3 ),
    absTol_( 1.0e-6 ),
    currentTime_( 0.0 ),
    oldTime_( 0.0 ),
    epsilon_( std::numeric_limits<realtype>::epsilon() ),
    sqrtEpsilon_( std::sqrt( std::numeric_limits<realtype>::epsilon() ) ),
    cubicEpsilon_( std::pow( std::numeric_limits<realtype>::epsilon(), 3.0 ) ),
    tiny_( std::sqrt( std::numeric_limits<realtype>::min() ) ),
    saveHistory_( false ),
    maxTol_( 0.1 ),
    minTol_( 10*std::numeric_limits<realtype>::epsilon() )
{ }



void OdeIntegratorBase::assignExplicitOde (
    ExplicitOde&            odeProblem,
    realtype                initialTime,
    const Vector<realtype> &initialState )
{
  odeProblem_ = &odeProblem;
  currentState_.reinit( odeProblem.numberOfEquations() );
  currentState_ = initialState;
  currentTime_ = initialTime;
}



void OdeIntegratorBase::printInfo( ) const
{
  using namespace std;

  logger << endl;
  logger << "ODE integrator information" << endl;
  logger << "--------------------------" << endl;
  logger << "*Solver name*:: " << solverName_ << endl;
  logger << "*Relative tolerance*:: " << relTol_ << endl;
  logger << "*Absolute tolerance*:: " << absTol_ << endl;
}



void
OdeIntegratorBase::setRelativeTolerance( realtype relTol )
{
  AssertThrow( minTol_ <= relTol && relTol <= maxTol_,
      ExcNotInAdmissibleRange( relTol, minTol_, maxTol_ )); 
  relTol_ = relTol;
}



realtype
OdeIntegratorBase::relativeTolerance() const
{
  return relTol_;
}



void
OdeIntegratorBase::setAbsoluteTolerance( realtype absTol )
{
  AssertThrow( 0.0 <= absTol && absTol <= maxTol_,
      ExcNotInAdmissibleRange( absTol, minTol_, maxTol_ )); 
  absTol_ = absTol;
}



realtype
OdeIntegratorBase::absoluteTolerance() const
{
  return absTol_;
}



void
OdeIntegratorBase::setSaveHistory( bool save )
{
  saveHistory_ = save;
}


bool
OdeIntegratorBase::saveHistory() const
{
    return saveHistory_;
}



void
OdeIntegratorBase::declareParameters( ParameterHandler & prm )
{
  prm.enter_subsection("ODE integrator");
      prm.declare_entry( "relative tolerance", "1.0e-08", Patterns::Double() );
      prm.declare_entry( "absolute tolerance", "1.0e-08", Patterns::Double() );
  prm.leave_subsection();
}



void
OdeIntegratorBase::getParameters( ParameterHandler & prm )
{
  prm.enter_subsection("ODE integrator");
      setRelativeTolerance( prm.get_double("relative tolerance") );
      setAbsoluteTolerance( prm.get_double("absolute tolerance") );
  prm.leave_subsection();
}

