#ifndef ODE_INTEGRATOR_BASE_H
#define ODE_INTEGRATOR_BASE_H

#include "../utils/Exceptions.h"
#include "../utils/Vector.h"
#include "../utils/ParameterHandler.h"

#include <sundials/sundials_types.h>

#include <iosfwd>
#include <string>

class ParameterHandler;
class IntegratorStatsBase;
class ExplicitOde;

class OdeIntegratorBase
{
    public:

        OdeIntegratorBase();
        virtual ~OdeIntegratorBase() {}

        virtual void assignExplicitOde(
                ExplicitOde& odeProblem,
                const realtype initialTime,
                const Vector<realtype> &initialState );

        virtual void integrateTo( const realtype tOut ) = 0;
        virtual IntegratorStatsBase& stats() = 0;
        virtual void updateHistory() = 0;

        void integrateForwardDt( const realtype dt );

        void setRelativeTolerance( realtype relTol );
        realtype relativeTolerance() const;
        void setAbsoluteTolerance( realtype absTol );
        realtype absoluteTolerance() const;

        realtype currentTime() const;
        const Vector<realtype>& currentState() const;
        const ExplicitOde * odeProblem() const;

        void setSaveHistory( bool save );
        bool saveHistory() const;
        
        virtual void printInfo() const;
        virtual void getParameters( ParameterHandler &prm );
        static void declareParameters( ParameterHandler & prm );

        DeclException3(ExcNotInAdmissibleRange,
                realtype, realtype, realtype,
                << "The real number " << arg1 << " has to be in the range [" << arg2 << "," << arg3 << "]" );

        DeclException2(ExcFinalTimeRange,
                realtype,realtype,
                << "Final integration time " << arg1 << " is lower than current time " << arg2
                << " or it is not sufficiently apart" );

        DeclException3(ExcEndTimeTooClose,
                realtype,realtype,realtype,
                << "Start and end times are not clearly distinguishable for the method\n \
                and the precision of the computer being used (time_end - time_start = "
                << arg2 - arg1 << "),\n should be at least " << arg3);

        DeclException1(ExcMaxStepsExceeded,
                int,
                << "Maximum number of steps exceeded during integration ( maxSteps = " << arg1 << " )");

        DeclException2(ExcTimestepTooSmall,
                realtype, realtype,
                << "Too small time stepsize for current precision would be needed to achieve "
                "requiered accuracy ( time stepsize = " << arg1 << ", minimum stepsize = " << arg2 << ")");


    protected:

        ExplicitOde     *odeProblem_;
        realtype         relTol_;
        realtype         absTol_;
        realtype         currentTime_;
        Vector<realtype> currentState_;
        realtype         oldTime_;

        const realtype epsilon_;
        const realtype sqrtEpsilon_;
        const realtype cubicEpsilon_;
        const realtype tiny_;

        std::string solverName_;

    private:

        bool saveHistory_;
        const realtype maxTol_;
        const realtype minTol_;
};


inline
void
OdeIntegratorBase::integrateForwardDt( const realtype dt )
{
    integrateTo( dt + currentTime() );
}


inline
realtype
OdeIntegratorBase::currentTime() const
{
    return currentTime_;
}


inline
const Vector<realtype>&
OdeIntegratorBase::currentState() const
{
    return currentState_;
}


inline
const ExplicitOde *
OdeIntegratorBase::odeProblem() const
{
    return odeProblem_;
}


#endif // ODE_INTEGRATOR_BASE_H
