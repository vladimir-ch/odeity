#ifndef RUNGE_KUTTA_BASE_H
#define RUNGE_KUTTA_BASE_H

#include "OdeIntegratorBase.h"


class RungeKuttaBase : public OdeIntegratorBase
{
    public:

        void assignExplicitOde( ExplicitOde& odeProblem, realtype initialTime, const Vector<realtype>& initialState );
        void integrateTo( realtype tOut );
        void setStabilizationFactor( realtype stabilizationBeta );
        void setFafetyFactor( realtype value );
        realtype getFafetyFactor() const;
        virtual void printInfo() const;

        static void declareParameters( ParameterHandler & prm )
        {
            OdeIntegratorBase::declareParameters( prm );

            prm.enter_subsection( "ODE integrator" );
            prm.enter_subsection( "Runge-Kutta" );
                prm.declare_entry( "safety", "0.9", Patterns::Double() );
            prm.leave_subsection();
            prm.leave_subsection();
        }
        void getParameters( ParameterHandler & prm )
        {
            OdeIntegratorBase::getParameters( prm );

            prm.enter_subsection( "ODE integrator" );
            prm.enter_subsection( "Runge-Kutta" );
                safety_ = prm.get_double( "safety" );
            prm.leave_subsection();
            prm.leave_subsection();
        }

    protected:

        int orderError_;

        Vector<realtype> newState_;
        Vector<realtype> localError_;
        Vector<realtype> newLocalError_;
        Vector<realtype> weights_;
        realtype newLocalErrorNorm_;
        realtype oldLocalErrorNorm_;
        realtype oldStepsize_;

        realtype minStepsize_;

        realtype stabilizationAlpha_;
        realtype stabilizationBeta_;

        realtype stepsize_;
        realtype endTime_;
        realtype newTime_;
        bool last_;  // endTime can be reached with current stepsize?

        realtype safety_;
        realtype decreaseFactor_;
        realtype increaseFactor_;

        virtual void estimateInitialStepsize() = 0;
        virtual void initializeIntegration() = 0;
        virtual void startStep() = 0;
        virtual void estimateError() = 0;
        virtual void calculateSolutionPoint() = 0;
        virtual void newStepsizeAfterAcceptedStep() = 0;
        virtual void newStepsizeAfterRejectedStep() = 0;
        virtual void prepareNextStep() = 0;

        RungeKuttaBase( int orderError );

        virtual void handleRejectedStep();
        virtual void calculateWeights( const Vector<realtype>& y);
        virtual bool testAccuracy();
        virtual void completeStep();
        virtual void performIntegrationStep();
};



/* ----------------------------- Inline functions and templates ---------------- */

inline
void RungeKuttaBase::setFafetyFactor( realtype value )
{
    safety_ = value;
}



inline
realtype RungeKuttaBase::getFafetyFactor() const
{
    return safety_;
}


#endif // RUNGE_KUTTA_BASE_H

