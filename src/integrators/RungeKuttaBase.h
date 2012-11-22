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

