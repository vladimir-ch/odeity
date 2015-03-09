#ifndef EXPLICIT_RUNGE_KUTTA_BASE_H
#define EXPLICIT_RUNGE_KUTTA_BASE_H

#include "RungeKuttaBase.h"
#include "IntegratorStats.h"
#include "../utils/SimpleArray.h"

#include <vector>
#include <fstream>

/* The method is based on the following scheme
 *  C(1)  |  A(1,1)   A(1,2)   ...  A(1,s)
 *  C(2)  |  A(2,1)   A(2,2)   ...  A(2,s)
 *   ..   |        ...
 *  C(s)  |  A(s,1)   A(s,1)   ...  A(s,s)
 *  -------------------------------------------------
 *        |  B'(1)    B'(2)    ...  B'(s)
 *  -------------------------------------------------
 *        |  B(1)     B(2)     ...  B(s)
 *
 * 	E(1) = B'(1) - B(1), ..., E(s) = B'(s) - B(s)
 */

class ExplicitRungeKuttaBase : public RungeKuttaBase
{
  public:

    void assignExplicitOde( ExplicitOde& odeProblem,
                            const realtype initialTime,
                            const Vector<realtype>& initialState );
    void printInfo() const;
    IntegratorStatsBase& stats();

  protected:

    int  stages_; // number of stages the method uses;
    bool fsal_;  // is "first-same-as-last"?

    SimpleArray<realtype>           a_;
    Vector<realtype>                b_;
    Vector<realtype>                c_;
    Vector<realtype>                e_;
    std::vector< Vector<realtype> > k_;

    realtype tooSmall_; // tooSmall_ - quantity used to determine when a step size is too small for
                        // the precision available

    RungeKuttaStats stats_;
    std::ofstream stepsizeHistoryFile_;


    ExplicitRungeKuttaBase( const int orderError, const int stages, const bool fsal );

    void startStep();
    void calculateSolutionPoint();
    void estimateError();
    void newStepsizeAfterRejectedStep();
    void newStepsizeAfterAcceptedStep();
    void completeStep();
    void prepareNextStep();

    void initializeIntegration() {}
    void estimateInitialStepsize();
    void updateHistory();

    void computeTooSmall();
};

#endif // EXPLICIT_RUNGE_KUTTA_BASE_H

