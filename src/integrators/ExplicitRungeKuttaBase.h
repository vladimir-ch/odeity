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

