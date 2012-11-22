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


#ifndef RUNGE_KUTTA_CHEBYSHEV_H
#define RUNGE_KUTTA_CHEBYSHEV_H

#include "RungeKuttaBase.h"
#include "IntegratorStats.h"

#include <fstream>

class ExplicitOde;

class RungeKuttaChebyshev : public RungeKuttaBase
{
  public:
    RungeKuttaChebyshev();

    void assignExplicitOde( ExplicitOde& odeProblem,
                            const realtype initialTime,
                            const Vector<realtype>& initialState );

    void printInfo() const;
    void setMaxSpecRadIterations( const int maxIterations );
    int getMaxSpecRadIterations() const;
    IntegratorStatsBase& stats();

    DeclException1( ExcMaxIterExceededSR,
        int,
        << "Maximum number of iterations (" << arg1 << ") exceeded during "
        "computation of the spectral radius of Jacobian " << std::endl );

  private:

    int maxStage_;
    int stages_;
    realtype spectralRadius_;
    bool computeSpectralRadius_;
    bool jacobianAtT_;
    int computeSRFlag_;
    realtype small_;
    int maxIterations_;

    Vector<realtype> fn_;
    Vector<realtype> fev_;
    Vector<realtype> temp1_;
    Vector<realtype> temp2_;
    Vector<realtype> yjm1_;
    Vector<realtype> yjm2_;
    Vector<realtype> eigenVector_;

    RungeKuttaChebyshevStats stats_;
    std::ofstream stepsizeHistoryFile_;
    std::ofstream stageHistoryFile_;

    void initializeIntegration();
    void estimateInitialStepsize();
    void estimateSpectralRadius();
    void startStep();
    void calculateSolutionPoint();
    void estimateError();
    void newStepsizeAfterRejectedStep();
    void handleRejectedStep();
    void newStepsizeAfterAcceptedStep();
    void completeStep();
    void prepareNextStep();
    void updateHistory();
};



inline
void
RungeKuttaChebyshev::setMaxSpecRadIterations( const int maxIterations )
{
  maxIterations_ = maxIterations_;
}



inline
int
RungeKuttaChebyshev::getMaxSpecRadIterations() const
{
  return maxIterations_;
}

OdeIntegratorBase * createRungeKuttaChebyshevSolver();

#endif // RUNGE_KUTTA_CHEBYSHEV_H

