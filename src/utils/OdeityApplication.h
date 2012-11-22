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

#ifndef ODEITYAPPLICATION_H
#define ODEITYAPPLICATION_H

#include "../functions/Function.h"
#include "../geometry/RectangularDomain.h"
#include "../geometry/RectangularGrid.h"
#include "../io/NetCDFWriter.h"
#include "../odesystem/MolOdeSystem.h"
#include "Factory.h"
#include "ParameterHandler.h"
#include "Vector.h"

#include <string>
#include <fstream>

class OdeIntegratorBase;

class OdeityApplication
{
    public:

        virtual ~OdeityApplication();

        int run();
        virtual void initialize();

    protected:

        OdeityApplication( int argc, char * argv[] );
        virtual MolOdeSystem<2>* createMolProblem() const = 0;
        virtual void writeGlobalAttributes();

        ParameterHandler prm;
        NetCDFWriter<2> * writer;
        MolOdeSystem<2>* molProblem;

    private:

        std::string configFileName_;
        std::string computationName_;
        std::string programName_;

        RectangularDomain<2> domain;
        RectangularGrid<2> grid;

        realtype finalTime;
        int numOutputPoints;
        bool precond, saveHistory, saveResults, useJacobian;

        Vector<realtype> initialState;
        Function<2> * initialCondition;
        Factory<Function<2>,std::string> initCondFactory;

        OdeIntegratorBase * solver;
        Factory<OdeIntegratorBase,std::string> solverFactory;

        std::ofstream logFile;
        std::ofstream historyFile_;

        void printHeader();

        void declareParameters();
        void getParameters();

        void registerInitialConditions();
        void registerOdeSolvers();

};


#endif // ODEITYAPPLICATION_H
