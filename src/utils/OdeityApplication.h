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
