#include "../functions/WavyCircleFunction.h"
#include "../functions/ZeroFunction.h"
#include "../functions/RandomFunction.h"
#include "../functions/SincFunction.h"
#include "../functions/MathevalFunction.h"
#include "../functions/TwoScaleSineFunction.h"
#include "../functions/RectangleFunction.h"
#include "../functions/TwoRectanglesFunction.h"
#include "../functions/ZigZagFunction.h"
#include "../geometry/Point.h"
#include "../integrators/DormandPrince45.h"
#include "../integrators/RungeKuttaMerson45.h"
#include "../integrators/RungeKutta23.h"
#include "../integrators/RungeKuttaChebyshev.h"
#include "../integrators/CVode.h"
#include "../integrators/IntegratorStats.h"
#include "../io/NetCDFWriter.h"
#include "../odesystem/AllenCahnEquation.h"
#include "../odesystem/CahnHilliardEquation.h"
#include "../odesystem/DegenerateCahnHilliardEquation.h"
#include "../odesystem/LoretiMarchEquation.h"

#include "OdeityApplication.h"
#include "JobIdentifier.h"
#include "ProgressDisplay.h"
#include "Timer.h"
#include "LogStream.h"

#include <iostream>
#include <iomanip>
#include <fstream>

OdeityApplication::OdeityApplication( int argc, char ** argv )
    :
        writer( 0 ),
        molProblem( 0 ),
        finalTime( 1.0 ),
        numOutputPoints( 10 ),
        precond( false ),
        saveHistory( false ),
        saveResults( true ),
        useJacobian( false ),
        initialCondition( 0 ),
        solver( 0 )
{
    // parse command line parameters
    // TODO: better message and usage printing
    AssertThrow( argc == 2, ExcMessage("Parameter file name expected on command line") );

    configFileName_ = argv[1];
    programName_ = argv[0];

    // register creator functions in the factories
    registerInitialConditions();
    registerOdeSolvers();

    declareParameters();
}



OdeityApplication::~OdeityApplication()
{
    delete initialCondition;
    delete solver;
    delete molProblem;
    delete writer;
}



void
OdeityApplication::initialize()
{
    getParameters();

    writer = new NetCDFWriter<2>( computationName_, "data.nc" );

    // make a copy of the parameters for later reference in the output directory
    {
        std::ofstream prmFile( std::string(writer->outputPath() + configFileName_).c_str() );
        prm.print_parameters( prmFile, ParameterHandler::Text );
    }

    // open log file and attach it to the log stream
    {
        std::string reportFileName = writer->outputPath() + "report.log";
        logFile.open( reportFileName.c_str() );
        logger.attach( logFile );
    }

    // initial condition
    initialState.reinit( grid.numberOfNodes() );
    for( unsigned int i = 0; i < grid.dimension( xDim ); i++ )
        for( unsigned int j = 0; j < grid.dimension( yDim ); j++ )
            initialState( grid.nodeIndex(i,j) ) = (*initialCondition)( grid(i,j) );

    solver->setSaveHistory( saveHistory );
    solver->assignExplicitOde( *molProblem, 0.0, initialState );

    writer->initialize( molProblem, solver );
    writeGlobalAttributes();
    writer->writeTimeStep();

    if ( solver->saveHistory() )
    {
        historyFile_.open( std::string( writer->outputPath() + "history.dat").c_str() );
        solver->stats().setHistoryFile( historyFile_ );
    }

}



void
OdeityApplication::declareParameters()
{
    prm.declare_entry( "computation name", "", Patterns::Anything() );
    prm.declare_entry( "final time", "10.0", Patterns::Double() );
    prm.declare_entry( "number of output points", "1", Patterns::Integer() );
    prm.declare_entry( "save results", "true", Patterns::Bool() );
    prm.declare_entry( "save history", "false", Patterns::Bool() );
    prm.declare_entry( "solver", "CVodeGMRES", Patterns::Selection("CVodeGMRES|CVodeBiCG|CVodeTFQMR|RK23|RKM45|DP45|RKC") );

    CVodeSpils::declareParameters( prm );
    RungeKuttaBase::declareParameters( prm );

    AllenCahnEquation::declareParameters( prm );
    CahnHilliardEquation::declareParameters( prm );
    DegenerateCahnHilliardEquation::declareParameters( prm );
    LoretiMarchEquation::declareParameters( prm );

    prm.enter_subsection( "Initial condition" );
        prm.declare_entry( "type", "Wavy circle", Patterns::Selection("Wavy circle|Random|Zero|Sinc|Matheval|Two-scale sine|Rectangle|Two rectangles|Zig-zag") );
        WavyCircleFunction::declareParameters( prm );
        RandomFunction<2>::declareParameters( prm );
        MathevalFunction<2>::declareParameters( prm );
        TwoScaleSineFunction::declareParameters( prm );
        RectangleFunction::declareParameters( prm );
        TwoRectanglesFunction::declareParameters( prm );
    prm.leave_subsection();

    RectangularGrid<2>::declareParameters( prm );
    RectangularDomain<2>::declareParameters( prm );
}



void
OdeityApplication::getParameters()
{
    prm.read_input( configFileName_ );

    computationName_ = prm.get( "computation name" );
    finalTime = prm.get_double("final time");
    numOutputPoints = prm.get_integer("number of output points");
    saveResults = prm.get_bool( "save results" );
    saveHistory = prm.get_bool("save history");
    if ( saveHistory )
        numOutputPoints = 1;

    solver = solverFactory.createObject( prm.get("solver") );
    solver->getParameters( prm );

    prm.enter_subsection( "Initial condition" );
        initialCondition = initCondFactory.createObject( prm.get( "type" ) );
        initialCondition->getParameters( prm );
    prm.leave_subsection();

    domain.getParameters( prm );

    // grid has to be attached to a domain before initializing it from prm file
    grid.attachToRectangularDomain( domain );
    grid.getParameters( prm );
    grid.update();

    molProblem = createMolProblem();
    molProblem->attachGrid( grid );
    molProblem->getParameters( prm );
}



int
OdeityApplication::run()
{
    printHeader();
    ProgressDisplay progDisp( numOutputPoints, std::cerr );
    Timer timer;
    realtype dt = finalTime / numOutputPoints;
    timer.start();
    for ( int i = 1; i <= numOutputPoints; i++)
    {
        solver->integrateTo( i*dt );
        if ( saveResults )
            writer->writeTimeStep();
        ++progDisp;
    }
    timer.stop();

    solver->stats().printInfo();
    logger << std::endl << "*Computational time*:: " << timer() << " seconds" << std::endl;

    writer->writeGlobalAtt( "computational_time", timer() );
    writer->writeGlobalAtt( "jobid", jobid() );

    return 0;
}



void
OdeityApplication::printHeader()
{
    using namespace std;

    std::string title = molProblem->name();
    title += " solved with ODEity";

    logger << title << endl;
    for( size_t i = 0; i < title.length(); i++ ) logger << "=";
    logger << endl;
    logger << "Vladimir Chalupecky <vladimir.chalupecky@fjfi.cvut.cz>" << endl;
    logger << endl << "Computation name: " << computationName_ << endl;
    logger << endl;
    logger << "Time interval" << endl;
    logger << "-------------" << endl;
    logger << "*Final time*:: " << finalTime << endl;
    logger << "*Number of output points*:: " << numOutputPoints << endl;

    domain.printInfo();
    grid.printInfo();
    molProblem->printInfo();
    solver->printInfo();

    logger << endl;
    logger << "Initial condition" << endl;
    logger << "-----------------" << endl;
    initialCondition->printInfo();
}



void OdeityApplication::writeGlobalAttributes()
{
    writer->writeGlobalAtt( "creator", "ODEity" );
    writer->writeGlobalAtt( "program", programName_ );
    writer->writeGlobalAtt( "final_time", finalTime );
    writer->writeGlobalAtt( "absolute_tolerance", solver->absoluteTolerance() );
    writer->writeGlobalAtt( "relative_tolerance", solver->relativeTolerance() );
}




void OdeityApplication::registerInitialConditions()
{
    initCondFactory.registerCreator( "Wavy circle", createWavyCircleFunction );
    initCondFactory.registerCreator( "Zero", createZeroFunction );
    initCondFactory.registerCreator( "Random", createRandomFunction );
    initCondFactory.registerCreator( "Sinc", createSincFunction );
    initCondFactory.registerCreator( "Matheval", createMathevalFunction );
    initCondFactory.registerCreator( "Two-scale sine", createTwoScaleSineFunction );
    initCondFactory.registerCreator( "Rectangle", createRectangleFunction );
    initCondFactory.registerCreator( "Two rectangles", createTwoRectanglesFunction );
    initCondFactory.registerCreator( "Zig-zag", createZigZagFunction );
}



void OdeityApplication::registerOdeSolvers()
{
    solverFactory.registerCreator( "CVodeGMRES", createCVodeGMRESSolver );
    solverFactory.registerCreator( "CVodeBiCG", createCVodeBiCGSolver );
    solverFactory.registerCreator( "CVodeTFQMR", createCVodeTFQMRSolver );
    solverFactory.registerCreator( "RK23", createRungeKutta23Solver );
    solverFactory.registerCreator( "RKM45", createRungeKuttaMerson45Solver );
    solverFactory.registerCreator( "DP45", createDormandPrince45Solver );
    solverFactory.registerCreator( "RKC", createRungeKuttaChebyshevSolver );
}

