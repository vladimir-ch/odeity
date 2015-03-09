#ifndef CVODE_H
#define CVODE_H

#include "OdeIntegratorBase.h"
#include "IntegratorStats.h"

#include <nvector/nvector_serial.h>       /* serial N_Vector types , fct . and macros */

#include <fstream>

extern "C"
int
ExplicitOdeRhs( double t, N_Vector nvY, N_Vector nvYdot, void * f_data );

extern "C"
int
ExplicitOdeJacobian( N_Vector nvV, N_Vector nvJv,
        realtype t, N_Vector nvY, N_Vector nvFy,
        void *jac_data, N_Vector nvTmp );


OdeIntegratorBase * createCVodeGMRESSolver();
OdeIntegratorBase * createCVodeBiCGSolver();
OdeIntegratorBase * createCVodeTFQMRSolver();

class CVodeBase : public OdeIntegratorBase
{
    public:

        // default constructor is protected
        ~CVodeBase();

        void integrateTo( const realtype tOut );
        void setRelativeTolerance( realtype relTol );
        void setAbsoluteTolerance( realtype absTol );
        virtual void updateHistory() = 0;

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

        void printInfo() const;

        DeclException1( ExcCVodeCreateError, char *, << arg1 );
        DeclException1( ExcCVodeInitError, int, << "CVodeInit failed with error code: " << arg1 );
        DeclException1( ExcCVodeReInitError, int, << "CVodeReInit failed with error code: " << arg1 );
        DeclException1( ExcCVodeSetUserDataError, int, << "CVodeSetUserData failed with error code: " << arg1 );
        DeclException1( ExcCVodeSetTolerancesError, int, << "CVodeSetTolerances failed with error code: " << arg1 );
        DeclException1( ExcCVodeError, int, << "CVode failed with error code: " << arg1 );

    protected:

        CVodeBase();

        void assignExplicitOde( 
                ExplicitOde& odeProblem,
                const realtype initialTime,
                const Vector<realtype> &initialState );

        void * cvodeMem_;
        bool initialized_;
        N_Vector nvCurrentState_;

        std::ofstream stepsizeHistoryFile_;
        std::ofstream orderHistoryFile_;

        int maxNumSteps_;

        // no copy constructor
        CVodeBase( const CVodeBase& );
};



class CVodeSpils : public CVodeBase
{
    public:
        virtual ~CVodeSpils();

        IntegratorStatsBase& stats();
        void updateHistory();

        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

        void printInfo() const;

        DeclException1( ExcCVSpgmrError, int, << "CVSpgmr failed with error code: " << arg1 );
        DeclException1( ExcCVSpbcgError, int, << "CVSpbcg failed with error code: " << arg1 );
        DeclException1( ExcCVSptfqmrError, int, << "CVStfqmr failed with error code: " << arg1 );

    protected:

        CVodeSpils();
        void assignExplicitOde( 
                ExplicitOde& odeProblem,
                const realtype initialTime,
                const Vector<realtype> &initialState );

        bool precond_;
        int krylovSubspaceDim_;

    private:

        CVodeSpilsStats stats_;

        CVodeSpils( const CVodeSpils& );
        CVodeSpils& operator = ( const CVodeSpils& );
};



class CVodeGMRES : public CVodeSpils
{
    public:
        CVodeGMRES();

        void assignExplicitOde( 
                ExplicitOde& odeProblem,
                const realtype initialTime,
                const Vector<realtype> &initialState );
};



class CVodeBiCG : public CVodeSpils
{
    public:
        CVodeBiCG();

        void assignExplicitOde( 
                ExplicitOde& odeProblem,
                const realtype initialTime,
                const Vector<realtype> &initialState );
};



class CVodeTFQMR : public CVodeSpils
{
    public:
        CVodeTFQMR();

        void assignExplicitOde( 
                ExplicitOde& odeProblem,
                const realtype initialTime,
                const Vector<realtype> &initialState );    
};


#endif // CVODE_H

