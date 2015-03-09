#ifndef INTEGRATOR_STATS_H
#define INTEGRATOR_STATS_H

#include <sundials/sundials_types.h>
#include <cvode/cvode.h>                  /* main integrator header file */
#include <cvode/cvode_spgmr.h>            /* prototypes & constants for CVSPGMR solver */
#include <cvode/cvode_spbcgs.h>           /* prototypes & constants for CVSPBCG solver */
#include <cvode/cvode_sptfqmr.h>          /* prototypes & constants for CVSTFQMR solver */
#include <cvode/cvode_band.h>
#include <cvode/cvode_bandpre.h>

#include <iosfwd>

class IntegratorStatsBase
{
    public:
        IntegratorStatsBase();

        virtual void printInfo() const = 0;
        void setHistoryFile( std::ofstream & historyFile );

    protected:

        void updateTimeStepsize( realtype time, realtype stepsize );
        virtual void reset();

        std::ofstream * historyFile_;
        realtype maxTimeStepsize_;
        realtype minTimeStepsize_;
};



class RungeKuttaStats : public IntegratorStatsBase
{
    public:

        void printInfo() const;

    protected:

        RungeKuttaStats();

        void updateHistory( realtype time, realtype stepsize );
        void reset();

        void setRhsEvaluations( int numRhsEval );
        void incRhsEvaluations();
        int getRhsEvaluations() const;

        void setAcceptedSteps( int numAccSteps );
        void incAcceptedSteps();
        int getAcceptedSteps() const;

        void setRejectedSteps( int numRejSteps );
        void incRejectedSteps();
        int getRejectedSteps() const;

        int rhsEvaluations_;
        int acceptedSteps_;
        int rejectedSteps_;

        friend class ExplicitRungeKuttaBase;
};



class RungeKuttaChebyshevStats : public RungeKuttaStats
{
    public:

        void printInfo() const;

    protected:

        int nspraditers_;
        int maxStage_;

        RungeKuttaChebyshevStats();

        void updateHistory( realtype time, realtype stepsize, int stage );
        void reset();

        friend class RungeKuttaChebyshev;
};



class CVodeStats : public IntegratorStatsBase
{
    protected:

        CVodeStats();

        void updateHistory( realtype time, realtype stepsize, int order );
        void reset();
        virtual void update( void * cvodeMem );
        void printInfo() const;

        // CVode stats
        long int nsteps;
        long int nfevals;
        long int nlinsetups;
        long int netfails;

        // NonlinSolvStats
        long int nniters;
        long int nncfails;

        int maxOrder_;
};



class CVodeBandStats : public CVodeStats
{
    public:

        void printInfo() const;

    protected:

        CVodeBandStats() { reset(); }

        void reset()
        {
            njevals = 0;
            nfevalsLS = 0;
        }

        void update( void * cvodeMem );

        long int njevals;
        long int nfevalsLS;

        friend class CVodeBand;
};



class CVodeSpilsStats : public CVodeStats
{
    public:

        void printInfo() const;

    protected:

        CVodeSpilsStats() { reset(); }

        void reset()
        {
            nliters = 0;
            nlcfails = 0;
            npevals = 0;
            npsolves = 0;
            njvevals = 0;
            nfevalsLS = 0;
            nfevalsBP = 0;
        }

        void update( void * cvodeMem );
        void update( void * cvodeMem, void * bpData );

        long int nliters;
        long int nlcfails;
        long int npevals;
        long int npsolves;
        long int njvevals;
        long int nfevalsLS;

        // BANDPRE stats
        long int nfevalsBP;

        friend class CVodeSpils;
};



/* ------------------ Inline functions -------------------- */

inline
void RungeKuttaStats::setRhsEvaluations( int numRhsEval )
{
    rhsEvaluations_ = numRhsEval;
}


inline
void RungeKuttaStats::incRhsEvaluations()
{
    ++rhsEvaluations_;
}

inline
int RungeKuttaStats::getRhsEvaluations() const
{
    return rhsEvaluations_;
}

inline
void RungeKuttaStats::setAcceptedSteps( int numAccSteps )
{
    acceptedSteps_ = numAccSteps;
}

inline
void RungeKuttaStats::incAcceptedSteps()
{
    ++acceptedSteps_;
}

inline
int RungeKuttaStats::getAcceptedSteps() const
{
    return acceptedSteps_;
}

inline
void RungeKuttaStats::setRejectedSteps( int numRejSteps )
{
    rejectedSteps_ = numRejSteps;
}

inline
void RungeKuttaStats::incRejectedSteps()
{
    ++rejectedSteps_;
}

inline
int RungeKuttaStats::getRejectedSteps() const
{
    return rejectedSteps_;
}


#endif // INTEGRATOR_STATS_H

