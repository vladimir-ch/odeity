#ifndef MATHEVAL_FUNCTION_H
#define MATHEVAL_FUNCTION_H

#include "Function.h"

#include <matheval.h>
#include <string>


/** \brief Encapsulates GNU matheval library
 *
 * Class interfacing GNU matheval library for parsing and evaluating functions
 * given in a string.
 */
template<int N>
class MathevalFunction : public Function<N>
{
    public:

        // ====================  LIFECYCLE     =======================================
        MathevalFunction () : f_( 0 ), varNames_( 0 ), varCount_( 0 )
        {}

        ~MathevalFunction () { evaluator_destroy( f_ ); }

        // ====================  OPERATIONS    =======================================
        static void declareParameters( ParameterHandler & prm );
        void getParameters( ParameterHandler & prm );

        void printInfo() const;

        // ====================  INQUIRY       =======================================
        realtype value ( const Point<N> & p, const unsigned int component = 0 );

    private:
        void * f_;
        std::string fString_;
        char **varNames_;
        int varCount_;

};

#include "MathevalFunction.impl.h"

template<int N>
Function<N> *
createMathevalFunction()
{
    return new MathevalFunction<N>();
}

#endif // MATHEVAL_FUNCTION_H

