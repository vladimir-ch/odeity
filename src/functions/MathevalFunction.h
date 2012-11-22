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

