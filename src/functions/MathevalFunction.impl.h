#include <iostream>

/** \brief Declare parameters in the ParameterHandler
 */
template<int N>
void MathevalFunction<N>::declareParameters ( ParameterHandler & prm )
{
    prm.enter_subsection( "Matheval" );
        prm.declare_entry( "f(x)", "sin(x)*sin(y)", Patterns::Anything() );
    prm.leave_subsection();

    return ;
}


/** \brief Get parameters from the ParameterHandler
 */
template<int N>
void MathevalFunction<N>::getParameters ( ParameterHandler & prm )
{
    prm.enter_subsection( "Matheval" );
        fString_ = prm.get( "f(x)" );
    prm.leave_subsection();

    f_ = evaluator_create( const_cast<char*>( fString_.c_str() ) );
    if ( f_ == 0 )
    {
        std::cerr << "Error in function expression: " << fString_ << std::endl;
        exit( 2 );
    }
    
    evaluator_get_variables( f_, &varNames_, &varCount_ );
    if ( varCount_ != N )
    {
        std::cerr << "Number of variables in function expression: " << fString_ << " does not agree with dimension of Point<" << N << ">" << std::endl;
        exit( 2 );
    }

    return;
}


/** \brief Evaluate the function at the given Point
 */
template<int N>
realtype MathevalFunction<N>::value ( const Point<N> & p, const unsigned int component )
{
    return evaluator_evaluate( f_, varCount_, varNames_, const_cast<double*>( p.values() ) );
}


/** \brief Print object information to LogStream
 */
template<int N>
void MathevalFunction<N>::printInfo() const
{
    using namespace std;

    logger << endl;
    logger << "Matheval function" << endl;
    logger << "~~~~~~~~~~~~~~~~~" << endl;
    logger << "*Function expression*:: f(x) = " << fString_ << endl;
}
