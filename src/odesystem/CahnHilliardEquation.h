#ifndef CAHN_HILLIARD_EQUATION_H
#define CAHN_HILLIARD_EQUATION_H

#include "MolOdeSystem.h"
#include "../utils/Vector.h"

class ParameterHandler;

class CahnHilliardEquation : public MolOdeSystem<2>
{
  public:

    CahnHilliardEquation();

    int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );
    int jacobian( const Vector<realtype>& v, Vector<realtype>& Jv,
        const realtype t, const Vector<realtype>& y, const Vector<realtype>& fy,
        const Vector<realtype>& tmp );

    std::string name() const;
    std::string componentName( int i ) const { return "phase_field"; }
    realtype interfaceWidth() const;
    realtype potentialCoefficient() const;

    void printInfo() const;
    static void declareParameters( ParameterHandler& prm );
    void getParameters( ParameterHandler& prm );

    void attachGrid( const RectangularGrid<2>& grid );

  private:

    realtype xi_, xiInv_, xiSqr_;
    realtype a_;
    Vector<realtype> w_;

    realtype f0( const realtype u );
    realtype f0deriv( const realtype u);

};


inline
std::string
CahnHilliardEquation::name() const
{
  return "Cahn-Hilliard equation with constant mobility";
}


inline
realtype
CahnHilliardEquation::interfaceWidth() const
{
  return xi_;
}


inline
realtype
CahnHilliardEquation::potentialCoefficient() const
{
    return a_;
}


#endif // CAHN_HILLIARD_EQUATION_H

