#ifndef DEG_CAHN_HILLIARD_EQUATION_H
#define DEG_CAHN_HILLIARD_EQUATION_H

#include "MolOdeSystem.h"
#include "../utils/Vector.h"

class ParameterHandler;

class DegenerateCahnHilliardEquation : public MolOdeSystem<2>
{
    public:

        DegenerateCahnHilliardEquation();

        int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );

        std::string name() const;
        std::string componentName( int i ) const { return "phase_field"; }

        realtype interfaceWidth() const;
        realtype alpha() const;
        realtype mobilityCoefficient() const;

        void printInfo() const;
        static void declareParameters( ParameterHandler& prm );
        void getParameters( ParameterHandler& prm );
        void attachGrid( const RectangularGrid<2>& grid );

    private:

        realtype xi_, xiPow2Inv_, xiPow2_;
        realtype xiPowAlpha_;
        realtype alpha_, beta_;
        realtype eps_;
        Vector<realtype> w_;
        Vector<realtype> mAvgX_;
        Vector<realtype> mAvgY_;

        realtype f0( const realtype u );
        realtype M( const realtype u );
};



inline
std::string
DegenerateCahnHilliardEquation::name() const
{
    return "Cahn-Hilliard equation with degenerate mobility";
}



inline
realtype
DegenerateCahnHilliardEquation::interfaceWidth() const
{
    return xi_;
}



inline
realtype
DegenerateCahnHilliardEquation::alpha() const
{
    return alpha_;
}



inline
realtype
DegenerateCahnHilliardEquation::mobilityCoefficient() const
{
    return beta_;
}



#endif // DEG_CAHN_HILLIARD_EQUATION_H


