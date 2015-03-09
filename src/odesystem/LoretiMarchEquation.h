#ifndef LORETI_MARCH_EQUATION_H
#define LORETI_MARCH_EQUATION_H

#include "MolOdeSystem.h"
#include "../utils/Vector.h"

class ParameterHandler;

class LoretiMarchEquation : public MolOdeSystem<2>
{
    public:

        LoretiMarchEquation();

        int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );

        std::string name() const;
        std::string componentName( int i ) const { return "phase_field"; }
        realtype interfaceWidth() const;
        realtype potentialCoefficient() const;

        void printInfo() const;
        static void declareParameters( ParameterHandler& prm );
        void getParameters( ParameterHandler& prm );

        void attachGrid( const RectangularGrid<2>& grid );

    private:

        realtype xi_, xiInv_, xiSqr_, twoOverXi_, twoOverXiPow3_;
        realtype a_;
        Vector<realtype> w_;

        realtype PsiDer( const realtype u);
        realtype PsiDerDer( const realtype u);

};


inline
std::string
LoretiMarchEquation::name() const
{
  return "Loreti-March equation for motion V=H-2H_{ss}-H^3";
}


inline
realtype
LoretiMarchEquation::interfaceWidth() const
{
  return xi_;
}


inline
realtype
LoretiMarchEquation::potentialCoefficient() const
{
    return a_;
}


#endif // LORETI_MARCH_EQUATION_H

