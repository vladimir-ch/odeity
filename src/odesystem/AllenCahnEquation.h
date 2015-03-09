#ifndef ALLEN_CAHN_EQUATION_H
#define ALLEN_CAHN_EQUATION_H

#include "MolOdeSystem.h"

class ParameterHandler;

class AllenCahnEquation : public MolOdeSystem<2>
{
    public:

        AllenCahnEquation();

        virtual int rhs( const realtype &t, const Vector<realtype> &Y, Vector<realtype> &Ydot );

        realtype interfaceWidth() const;
        realtype forcingTerm() const;

        std::string name() const;
        void printInfo() const;
        static void declareParameters( ParameterHandler& prm );
        void getParameters( ParameterHandler& prm );

        std::string componentName( int i ) const { return std::string("phase_field"); }

    private:

        realtype xi_, xiSqrInv_;
        realtype F_;
};


inline
std::string
AllenCahnEquation::name() const
{
    return "Allen-Cahn equation";
}


inline
realtype
AllenCahnEquation::interfaceWidth() const
{
    return xi_;
}



inline
realtype
AllenCahnEquation::forcingTerm() const
{
    return F_;
}



#endif // ALLEN_CAHN_EQUATION_H

