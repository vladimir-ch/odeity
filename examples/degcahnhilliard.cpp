#include <io/NetCDFWriter.h>
#include <utils/OdeityApplication.h>
#include <odesystem/DegenerateCahnHilliardEquation.h>


class DegenerateCahnHilliardApplication : public OdeityApplication
{
    public:

        DegenerateCahnHilliardApplication( int argc, char * argv[] )
            :
                OdeityApplication( argc, argv )
        {}

    private:

        MolOdeSystem<2>* createMolProblem() const
        {
            return new DegenerateCahnHilliardEquation();
        }


        void writeGlobalAttributes()
        {
            OdeityApplication::writeGlobalAttributes();

            DegenerateCahnHilliardEquation * eq = (DegenerateCahnHilliardEquation*) molProblem;
            writer->writeGlobalAtt( "xi", eq->interfaceWidth() );
            writer->writeGlobalAtt( "alpha", eq->alpha() );
            writer->writeGlobalAtt( "beta", eq->mobilityCoefficient() );
        }
};


int main( int argc, char * argv[] )
{
  DegenerateCahnHilliardApplication * app = new DegenerateCahnHilliardApplication( argc, argv );

  app->initialize();
  return app->run();
}

