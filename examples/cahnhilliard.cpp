#include <io/NetCDFWriter.h>
#include <odesystem/CahnHilliardEquation.h>
#include <utils/OdeityApplication.h>


class CahnHilliardApplication : public OdeityApplication
{
    public:

        CahnHilliardApplication( int argc, char * argv[] )
            :
                OdeityApplication( argc, argv )
        { }


    private:

        MolOdeSystem<2>* createMolProblem() const
        {
            return new CahnHilliardEquation();
        }

        void writeGlobalAttributes()
        {
            OdeityApplication::writeGlobalAttributes();

            CahnHilliardEquation * eq = (CahnHilliardEquation*) molProblem;
            writer->writeGlobalAtt( "xi", eq->interfaceWidth() );
            writer->writeGlobalAtt( "a", eq->potentialCoefficient() );
        }
};


int main( int argc, char * argv[] )
{
  CahnHilliardApplication * app = new CahnHilliardApplication( argc, argv );

  app->initialize();
  return app->run();
}
