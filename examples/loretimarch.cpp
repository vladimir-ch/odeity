#include <io/NetCDFWriter.h>
#include <odesystem/LoretiMarchEquation.h>
#include <utils/OdeityApplication.h>


class LoretiMarchApplication : public OdeityApplication
{
    public:

        LoretiMarchApplication( int argc, char * argv[] )
            :
                OdeityApplication( argc, argv )
        { }


    private:

        MolOdeSystem<2>* createMolProblem() const
        {
            return new LoretiMarchEquation();
        }

        void writeGlobalAttributes()
        {
            OdeityApplication::writeGlobalAttributes();

            LoretiMarchEquation * eq = (LoretiMarchEquation*) molProblem;
            writer->writeGlobalAtt( "xi", eq->interfaceWidth() );
            writer->writeGlobalAtt( "a", eq->potentialCoefficient() );
        }
};


int main( int argc, char * argv[] )
{
  LoretiMarchApplication * app = new LoretiMarchApplication( argc, argv );

  app->initialize();
  return app->run();
}
