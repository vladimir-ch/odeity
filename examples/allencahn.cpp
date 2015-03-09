#include <io/NetCDFWriter.h>
#include <odesystem/AllenCahnEquation.h>
#include <utils/OdeityApplication.h>


using namespace std;

class AllenCahnApplication : public OdeityApplication
{
    public:

        AllenCahnApplication( int argc, char * argv[] )
            :
                OdeityApplication( argc, argv )
        {}

    private:

        MolOdeSystem<2>* createMolProblem() const
        {
            return new AllenCahnEquation();
        }


        void writeGlobalAttributes()
        {
            OdeityApplication::writeGlobalAttributes();

            AllenCahnEquation * eq = (AllenCahnEquation*) molProblem;
            writer->writeGlobalAtt( "xi", eq->interfaceWidth() );
            writer->writeGlobalAtt( "F", eq->forcingTerm() );
        }
};


int main( int argc, char * argv[])
{
    try
    {
        AllenCahnApplication * app = new AllenCahnApplication( argc, argv );
        app->initialize();
        app->run();
    }
    catch( std::exception &exc )
    {
        std::cerr << std::endl << std::endl
                  << "--------------------------------------------------"
                  << std::endl
                  << "Exception encountered: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting" << std::endl
                  << "--------------------------------------------------"
                  << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << std::endl << std::endl
                  << "--------------------------------------------------"
                  << std::endl
                  << "Unknown exception encountered" << std::endl
                  << "Aborting" << std::endl
                  << "--------------------------------------------------"
                  << std::endl;

        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
