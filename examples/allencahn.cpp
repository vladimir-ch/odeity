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
