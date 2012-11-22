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

#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <ctime>

#include "NetCDFWriter.h"
#include "../geometry/RectangularGrid.h"
#include "../integrators/OdeIntegratorBase.h"
#include "../odesystem/MolOdeSystem.h"
#include "../utils/JobIdentifier.h"
#include "../utils/Utilities.h"
#include "../utils/Exceptions.h"



template<int N>
NetCDFWriter<N>::NetCDFWriter( const std::string & prefix, const std::string & filename )
:
    odeSystem_( 0 ),
    solver_( 0 ),
    err( NcError::silent_nonfatal ),
    fileName_( filename ),
    prefix_( prefix ),
    recCounter( 0 )
{
    prepareOutputDirectory();
}



template<int N>
void
NetCDFWriter<N>::initialize( const MolOdeSystem<N> * odeSystem, const OdeIntegratorBase * solver )
{
    odeSystem_ = odeSystem;
    solver_ = solver;

    // open NetCDF file
    dataFile = new NcFile( (dirName_+fileName_).c_str(), NcFile::Replace );
    if ( dataFile == 0 || not dataFile->is_valid() )
        exit( 2 );

    std::vector<std::string> dimNames;
    dimNames.push_back("time");
    dimNames.push_back("xdim");
    dimNames.push_back("ydim");
    dimNames.push_back("zdim");

    // write one unlimited time dimension
    if ( !( dims[0] = dataFile->add_dim( dimNames[0].c_str() ) ) )
        exit( 2 );
    // write other spatial dimensions
    for ( int i = 1; i < N+1; ++i )
        if ( !( dims[i] = dataFile->add_dim( dimNames[i].c_str(), odeSystem_->grid()->dimension(i-1) ) ) )
            exit( 2 );

    // auxiliary dimensions for OpenDX fields
    if ( !( ndelta = dataFile->add_dim( "ndelta", 2 ) ) )
        exit( 2 );
    if ( !( ndim = dataFile->add_dim( "ndim", N ) ) )
        exit( 2 );

    // auxiliary positions for OpenDX fields
    if ( !( pos = dataFile->add_var( "pos", ncFloat, ndim, ndelta ) ) )
        exit( 2 );
 
    // define variables
    NcVar * var;
    for ( int i = 0; i < odeSystem_->numberOfComponents(); ++i )
    {
        if ( ! (var = dataFile->add_var( odeSystem_->componentName( i ).c_str(), ncDouble, N+1, (const NcDim**)dims ) ) ) 
            exit( 2 );
        if ( ! var->add_att( "field", (odeSystem_->componentName( i ) + std::string(", scalar, series")).c_str() ) )
            exit( 2 );
        if ( ! var->add_att( "positions", "pos, regular" ) )
            exit( 2 );
        vars.push_back( var );
    }


    float regPos[N][2];
    for ( int i = 0; i < N; ++i )
    {
        regPos[i][0] = 0.0; // TODO: should be spatialStep/2.0 in case we use cell-centered scheme => add a method to RectangularGrid telling us if that is so?
        regPos[i][1] = odeSystem_->grid()->spatialStep(i);
    }
    if ( ! pos->put( &regPos[0][0], N, 2 ) )
        exit( 2 );

    dataFile->sync();

    tmpStorage_.reinit( odeSystem_->grid()->numberOfNodes() );
}



template<int N>
NetCDFWriter<N>::~NetCDFWriter()
{
    if ( dataFile )
    {
        time_t curTime = time( 0 );
        tm * ltime = localtime( &curTime );
        char timeStr[20];
        strftime( timeStr, 20, "%F %T\0", ltime );

        if ( !dataFile->add_att( "finished_datetime", timeStr ) )
            exit( 2 );

        dataFile->sync();
        dataFile->close(); // maybe not necessary

        delete dataFile;
    }
}



// TODO: I don't like this, isn't there a better way how to write this and
// other attributes more generically?
template<int N>
void
NetCDFWriter<N>::writeGlobalAtt( const std::string & name, ncbyte value )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! dataFile->add_att( name.c_str(), value ) )
        exit( 2 );
    dataFile->sync();
}



template<int N>
void
NetCDFWriter<N>::writeGlobalAtt( const std::string & name, char value )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! dataFile->add_att( name.c_str(), value ) )
        exit( 2 );
    dataFile->sync();
}



template<int N>
void
NetCDFWriter<N>::writeGlobalAtt( const std::string & name, short value )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! dataFile->add_att( name.c_str(), value ) )
        exit( 2 );
    dataFile->sync();
}



template<int N>
void
NetCDFWriter<N>::writeGlobalAtt( const std::string & name, int value )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! dataFile->add_att( name.c_str(), value ) )
        exit( 2 );
    dataFile->sync();
}



template<int N>
void
NetCDFWriter<N>::writeGlobalAtt( const std::string & name, float value )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! dataFile->add_att( name.c_str(), value ) )
        exit( 2 );
    dataFile->sync();
}



template<int N>
void
NetCDFWriter<N>::writeGlobalAtt( const std::string & name, double value )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! dataFile->add_att( name.c_str(), value ) )
        exit( 2 );
    dataFile->sync();
}



template<int N>
void
NetCDFWriter<N>::writeGlobalAtt( const std::string & name, const std::string & value )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! dataFile->add_att( name.c_str(), value.c_str() ) )
        exit( 2 );
    dataFile->sync();
}



template<int N>
void
NetCDFWriter<N>::writeTimeStep( )
{
    Assert( dataFile != 0, ExcNotInitialized() );

    if ( ! vars[0]->put_rec( solver_->currentState().data(), recCounter ) )
        exit( 2 );

    // TODO: make it generic for any number of components => need to implement
    // mapping in MolOdeSystem between solution Vectors and Vectors containing
    // just one component
    // for ( int i = 0; i < vars.size(); ++i )
    // {
    //     odeSystem_->
    //   if ( !pressureVar->put_rec( macroTemp.data(), recCounter ) )
    //     exit( 2 );

    // if ( !macConVar->put_rec( macroTemp.data(), recCounter ) )
        // exit( 2 );

    dataFile->sync();
    recCounter++;
}



template<int N>
void
NetCDFWriter<N>::prepareOutputDirectory()
{
    dirName_ = "./results/";
    mkdir( dirName_.c_str(), 0700 );
    if ( prefix_.size() != 0 )
        dirName_ += prefix_ + "-";
    dirName_ += jobid() + "/";
    mkdir( dirName_.c_str(), 0700 );
}


