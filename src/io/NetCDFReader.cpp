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

#include <iomanip>
#include <iostream>

#include "NetCDFReader.h"

using namespace std;

NetCDFReader::NetCDFReader( const string & fileName, const string & variableName ) :
    fileName_( fileName ),
    variableName_( variableName ),
    err( NcError::verbose_nonfatal ),
    dataFile( fileName.c_str(), NcFile::ReadOnly )
{
  // verify if file is valid
  if( not dataFile.is_valid() )
  {
    cerr << "File " << fileName_ << " cannot be opened" << endl;
    exit( EXIT_FAILURE );
  }

  if ( not (creatorAtt = dataFile.get_att("creator")) )
  {
    cerr << "Attribute \"creator\" cannot be read" << endl;
    exit( EXIT_FAILURE );
  }
  
  string creatorAttStr( creatorAtt->as_string(0) );
  if ( creatorAttStr != "ODEity" )
  {
    cerr << "NetCDF file " << fileName_ << " was not written by ODEity" << endl;
    exit( EXIT_FAILURE );
  }

  // read variable
  if ( not (var = dataFile.get_var(variableName_.c_str()) ) )
  {
    cerr << "Variable " << variableName_ << " cannot be read" << endl;
    exit( EXIT_FAILURE );
  }
  if ( not (posVar = dataFile.get_var("pos")) )
  {
    cerr << "Variable pos cannot be read" << endl;
    exit( EXIT_FAILURE );
  }

  // read dimensions
  if ( not (xdim = dataFile.get_dim("xdim")) )
  {
    cerr << "Dimension xdim cannot be read" << endl;
    exit( EXIT_FAILURE );
  }
  if ( not (ydim = dataFile.get_dim("ydim")) )
  {
    cerr << "Dimension ydim cannot be read" << endl;
    exit( EXIT_FAILURE );
  }
  if ( not (timedim = dataFile.rec_dim() ) )
  {
    cerr << "Unlimited dimension time cannot be read" << endl;
    exit( EXIT_FAILURE );
  }

  // these two variables are not needed right now, but maybe in the future, when
  // the implementation of the reader is more generic and dimension-independent
  if ( not (ndim = dataFile.get_dim("ndim")) )
  {
    cerr << "Dimension ndim cannot be read" << endl;
    exit( EXIT_FAILURE );
  }
  if ( not (ndelta = dataFile.get_dim("ndelta")) )
  {
    cerr << "Dimension ndelta cannot be read" << endl;
    exit( EXIT_FAILURE );
  }

  // get pos variable
  if ( not posVar->get( &pos[0][0], 2, 2) )
  {
    cerr << "Cannot get pos variable" << endl;
    exit( EXIT_FAILURE );
  }
}



void
NetCDFReader::getGridDimensions( int & dimX, int & dimY ) const
{
  dimX = xdim->size();
  dimY = ydim->size();
}



void
NetCDFReader::getGridSpacing( float & hX, float & hY ) const
{
  hX = pos[0][1];
  hY = pos[1][1];
}



std::string
NetCDFReader::fileName() const
{
  return fileName_;
}



void
NetCDFReader::readTimeStep( double * data, int timeStep ) const
{
  if (!var->set_cur(timeStep, 0, 0))
  {
    cerr << "Cannot set record " << timeStep << " from variable" << variableName_ << endl;
    exit( EXIT_FAILURE );
  }

  if (!var->get( data, 1, xdim->size(), ydim->size() ))
  {
    cerr << "Cannot read record " << timeStep << " from variable" << variableName_ << endl;
    exit( EXIT_FAILURE );
  }
}



int
NetCDFReader::numberOfTimeSteps() const
{
  return timedim->size();
}



double
NetCDFReader::finalTime() const
{
  return dataFile.get_att("final_time")->as_double(0);
}


