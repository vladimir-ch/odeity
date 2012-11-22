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

#ifndef NETCDF_READER_H
#define NETCDF_READER_H

#include <nvector/nvector_serial.h>

#include <string>
#include <netcdfcpp.h>
#include <iostream>

// TODO: for now the reader is implemented to read 2D fields only. Make it
// generic to read 1D, 2D and 3D fields

class NetCDFReader
{
  public:
    NetCDFReader( const std::string & fileName, const std::string & variableName );

    void getGridDimensions( int & dimX, int & dimY ) const;
    void getGridSpacing( float & hX, float & hY ) const;
    int numberOfTimeSteps() const;
    double finalTime() const;
    void readTimeStep( double * data, int timeStep ) const;

    double getAttributeDouble( const std::string & attName )
    {
        NcAtt * att;
        
        if ( not ( att = dataFile.get_att( attName.c_str() ) ) )
        {
            std::cerr << "Attribute " << attName << " cannot be read" << std::endl;
            exit( EXIT_FAILURE );
        }

        return att->as_double(0);
    }

    std::string fileName() const;
    std::string variableName() const;

  private:
    std::string fileName_;
    std::string variableName_;

    NcError err;
    NcFile  dataFile;
 
    NcDim  *xdim;
    NcDim  *ydim;
    NcDim  *timedim;

    NcDim  *ndim;
    NcDim  *ndelta;

    NcVar  *posVar;
    NcVar  *var;

    NcAtt  *creatorAtt;

    float  pos[2][2];
};


#endif

