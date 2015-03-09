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

