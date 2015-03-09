#ifndef NETCDFWRITER_H
#define NETCDFWRITER_H

#include "../odesystem/MolOdeSystem.h"
#include "../utils/Vector.h"

#include <netcdfcpp.h>
#include <string>
#include <vector>

class OdeIntegratorBase;

template<int N>
class NetCDFWriter
{
    public:

        NetCDFWriter( const std::string & prefix, const std::string & filename );
        ~NetCDFWriter();

        void initialize( const MolOdeSystem<N> * app, const OdeIntegratorBase * solver ); 
        void writeTimeStep();

        void writeGlobalAtt( const std::string & name, ncbyte value );
        void writeGlobalAtt( const std::string & name, char value );
        void writeGlobalAtt( const std::string & name, short value );
        void writeGlobalAtt( const std::string & name, int value );
        void writeGlobalAtt( const std::string & name, float value );
        void writeGlobalAtt( const std::string & name, double value );
        void writeGlobalAtt( const std::string & name, const std::string & value );

        std::string outputPath() const { return dirName_; }

    private:

        void prepareOutputDirectory();

        const MolOdeSystem<N> * odeSystem_;
        const OdeIntegratorBase * solver_;
        NcError err;
        NcFile *dataFile;
        std::string dirName_;
        std::string fileName_;
        std::string prefix_;
        std::string dateTime_;

        int recCounter;
        NcDim  *timedim;
        NcDim  *ndim;
        NcDim  *ndelta;
        NcVar  *pos;
        NcDim  *dims[N+1];
        std::vector<NcVar*> vars;

        Vector<double> tmpStorage_;
};

#include "NetCDFWriter.impl.h"

#endif // NETCDFWRITER_H

