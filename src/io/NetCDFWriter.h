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

