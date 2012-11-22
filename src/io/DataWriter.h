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

#ifndef DATA_WRITER_H
#define DATA_WRITER_H

#include "../utils/Vector.h"

#include <sundials/sundials_types.h>

#include <string>
#include <vector>


template <int dim> class RectangularGrid;

class DataWriter
{
  public:
    DataWriter( const std::string& fieldName )
      :
        fieldName_( fieldName )
    {}

    void addTimeStep( const realtype time, const Vector<realtype>& solution )
    {
      times_.push_back( time );
      solutions_.push_back( Vector<realtype>( solution ) );
    }

    template<int dim>
    void writeDx( const std::string& fileNameBase, const RectangularGrid<dim>& grid );

  private:

    std::vector<Vector<realtype> > solutions_;
    std::vector<realtype>          times_;
    std::string                    fieldName_;

    template<int dim>
    void writeDxHeader( const std::string& fileNameBase, const RectangularGrid<dim>& grid );
};


#endif // DATA_WRITER_H
