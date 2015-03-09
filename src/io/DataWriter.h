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
