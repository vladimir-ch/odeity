#include "DataWriter.h"
#include "../geometry/RectangularGrid.h"
#include "../utils/Exceptions.h"

#include <fstream>


template<>
void
DataWriter::writeDxHeader( const std::string& fileNameBase, const RectangularGrid<2>& grid )
{
  using namespace std;

  string headerFilename = fileNameBase + string(".general");

  ofstream headerFile( headerFilename.c_str() );
  AssertThrow( headerFile.is_open(), ExcIO() );

  headerFile << "file = ./" << fileNameBase << ".data" << endl;
  headerFile << "grid = " << grid.dimension( xDim ) << " x " << grid.dimension( yDim ) << endl;
  headerFile << "format = ascii" << endl;
  headerFile << "interleaving = record" << endl; // not required for scalar data
  headerFile << "majority = row" << endl;
  headerFile << "header = lines 2" << endl;
  headerFile << "series = " << times_.size() << ", 1, 1, separator = lines 1" << endl;
  headerFile << "field = " << fieldName_ << endl;
  headerFile << "structure = scalar" << endl;
  headerFile << "type = float" << endl;
  headerFile << "dependency = positions" << endl;
  headerFile << "positions = regular, regular, "
    << grid(0,0)(xDim) << ", " << grid.spatialStep(xDim) << ", "
    << grid(0,0)(yDim) << ", " << grid.spatialStep(yDim) << endl;
  headerFile << endl << "end" << endl;
}



template<>
void
DataWriter::writeDx( const std::string& fileNameBase, const RectangularGrid<2>& grid )
{
  using namespace std;

  writeDxHeader( fileNameBase, grid );

  string dataFilename = fileNameBase + string(".data");
  ofstream dataFile( dataFilename.c_str() );
  AssertThrow( dataFile.is_open(), ExcIO() );

  dataFile << "#Time evolution of " << fieldName_ << endl;

  dataFile.setf(ios::scientific, ios::floatfield);
  dataFile.precision(7);

  for ( size_t step = 0; step < times_.size(); step++ )
  {
    dataFile << "Time step " << step << " at t = " << times_[step] << endl;
    for ( unsigned int x = 0; x < grid.dimension(xDim); x++ )
    {
      for ( unsigned int y = 0; y < grid.dimension(yDim); y++ )
        dataFile << solutions_[step](grid.nodeIndex(x,y)) << " ";
      dataFile << endl;
    }
  }
}



