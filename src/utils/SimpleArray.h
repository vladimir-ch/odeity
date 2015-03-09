#ifndef SIMPLE_ARRAY_H
#define SIMPLE_ARRAY_H

#include "Exceptions.h"

#include <ostream>

template < typename T >
class SimpleArray {
  public:
    SimpleArray( int rows, int cols )
      :
        rows_(0),
        cols_(0),
        data_(0)
    {
      resize( rows, cols );
    }

    ~SimpleArray()
    {
      delete [] data_;
    }

    void resize( int rows, int cols )
    {
      if ( rows != rows_ || cols != cols_ )
      {
        if ( data_ )
          delete [] data_;
        data_ = new T[rows*cols];
        if ( data_ == 0 )
          throw ExcOutOfMemory();
        rows_ = rows;
        cols_ = cols;
        for ( int i = 0 ; i < rows*cols; i++ )
          data_[i] = 0;
      }
    }

    void operator=(T value)
    {
      Assert (rows_!=0 && cols_!=0, ExcEmptyObject());

      for ( int i = 0; i < rows_*cols_; i++ )
        data_[i] = value;
    }


    T& operator()( int i, int j )
    {
      Assert (rows_!=0 && cols_!=0, ExcEmptyObject());
      Assert(0<=i && i<rows_,ExcIndexRange(i,0,rows_));
      Assert(0<=j && j<rows_,ExcIndexRange(j,0,cols_));

      return data_[i*cols_+j];
    }

    const T& operator()( int i, int j ) const
    {
      Assert (rows_!=0 && cols_!=0, ExcEmptyObject());
      Assert(0<=i && i<rows_,ExcIndexRange(i,0,rows_));
      Assert(0<=j && j<rows_,ExcIndexRange(j,0,cols_));

      return data_[i*cols_+j];
    }

    int rows() const
    {
      return rows_;
    }

    int cols() const
    {
      return cols_;
    }

  protected:

  private:
    int rows_;
    int cols_;
    T *data_;
};

template < typename T > 
std::ostream& operator<<( std::ostream& out, const SimpleArray<T>& a)
{
  out << "SimpleArray<T>: rows = " << a.rows() << ", cols = " << a.cols() << std::endl;
  for ( int i = 0; i < a.rows(); i++ )
  {
    for ( int j = 0; j < a.cols(); j++ )
      out << a(i,j) << " ";

    out << std::endl;
  }
  return out;
}

#endif // SIMPLE_ARRAY_H

