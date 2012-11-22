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

