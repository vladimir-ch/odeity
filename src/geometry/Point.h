#ifndef POINT_H
#define POINT_H

#include "../utils/Exceptions.h"

#include <sundials/sundials_types.h>

#include <ostream>
#include <cmath>


template<int dim>
class Point
{
  public:

    static const unsigned int dimension = dim;

    typedef realtype array_type[dim];

    Point();
    explicit Point( const bool initialize );
    Point( const array_type& initializer );

    Point( const Point<dim>& );

    explicit Point( const realtype x );
    Point( const realtype x, realtype y );
    Point( const realtype x, realtype y, realtype z );

    realtype operator () ( const unsigned int index ) const;
    realtype& operator () ( const unsigned int index );
    realtype operator [] ( const unsigned int index ) const;
    realtype& operator [] ( const unsigned int index );

    Point<dim> & operator = ( const Point<dim> & );
    bool operator == ( const Point<dim> & ) const;
    bool operator != ( const Point<dim> & ) const;
    Point<dim> & operator += ( const Point<dim> & );
    Point<dim> & operator -= ( const Point<dim> & );
    Point<dim> & operator *= ( const double factor );
    Point<dim> & operator /= ( const double factor );

    Point<dim> operator + ( const Point<dim>& ) const;
    Point<dim> operator - ( const Point<dim>& ) const;
    Point<dim> operator - () const;
    realtype operator * ( const Point<dim>& ) const;
    Point<dim> operator * ( const realtype ) const;
    Point<dim> operator / ( const realtype ) const;

    void zeros();
    realtype norm_square() const;
    realtype norm() const;
    realtype square() const;
    realtype distance( const Point<dim>& p ) const;

    const realtype * values() const { return values_; }
    
    DeclException0( ExcInvalidConstructorCalled );

  private:

    realtype values_[dim];
};

/* ---------------------- Inline functions ----------------------------------*/

template< int dim >
inline
Point<dim>::Point()
{}

template< int dim >
inline
Point<dim>::Point( const bool initialize )
{
  if ( initialize )
    for ( unsigned int i = 0; i < dim; i++ )
      values_[i] = 0.0;
}

template< int dim >
inline
Point<dim>::Point( const Point<dim>& p )
{
  for ( int i = 0; i < dim; i++ )
    values_[i] = p.values_[i];
}

template< int dim >
inline
Point<dim>::Point( const realtype x )
{
  Assert( dim == 1, ExcInvalidConstructorCalled() );
  values_[0] = x;
}

template< int dim >
inline
Point<dim>::Point( const realtype x, const realtype y )
{
  Assert( dim == 2, ExcInvalidConstructorCalled() );
  values_[0] = x;
  values_[1] = y;
}

template< int dim >
inline
Point<dim>::Point( const realtype x, const realtype y, const realtype z )
{
  Assert( dim == 3, ExcInvalidConstructorCalled() );
  values_[0] = x;
  values_[1] = y;
  values_[2] = z;
}

template< int dim >
inline
realtype
Point<dim>::operator()( const unsigned int index ) const
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );

  return values_[ index ];
}

template< int dim >
inline
realtype&
Point<dim>::operator()( const unsigned int index )
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );

  return values_[ index ];
}

template< int dim >
inline
realtype
Point<dim>::operator[]( const unsigned int index ) const
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );

  return values_[ index ];
}

template< int dim >
inline
realtype&
Point<dim>::operator[]( const unsigned int index )
{
  Assert( index < dim, ExcIndexRange( index, 0, dim ) );

  return values_[ index ];
}



template<>
inline
Point<1> &
Point<1>::operator = ( const Point<1> & p )
{
  values_[0] = p.values_[0];
  return *this;
}

template<>
inline
Point<2> &
Point<2>::operator = ( const Point<2> & p )
{
  values_[0] = p.values_[0];
  values_[1] = p.values_[1];
  return *this;
}

template<>
inline
Point<3> &
Point<3>::operator = ( const Point<3> & p )
{
  values_[0] = p.values_[0];
  values_[1] = p.values_[1];
  values_[2] = p.values_[2];
  return *this;
}

template< int dim >
inline
Point<dim> &
Point<dim>::operator = ( const Point<dim> & p )
{
  for( unsigned int i = 0; i < dim; ++i )
    values_[i] = p.values_[i];

  return *this;
}


template< int dim >
inline
bool
Point<dim>::operator == ( const Point<dim> & p ) const
{
  for ( unsigned int i = 0; i < dim; i++ )
    if ( values_[i] != p.values_[i] )
      return false;
  return true;
}


template< int dim >
inline
bool
Point<dim>::operator != ( const Point<dim> & p ) const
{
  return !((*this) == p);
}


template< int dim >
inline
Point<dim> &
Point<dim>::operator += ( const Point<dim> & p )
{
  for ( unsigned int i = 0; i < dim; i++ )
    values_[i] += p.values_[i];
  return *this;
}


template< int dim >
inline
Point<dim> &
Point<dim>::operator -= ( const Point<dim> & p )
{
  for ( unsigned int i = 0; i < dim; i++ )
    values_[i] -= p.values_[i];
  return *this;
}


template< int dim >
inline
Point<dim> &
Point<dim>::operator *= ( const double factor )
{
  for ( unsigned int i = 0; i < dim; i++ )
    values_[i] *= factor;
  return *this;
}


template< int dim >
inline
Point<dim> &
Point<dim>::operator /= ( const double factor )
{
  for ( unsigned int i = 0; i < dim; i++ )
    values_[i] /= factor;
  return *this;
}



template< int dim >
inline
Point<dim>
Point<dim>::operator + ( const Point<dim> & p ) const
{
  return (Point<dim>(*this) += p);
}

template< int dim >
inline
Point<dim>
Point<dim>:: operator - ( const Point<dim> & p ) const
{
  return (Point<dim>(*this) -= p);
}

template< int dim >
inline
Point<dim>
Point<dim>:: operator - () const
{
  Point<dim> result;

  for( unsigned int i = 0; i < dim; ++i )
    result.values_[i] = -this->values_[i];

  return result;
}


template<>
inline
realtype
Point<1>:: operator * ( const Point<1> & p ) const
{
  return values_[0]*p.values_[0];
}

template<>
inline
realtype
Point<2>:: operator * ( const Point<2> & p ) const
{
  return values_[0]*p.values_[0] +
         values_[1]*p.values_[1];
}

template<>
inline
realtype
Point<3>:: operator * ( const Point<3> & p ) const
{
  return values_[0]*p.values_[0] +
         values_[1]*p.values_[1] +
         values_[2]*p.values_[2];
}

template< int dim >
inline
realtype
Point<dim>:: operator * ( const Point<dim> & p ) const
{
  realtype sum = 0.0;
  for ( unsigned int i = 0; i < dim; i++ )
    sum += values_[i] * p.values_[i];
  return sum;
}

template< int dim >
inline
Point<dim>
Point<dim>:: operator * ( const realtype factor ) const
{
  return (Point<dim>(*this) *= factor);
}

template< int dim >
inline
Point<dim>
Point<dim>:: operator / ( const realtype factor ) const
{
  return (Point<dim>(*this) /= factor);
}

template< int dim >
inline
void
Point<dim>::zeros()
{
  for( unsigned int i = 0; i < dim; i++ )
    values_[i] = 0.0;
}

template< int dim >
inline
realtype
Point<dim>::norm_square() const
{
  double q = 0;
  for( unsigned int i = 0; i < dim; ++i )
    q += values_[i] * values_[i];
  return q;
}

template< int dim >
inline
realtype
Point<dim>::norm() const
{
  return std::sqrt( norm_square() );
}

template< int dim >
inline
realtype
Point<dim>::square() const
{
  return norm_square();
}

template< int dim >
inline
realtype
Point<dim>::distance( const Point<dim> & p ) const
{
  double sum = 0;

  for( unsigned int i = 0; i < dim; ++i )
  {
    const double diff = values_[i] - p.values_[i];
    sum += diff*diff;
  }
  
  return std::sqrt(sum);
}

/* ---------------------- Global functions: Point ---------------------------*/

template <int dim>
inline
Point<dim> operator * (const double factor, const Point<dim> & p)
{
  return p*factor;
}

template <int dim>
inline
std::ostream & operator << ( std::ostream &out, const Point<dim> & p)
{
  for( unsigned int i = 0; i < dim-1; ++i )
    out << p(i) << ' ';
  out << p(dim-1);

  return out;
}

inline
std::ostream & operator << (std::ostream &out, const Point<1> & p)
{
  out << p(0);

  return out;
}

#endif // POINT_H
