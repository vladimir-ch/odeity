//---------------------------------------------------------------------------
//    vector.templates.h,v 1.52 2005/03/29 00:04:22 guido Exp
//    Version: Version-5-2-0
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef VECTOR_TEMPLATES_H
#define VECTOR_TEMPLATES_H

#include "Vector.h"

#include <sundials/sundials_types.h>

#include <cmath>
#include <algorithm>
#include <iostream>

/*
   Note that in this file, we use std::fabs, std::sqrt, etc
   everywhere. The reason is that we want to use those version of these
   functions that take a variable of the template type "Number", rather
   than the C standard function which accepts and returns a double. The
   C++ standard library therefore offers overloaded versions of these
   functions taking floats, or long doubles, with the importance on the
   additional accuracy when using long doubles.
   */

namespace internal
{
  namespace VectorHelper
  {
    template <typename Number>
      inline Number sqr (const Number x)
      {
        return x*x;
      }
  }
}


  template <typename Number>
Vector<Number>::Vector (const Vector<Number>& v)
  :
    dim(v.size()),
    maxdim(v.size()),
    val(0),
    dataOwner(true)
{
  if (dim != 0)
  {
    val = new Number[maxdim];
    Assert (val != 0, ExcOutOfMemory());
    std::copy (v.begin(), v.end(), begin());
  }
}


template< typename Number>
Vector<Number>::Vector( N_Vector v )
  :
    dim( NV_LENGTH_S(v) ),
    maxdim( NV_LENGTH_S(v) ),
    val( reinterpret_cast<Number*>(NV_DATA_S(v)) ),
    dataOwner(false)
{}


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG

template <typename Number>
  template <typename OtherNumber>
Vector<Number>::Vector (const Vector<OtherNumber>& v)
  :
    dim(v.size()),
    maxdim(v.size()),
    val(0),
    dataOwner(true)
{
  if (dim != 0)
  {
    val = new Number[maxdim];
    Assert (val != 0, ExcOutOfMemory());
    std::copy (v.begin(), v.end(), begin());
  }
}

#endif



template <typename Number>
  void
Vector<Number>::swap (Vector<Number> &v)
{
  Assert( this->isDataOwner() && v.isDataOwner(),
     ExcMessage("Trying to swap data of Vectors, that are not its owners") );

  std::swap (dim,    v.dim);
  std::swap (maxdim, v.maxdim);
  std::swap (val,    v.val);
}



template <typename Number>
bool
Vector<Number>::all_zero () const
{
  Assert (dim!=0, ExcEmptyObject());

  const_iterator p = begin(),
                 e = end();
  while (p!=e)
    if (*p++ != 0.0)
      return false;
  return true;
}



template <typename Number>
bool
Vector<Number>::is_non_negative () const
{
  Assert (dim!=0, ExcEmptyObject());

  const_iterator p = begin(),
                 e = end();
  while (p!=e)
    if (*p++ < 0.0)
      return false;
  return true;
}



template <typename Number>
template <typename Number2>
Number Vector<Number>::operator * (const Vector<Number2>& v) const
{
  Assert (dim!=0, ExcEmptyObject());

  if (this == reinterpret_cast<const Vector<Number>*>(&v))
    return norm_sqr();

  Assert (dim == v.size(), ExcDimensionMismatch(dim, v.size()));

  Number sum0 = 0,
         sum1 = 0,
         sum2 = 0,
         sum3 = 0;

  // use modern processors better by
  // allowing pipelined commands to be
  // executed in parallel
  const_iterator ptr  = begin(),
                 eptr = ptr + (dim/4)*4;
  typename Vector<Number2>::const_iterator vptr = v.begin();
  while (ptr!=eptr)
  {
    sum0 += (*ptr++ * *vptr++);
    sum1 += (*ptr++ * *vptr++);
    sum2 += (*ptr++ * *vptr++);
    sum3 += (*ptr++ * *vptr++);
  };
//   add up remaining elements
  while (ptr != end())
    sum0 += *ptr++ * *vptr++;

  return sum0+sum1+sum2+sum3;
}


template <typename Number>
Number Vector<Number>::norm_sqr () const
{
  Assert (dim!=0, ExcEmptyObject());

  Number sum0 = 0,
         sum1 = 0,
         sum2 = 0,
         sum3 = 0;

  // use modern processors better by
  // allowing pipelined commands to be
  // executed in parallel
  const_iterator ptr  = begin(),
                 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
  {
    sum0 += internal::VectorHelper::sqr(*ptr++);
    sum1 += internal::VectorHelper::sqr(*ptr++);
    sum2 += internal::VectorHelper::sqr(*ptr++);
    sum3 += internal::VectorHelper::sqr(*ptr++);
  };
  // add up remaining elements
  while (ptr != end())
    sum0 += internal::VectorHelper::sqr(*ptr++);

  return sum0+sum1+sum2+sum3;
}

template <typename Number>
Number Vector<Number>::rms_norm( const Vector<Number>& w ) const
{
  Assert( dim != 0, ExcEmptyObject() );
  Assert( dim == w.dim, ExcDimensionMismatch( dim, w.dim ) );

  Number sum0 = 0,
         sum1 = 0,
         sum2 = 0,
         sum3 = 0;

  // use modern processors better by
  // allowing pipelined commands to be
  // executed in parallel
  const_iterator ptr = begin(),
                 eptr = ptr + (dim/4)*4;
  const_iterator wptr = w.begin();

  while( ptr != eptr )
  {
    sum0 += internal::VectorHelper::sqr(*ptr++ * *wptr++);
    sum1 += internal::VectorHelper::sqr(*ptr++ * *wptr++);
    sum2 += internal::VectorHelper::sqr(*ptr++ * *wptr++);
    sum3 += internal::VectorHelper::sqr(*ptr++ * *wptr++);
  };
  // add up remaining elements
  while (ptr != end())
    sum0 += internal::VectorHelper::sqr(*ptr++ * *wptr++);

  return std::sqrt( (sum0+sum1+sum2+sum3) / size() );
}


template <typename Number>
Number Vector<Number>::mean_value () const
{
  Assert (dim!=0, ExcEmptyObject());

  Number sum0 = 0,
         sum1 = 0,
         sum2 = 0,
         sum3 = 0;

  // use modern processors better by
  // allowing pipelined commands to be
  // executed in parallel
  const_iterator ptr  = begin(),
                 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
  {
    sum0 += *ptr++;
    sum1 += *ptr++;
    sum2 += *ptr++;
    sum3 += *ptr++;
  };
  // add up remaining elements
  while (ptr != end())
    sum0 += *ptr++;

  return (sum0+sum1+sum2+sum3)/size();
}



template <typename Number>
Number Vector<Number>::l1_norm () const
{
  Assert (dim!=0, ExcEmptyObject());

  Number sum0 = 0,
         sum1 = 0,
         sum2 = 0,
         sum3 = 0;

  // use modern processors better by
  // allowing pipelined commands to be
  // executed in parallel
  const_iterator ptr  = begin(),
                 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
  {
    sum0 += std::fabs(*ptr++);
    sum1 += std::fabs(*ptr++);
    sum2 += std::fabs(*ptr++);
    sum3 += std::fabs(*ptr++);
  };
  // add up remaining elements
  while (ptr != end())
    sum0 += std::fabs(*ptr++);

  return sum0+sum1+sum2+sum3;
}


template <typename Number>
Number Vector<Number>::l2_norm () const
{
  return std::sqrt(norm_sqr());
}


template <typename Number>
Number Vector<Number>::lp_norm (const Number p) const
{
  Assert (dim!=0, ExcEmptyObject());

  Number sum0 = 0,
         sum1 = 0,
         sum2 = 0,
         sum3 = 0;

  // use modern processors better by
  // allowing pipelined commands to be
  // executed in parallel
  const_iterator ptr  = begin(),
                 eptr = ptr + (dim/4)*4;
  while (ptr!=eptr)
  {
    sum0 += std::pow(std::fabs(*ptr++), p);
    sum1 += std::pow(std::fabs(*ptr++), p);
    sum2 += std::pow(std::fabs(*ptr++), p);
    sum3 += std::pow(std::fabs(*ptr++), p);
  };
  // add up remaining elements
  while (ptr != end())
    sum0 += std::pow(std::fabs(*ptr++), p);

  return std::pow(sum0+sum1+sum2+sum3,
      static_cast<Number>(1./p));
}


template <typename Number>
Number Vector<Number>::linfty_norm () const
{
  Assert (dim!=0, ExcEmptyObject());

  Number max0=0.,
         max1=0.,
         max2=0.,
         max3=0.;
  for (unsigned int i=0; i<(dim/4); ++i) 
  {
    if (max0<std::fabs(val[4*i]))   max0=std::fabs(val[4*i]);
    if (max1<std::fabs(val[4*i+1])) max1=std::fabs(val[4*i+1]);
    if (max2<std::fabs(val[4*i+2])) max2=std::fabs(val[4*i+2]);
    if (max3<std::fabs(val[4*i+3])) max3=std::fabs(val[4*i+3]);
  };
  // add up remaining elements
  for (unsigned int i=(dim/4)*4; i<dim; ++i)
    if (max0<std::fabs(val[i]))
      max0 = std::fabs(val[i]);

  return std::max (std::max(max0, max1),
      std::max(max2, max3));
}


  template <typename Number>
Vector<Number>& Vector<Number>::operator += (const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyObject());

  add (v);
  return *this;
}


  template <typename Number>
Vector<Number>& Vector<Number>::operator -= (const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));

  iterator i_ptr = begin(),
           i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ -= *v_ptr++;

  return *this;
}


  template <typename Number>
void Vector<Number>::add (const Number v)
{
  Assert (dim!=0, ExcEmptyObject());

  iterator i_ptr = begin(),
           i_end = end();
  while (i_ptr!=i_end)
    *i_ptr++ += v;
}


  template <typename Number>
void Vector<Number>::add (const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));

  iterator i_ptr = begin(),
           i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += *v_ptr++;
}


  template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));

  iterator i_ptr = begin(),
           i_end = end();
  const_iterator v_ptr = v.begin();
  while (i_ptr!=i_end)
    *i_ptr++ += a * *v_ptr++;
}


  template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number>& v,
    const Number b, const Vector<Number>& w)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionMismatch(dim, w.dim));

  iterator i_ptr = begin();
  const_iterator v_ptr = v.begin(),
                 w_ptr = w.begin();
  int N = dim;
  for( int i = 0; i < N; i++)
    i_ptr[i] += a*v_ptr[i] + b*w_ptr[i];
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  
  iterator i_ptr = begin(),
           i_end = end();
  const_iterator v_ptr = v.begin();
  // int N = dim;
  for (; i_ptr!=i_end; ++i_ptr)
    *i_ptr = x * *i_ptr  + *v_ptr++;
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
    const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  
  iterator i_ptr = begin();
  const_iterator v_ptr = v.begin();
  int N = dim;
  for( int i = 0; i < N; i++)
    i_ptr[i] = x*i_ptr[i] + a*v_ptr[i];
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
    const Vector<Number>& v, const Number b,
    const Vector<Number>& w)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionMismatch(dim, w.dim));

  iterator i_ptr = begin();
  const_iterator v_ptr = v.begin(),
                 w_ptr = w.begin();
  int N = dim;
  for( int i = 0; i < N; i++ )
    i_ptr[i] = x*i_ptr[i] + a*v_ptr[i] + b*w_ptr[i];
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
    const Vector<Number>& v, const Number b,
    const Vector<Number>& w, const Number c,
    const Vector<Number>& y)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));
  Assert (dim == w.dim, ExcDimensionMismatch(dim, w.dim));
  Assert (dim == y.dim, ExcDimensionMismatch(dim, y.dim));

  iterator i_ptr = begin();
  const_iterator v_ptr = v.begin(),
                 w_ptr = w.begin(),
                 y_ptr = y.begin();
  int N = dim;
  for( int i = 0; i<N; i++ )
    i_ptr[i] = x*i_ptr[i] + a*v_ptr[i] + b*w_ptr[i] + c*y_ptr[i];
}



  template <typename Number>
void Vector<Number>::scale (const Number factor)
{
  Assert (dim!=0, ExcEmptyObject());

  iterator             ptr  = begin();
  const const_iterator eptr = end();
  while (ptr!=eptr)
    *ptr++ *= factor;
}



template <typename Number>
  template <typename Number2>
void Vector<Number>::scale (const Vector<Number2> &s)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == s.dim, ExcDimensionMismatch(dim, s.dim));

  iterator             ptr  = begin();
  const const_iterator eptr = end();
  typename Vector<Number2>::const_iterator sptr = s.begin();
  while (ptr!=eptr)
    *ptr++ *= *sptr++;
}


  template <typename Number>
void Vector<Number>::invabsequ( const Number a, const Vector<Number>& v, const Number b )
{
  Assert( dim != 0, ExcEmptyObject() );
  Assert( dim == v.dim, ExcDimensionMismatch( dim, v.dim) );

  iterator i_ptr = begin();
  const_iterator v_ptr = v.begin();
  int N = dim;
  for( int i = 0; i < N; i++ )
    i_ptr[i] = 1.0 / ( a * std::fabs( v_ptr[i] ) + b );
}

  template <typename Number>
void Vector<Number>::maxabs( const Vector<Number>& v1, const Vector<Number>& v2)
{
  Assert( dim != 0, ExcEmptyObject() );
  Assert( dim == v1.dim, ExcDimensionMismatch( dim, v1.dim) );
  Assert( dim == v2.dim, ExcDimensionMismatch( dim, v2.dim) );

  iterator i_ptr = begin(),
           i_end = end();
  const_iterator v1_ptr = v1.begin();
  const_iterator v2_ptr = v2.begin();

  while( i_ptr != i_end )
    *i_ptr++ = std::max( std::fabs( *v1_ptr++ ), std::fabs( *v2_ptr++ ) );
}


template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number>& u,
    const Number b, const Vector<Number>& v)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == u.dim, ExcDimensionMismatch(dim, u.dim));
  Assert (dim == v.dim, ExcDimensionMismatch(dim, v.dim));

  iterator i_ptr = begin();
  const_iterator u_ptr = u.begin(),
                 v_ptr = v.begin();
  int N = dim;
  for( int i = 0; i < N; i++ )
    i_ptr[i] = a*u_ptr[i] + b*v_ptr[i];
}


template <typename Number>
  template <typename Number2>
void Vector<Number>::equ (const Number a, const Vector<Number2>& u)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == u.dim, ExcDimensionMismatch(dim, u.dim));
  iterator i_ptr = begin(),
           i_end = end();
  typename Vector<Number2>::const_iterator u_ptr = u.begin();
  while (i_ptr!=i_end)
    *i_ptr++ = a * *u_ptr++;
}



template <typename Number>
void Vector<Number>::ratio (const Vector<Number> &a, const Vector<Number> &b)
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (a.dim == b.dim, ExcDimensionMismatch (a.dim, b.dim));
  Assert (dim == a.dim, ExcDimensionMismatch (dim, a.dim));

  iterator i_ptr = begin();
  const_iterator a_ptr = a.begin(),
                 b_ptr = b.begin();
  int N = dim;
  for( int i = 0; i < N; i++ )
    i_ptr[i] = a_ptr[i] / b_ptr[i];
}



template <typename Number>
  Vector<Number>&
Vector<Number>::operator = (const Vector<Number>& v)
{
  Assert (dim == v.dim, ExcDimensionMismatch (dim, v.dim));
  if (dim!=0)
    std::copy (v.begin(), v.end(), begin());

  return *this;
}



template <typename Number>
template <typename Number2>
  Vector<Number>&
Vector<Number>::operator = (const Vector<Number2>& v)
{
  Assert (dim == v.dim, ExcDimensionMismatch (dim, v.dim));
  if (dim!=0)
    std::copy (v.begin(), v.end(), begin());

  return *this;
}



template <typename Number>
template <typename Number2>
bool
Vector<Number>::operator == (const Vector<Number2>& v) const
{
  Assert (dim!=0, ExcEmptyObject());
  Assert (dim == v.size(), ExcDimensionMismatch(dim, v.size()));

  for (unsigned int i=0; i<dim; ++i)
    if (val[i] != v.val[i])
      return false;

  return true;
}



template <typename Number>
void Vector<Number>::print (const char* format) const
{
  Assert (dim!=0, ExcEmptyObject());
  if (!format) format = " %5.2f";
  for (unsigned int j=0;j<size();j++)
    std::printf (format, val[j]);
  std::printf ("\n");
}



template <typename Number>
void Vector<Number>::print (std::ostream      &out,
    const unsigned int precision,
    const bool         scientific,
    const bool         across) const
{
  Assert (dim!=0, ExcEmptyObject());
  AssertThrow (out, ExcIO());

  std::ios::fmtflags old_flags = out.flags();
  unsigned int old_precision = out.precision (precision);

  out.precision (precision);
  if (scientific)
    out.setf (std::ios::scientific, std::ios::floatfield);
  else
    out.setf (std::ios::fixed, std::ios::floatfield);

  if (across)
    for (unsigned int i=0; i<size(); ++i)
      out << static_cast<double>(val[i]) << ' ';
  else
    for (unsigned int i=0; i<size(); ++i)
      out << static_cast<double>(val[i]) << std::endl;
  out << std::endl;

  AssertThrow (out, ExcIO());
  // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
}


#endif
