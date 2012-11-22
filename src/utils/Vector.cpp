//---------------------------------------------------------------------------
//    vector.cc,v 1.24 2005/03/29 00:47:24 guido Exp
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

#include "Vector.templates.h"

template class Vector<double>;
template Vector<double>::Vector (const Vector<float> &);
template Vector<double>& Vector<double>::operator=<float>(const Vector<float>&);
template bool Vector<double>::operator==<double>(const Vector<double>&) const;
template bool Vector<double>::operator==<float>(const Vector<float>&) const;
template double Vector<double>::operator*<float>(const Vector<float>&) const;
template double Vector<double>::operator*<double>(const Vector<double>&) const;
template void Vector<double>::equ<double>(const double, const Vector<double>&);
template void Vector<double>::equ<float>(const double, const Vector<float>&);
template void Vector<double>::scale<double>(const Vector<double>&);
template void Vector<double>::scale<float>(const Vector<float>&);

template class Vector<float>;
template Vector<float>::Vector (const Vector<double> &);
template Vector<float>& Vector<float>::operator=<double>(const Vector<double>&);
template bool Vector<float>::operator==<double>(const Vector<double>&) const;
template bool Vector<float>::operator==<float>(const Vector<float>&) const;
template float Vector<float>::operator*<float>(const Vector<float>&) const;
template float Vector<float>::operator*<double>(const Vector<double>&) const;
template void Vector<float>::equ<double>(const float, const Vector<double>&);
template void Vector<float>::equ<float>(const float, const Vector<float>&);
template void Vector<float>::scale<double>(const Vector<double>&);
template void Vector<float>::scale<float>(const Vector<float>&);
