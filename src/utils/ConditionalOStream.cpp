//---------------------------------------------------------------------------
//    conditional_ostream.cc,v 1.5 2005/03/29 00:49:28 guido Exp
//    Version: Version-5-2-0
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include "ConditionalOStream.h"


ConditionalOStream::ConditionalOStream(std::ostream &stream,
                                       const bool    active)
                :
    outputStream_(stream),
    activeFlag_(active)
{}


void ConditionalOStream::setCondition(bool flag)
{
  activeFlag_ = flag;
}


bool ConditionalOStream::isActive() const
{
  return activeFlag_;
}
