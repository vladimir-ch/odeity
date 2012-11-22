//---------------------------------------------------------------------------
//    job_identifier.cc,v 1.20 2005/03/29 00:49:28 guido Exp
//    Version: Version-5-2-0
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include "JobIdentifier.h"
#include "Utilities.h"

#include <ctime>
#include <sys/time.h>
#include <unistd.h>

JobIdentifier jobid;

JobIdentifier::JobIdentifier()
{
    using namespace Utilities;

    time_t curTime = std::time( 0 );
    tm * t = localtime( &curTime );

    char name[100];
    gethostname( name, 99 );

    // computername-date-time;
    id_ += std::string( name ) + std::string("-");
    id_ += toString( t->tm_year + 1900 ) + intToString( t->tm_mon+1, 2 ) + intToString( t->tm_mday, 2 ) + "-";
    id_ += intToString( t->tm_hour, 2 ) + intToString( t->tm_min, 2 ) + intToString( t->tm_sec, 2 );
}



const std::string
JobIdentifier::operator ()() const
{
    return id_;
}


