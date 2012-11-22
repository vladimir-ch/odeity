//---------------------------------------------------------------------------
//    log.cc,v 1.46 2005/04/08 09:16:43 guido Exp
//    Version: Version-5-2-0
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include "LogStream.h"
#include "JobIdentifier.h"

// include sys/resource.h for rusage(). Mac OS X needs sys/time.h then
// as well (strange), so include that, too.
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

LogStream logger;

LogStream::LogStream()
    :
        stdOut_(&std::cerr), file_(0), wasEndl_(true),
        stdDepth_(10000), fileDepth_(10000),
        printUtime_(false), diffUtime_(false),
        lastTime_ (0.), doubleThreshold_(0.), oldCerr_(0)
{
    prefixes_.push("ODEity:");
    stdOut_->setf(std::ios::showpoint | std::ios::left);
}

LogStream::~LogStream()
{
    if (oldCerr_)
        std::cerr.rdbuf(oldCerr_);
}

void LogStream::attach(std::ostream& o)
{
    file_ = &o;
    o.setf(std::ios::showpoint | std::ios::left);
    //  o << jobid();
}

void LogStream::detach ()
{
    file_ = 0;
}

void LogStream::logCerr ()
{
    if (oldCerr_ == 0)
    {
        oldCerr_ = std::cerr.rdbuf(file_->rdbuf());
    } else {
        std::cerr.rdbuf(oldCerr_);
        oldCerr_ = 0;
    }
}

std::ostream& LogStream::getConsole()
{
    return *stdOut_;
}

std::ostream& LogStream::getFileStream()
{
    Assert(file_, ExcNoFileStreamGiven());
    return *file_;
}

const std::string& LogStream::getPrefix() const
{
    return prefixes_.top();
}

void LogStream::push (const std::string& text)
{
    std::string pre=prefixes_.top();
    pre += text;
    pre += std::string(":");
    prefixes_.push(pre);
}

void LogStream::pop ()
{
    if (prefixes_.size() > 1)
        prefixes_.pop();
}

unsigned int LogStream::depthConsole (const unsigned n)
{
    const unsigned int h = stdDepth_;
    stdDepth_ = n;
    return h;
}

unsigned int LogStream::depthFile (const unsigned n)
{
    const unsigned int h = fileDepth_;
    fileDepth_ = n;
    return h;
}

void LogStream::thresholdDouble (const double t)
{
    doubleThreshold_ = t;
}

bool LogStream::logExecutionTime (const bool flag)
{
    const bool h = printUtime_;
    printUtime_ = flag;
    return h;
}

bool LogStream::logTimeDifferences (const bool flag)
{
    const bool h = diffUtime_;
    diffUtime_ = flag;
    return h;
}

void LogStream::printLineHead()
{
    rusage usage;
    double utime = 0.;
    if (printUtime_)
    {
        getrusage(RUSAGE_SELF, &usage);
        utime = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
        if (diffUtime_)
        {
            double diff = utime - lastTime_;
            lastTime_ = utime;
            utime = diff;
        }
    }

    const std::string& head = getPrefix();

    if (prefixes_.size() <= stdDepth_)
    {
        if (printUtime_)
        {
            int p = stdOut_->width(5);
            *stdOut_ << utime << ':';
            stdOut_->width(p);
        }
        *stdOut_ <<  head << ':';
    }

    if (file_ && (prefixes_.size() <= fileDepth_))
    {
        if (printUtime_)
        {
            int p = file_->width(6);
            *file_ << utime << ':';
            file_->width(p);
        }
        *file_ << head << ':';
    }
}
