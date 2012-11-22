//---------------------------------------------------------------------------
//    exceptions.cc,v 1.44 2005/08/25 07:26:17 hartmann Exp
//    Version: Version-5-2-0
//
//    Copyright (C) 1998, 2000, 2001, 2002, 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include "Exceptions.h"

#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <execinfo.h>

namespace exceptions
{
  std::string additionalAssertOutput;

  void setAdditionalAssertOutput (const char * const p)
  {
    additionalAssertOutput = p;
  }

  bool showStacktrace = true;

  void suppressStacktraceInExceptions ()
  {
    showStacktrace = false;
  }

  bool abortOnException = true;

  void disableAbortOnException ()
  {
    abortOnException = false;
  }
}


ExceptionBase::ExceptionBase ()
                :
    file_(""), line_(0), function_(""), cond_(""), exc_(""),
    stacktrace_ (0),
    nStacktraceFrames_ (0)
{}

ExceptionBase::ExceptionBase (const char* f, const int l, const char *func,
            const char* c, const char *e)
                :
    file_(f), line_(l), function_(func), cond_(c), exc_(e),
    stacktrace_ (0),
    nStacktraceFrames_ (0)
{}

ExceptionBase::ExceptionBase (const ExceptionBase &exc)
                :
    std::exception (exc),
    file_(exc.file_), line_(exc.line_),
    function_(exc.function_), cond_(exc.cond_), exc_(exc.exc_),
    // don't copy stacktrace to
    // avoid double de-allocation
    // problem
    stacktrace_ (0),
    nStacktraceFrames_ (0)
{}

ExceptionBase::~ExceptionBase () throw ()
{
  if( stacktrace_ != 0 )
  {
    free( stacktrace_ );
    stacktrace_ = 0;
  }
}

void ExceptionBase::setFields (const char* f,
        const int l,
        const char *func,
        const char *c,
        const char *e)
{
  file_ = f;
  line_ = l;
  function_ = func;
  cond_ = c;
  exc_  = e;

  // get a stacktrace how we got here
  void * array[25];
  nStacktraceFrames_ = backtrace(array, 25);
  stacktrace_ = backtrace_symbols(array, nStacktraceFrames_);
}

void ExceptionBase::printStackTrace (std::ostream &out) const
{
  if( nStacktraceFrames_ == 0 )
    return;

  if( exceptions::showStacktrace == false )
    return;


  // if there is a stackframe stored, print it
  out << std::endl;
  out << "Stacktrace:" << std::endl
      << "-----------" << std::endl;

  // print the stacktrace. first
  // omit all those frames that have
  // ExceptionBase or
  // deal_II_exceptions in their
  // names, as these correspond to
  // the exception raising mechanism
  // themselves, rather than the
  // place where the exception was
  // triggered
  int frame = 0;
  while ((frame < nStacktraceFrames_)
    &&
    ((std::string(stacktrace_[frame]).find("ExceptionBase") != std::string::npos)
    ||
    (std::string(stacktrace_[frame]).find("deal_II_exceptions") != std::string::npos)))
    ++frame;

  // output the rest
  const unsigned int first_significant_frame = frame;
  for (; frame < nStacktraceFrames_; ++frame)
    out << '#' << frame - first_significant_frame
        << "  "
        << stacktrace_[frame]
        << std::endl;
}

void ExceptionBase::printExcData (std::ostream &out) const
{
  out << "An error occurred in line <" << line_
      << "> of file <" << file_
      << "> in function" << std::endl
      << "    " << function_ << std::endl
      << "The violated condition was: "<< std::endl
      << "    " << cond_ << std::endl
      << "The name and call sequence of the exception was:" << std::endl
      << "    " << exc_  << std::endl
      << "Additional Information: " << std::endl;
}

void ExceptionBase::printInfo (std::ostream &out) const
{
  out << "(none)" << std::endl;
}


const char * ExceptionBase::what () const throw ()
{
  // if we say that this function
  // does not throw exceptions, we
  // better make sure it does not
  try
    {
      // have a place where to store the
      // description of the exception as
      // a char *
      //
      // this thing obviously is not
      // multi-threading safe, but we
      // don't care about that for now
      //
      // we need to make this object
      // static, since we want to return
      // the data stored in it and
      // therefore need a lifetime which
      // is longer than the execution
      // time of this function
      static std::string description;
      // convert the messages printed by
      // the exceptions into a
      // std::string
      std::ostringstream converter;

      converter << "--------------------------------------------------------"
                << std::endl;
      // put general info into the std::string
      printExcData (converter);
      // put in exception specific data
      printInfo (converter);
      printStackTrace (converter);
      converter << "--------------------------------------------------------"
                << std::endl;

      description = converter.str();

      return description.c_str();
    }
  catch (std::exception &exc) 
    {
      std::cerr << "*** Exception encountered in exception handling routines ***"
                << std::endl
                << "*** Message is "   << std::endl
                << exc.what ()         << std::endl
                << "*** Aborting! ***" << std::endl;

      std::abort ();
    }
  catch (...)
    {
      std::cerr << "*** Exception encountered in exception handling routines ***"
                << std::endl
                << "*** Aborting! ***" << std::endl;

      std::abort ();
    }
  return 0;
}

namespace exceptions
{
  namespace internals
  {

            /**
              * Number of exceptions dealt
              * with so far. Zero at program
              * start. Messages are only
              * displayed if the value is
              * zero.
              */
    unsigned int nTreatedExceptions;
    ExceptionBase *lastException;

    void issueErrorAssert (const char *file,
          int         line,
          const char *function,
          const char *cond,
          const char *exc_name,
          ExceptionBase &e)
    {
              // fill the fields of the
              // exception object
      e.setFields (file, line, function, cond, exc_name);

              // if no other exception has
              // been displayed before, show
              // this one
      if (nTreatedExceptions == 0)
      {
          std::cerr << "--------------------------------------------------------"
                    << std::endl;
            // print out general data
          e.printExcData (std::cerr);
            // print out exception
            // specific data
          e.printInfo (std::cerr);
          e.printStackTrace (std::cerr);
          std::cerr << "--------------------------------------------------------"
                    << std::endl;

            // if there is more to say,
            // do so
          if (!additionalAssertOutput.empty())
            std::cerr << additionalAssertOutput << std::endl;
      }
      else
      {
            // if this is the first
            // follow-up message,
            // display a message that
            // further exceptions are
            // suppressed
        if (nTreatedExceptions == 1)
          std::cerr << "******** More assertions fail but messages are suppressed! ********"
                    << std::endl;
      };

              // increase number of treated
              // exceptions by one
      nTreatedExceptions++;
      lastException = &e;

              // abort the program now since
              // something has gone horribly
              // wrong. however, there is one
              // case where we do not want to
              // do that, namely when another
              // exception, possibly thrown
              // by AssertThrow is active,
              // since in that case we will
              // not come to see the original
              // exception. in that case
              // indicate that the program is
              // not aborted due to this
              // reason.
      if (std::uncaught_exception() == true)
      {
            // only display message once
        if (nTreatedExceptions <= 1)
          std::cerr << "******** Program is not aborted since another exception is active! ********"
                    << std::endl;
      }
      else if(exceptions::abortOnException == true)
        std::abort ();
    }

    void abort ()
    {
      if (exceptions::abortOnException == true)
        std::abort ();
    }
  }
}


// from the aclocal file:
// Newer versions of gcc have a very nice feature: you can set
// a verbose terminate handler, that not only aborts a program
// when an exception is thrown and not caught somewhere, but
// before aborting it prints that an exception has been thrown,
// and possibly what the std::exception::what() function has to
// say. Since many people run into the trap of not having a
// catch clause in main(), they wonder where that abort may be
// coming from. The terminate handler then at least says what is
// missing in their program.
namespace __gnu_cxx
{
  extern void __verbose_terminate_handler ();
}

namespace
{
  struct preload_terminate_dummy
  {
    preload_terminate_dummy()
      { std::set_terminate(__gnu_cxx::__verbose_terminate_handler); }
  };

  static preload_terminate_dummy dummy;
}
