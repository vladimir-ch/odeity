//---------------------------------------------------------------------------
//      utilities.cc,v 1.7 2005/09/12 01:20:51 wolf Exp   
//    Version: Version-5-2-0
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include "Utilities.h"
#include "Exceptions.h"

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cerrno>
#include <cmath>
#include <unistd.h>
#include <sys/types.h>
#include <sstream>
#include <limits>


#include <sys/stat.h>
#include <dirent.h>

namespace Utilities
{

  DeclException2 (ExcInvalidNumber2StringConversersion,
      unsigned int, unsigned int,
      << "When trying to convert " << arg1
      << " to a string with " << arg2 << " digits");
  DeclException1 (ExcInvalidNumber,
      unsigned int,
      << "Invalid number " << arg1);
  DeclException1 (ExcCantConvertString,
      std::string,
      << "Can't convert the string " << arg1
                  << " to the desired type");

  std::string
  intToString (const unsigned int i,
               const unsigned int digits)
  {
    AssertThrow ( ! ((digits==1 && i>=10)   ||
        (digits==2 && i>=100)  ||
        (digits==3 && i>=1000) ||
        (digits==4 && i>=10000)||
        (i>=100000)),
        ExcInvalidNumber2StringConversersion(i, digits));

    std::string s;
    switch (digits) 
      {
        case 4:
          s += '0' + i/1000;
        case 3:
          s += '0' + (i%1000)/100;
        case 2:
          s += '0' + (i%100)/10;
        case 1:
          s += '0' + i%10;
          break;
        default:
          s += "invalid digits information";
      }
    return s;
  }

  unsigned int
  neededDigits (const unsigned int max_number)
  {
    if (max_number < 10)
      return 1;
    if (max_number < 100)
      return 2;
    if (max_number < 1000)
      return 3;
    if (max_number < 10000)
      return 4;
    if (max_number < 100000)
      return 5;
    if (max_number < 1000000)
      return 6;
    AssertThrow (false, ExcInvalidNumber(max_number));
    return 0;
  }


  int
  stringToInt (const std::string &s)
  {
    std::istringstream ss(s);
    static const int max_int = std::numeric_limits<int>::max();
    int i = max_int;
    ss >> i;

                                    // check for errors
    AssertThrow (i != max_int, ExcCantConvertString (s));

    return i;
  }

  std::vector<int>
  stringToInt (const std::vector<std::string> &s)
  {
    std::vector<int> tmp (s.size());
    for (unsigned int i=0; i<s.size(); ++i)
      tmp[i] = stringToInt (s[i]);
    return tmp;
  }

  std::vector<std::string>
  splitStringList (const std::string &s,
                   const char delimiter)
  {
    std::string tmp = s;
    std::vector<std::string> split_list;
    split_list.reserve (std::count (tmp.begin(), tmp.end(), delimiter)+1);

            // split the input list
    while (tmp.length() != 0)
      {
        std::string name;
        name = tmp;

        if (name.find(delimiter) != std::string::npos)
          {
            name.erase (name.find(delimiter), std::string::npos);
            tmp.erase (0, tmp.find(delimiter)+1);
          }
        else
          tmp = "";

        while ((name.length() != 0) && (name[0] == ' '))
          name.erase (0,1);

        while (name[name.length()-1] == ' ')
          name.erase (name.length()-1, 1);

        split_list.push_back (name);
      }

    return split_list;
  }

  std::vector<std::string>
  breakTextIntoLines (const std::string &original_text,
                      const unsigned int width)
  {
    std::string text = original_text;
    std::vector<std::string> lines;

                                    // remove trailing spaces
    while ((text.length() != 0) && (text[text.length()-1] == ' '))
      text.erase(text.length()-1,1);

                                    // then split the text into lines
    while (text.length() != 0)
      {
                                    // in each iteration, first remove
                                    // leading spaces
        while ((text.length() != 0) && (text[0] == ' '))
          text.erase(0, 1);

                                    // if we can fit everything into one
                                    // line, then do so. otherwise, we have
                                    // to keep breaking
        if (text.length() < width)
          {
            lines.push_back (text);
            text = "";
          }
        else
          {
                                    // starting at position width, find the
                                    // location of the previous space, so
                                    // that we can break around there
            int location = std::min<int>(width,text.length()-1);
            for (; location>=0; --location)
              if (text[location] == ' ')
                break;

                                    // if there are no spaces, then try if
                                    // there are spaces coming up
            if (location == 0)
              for (location = std::min<int>(width,text.length()-1);
                  location<static_cast<int>(text.length());
                  ++location)
                if (text[location] == ' ')
                  break;

                                    // now take the text up to the found
                                    // location and put it into a single
                                    // line, and remove it from 'text'
            lines.push_back (std::string (text, 0, location));
            text.erase (0, location);
          }
      }

    return lines;
  }

  double
  generateNormalRandomNumber (const double a,
                              const double sigma)
  {
            // if no noise: return now
    if (sigma == 0.0)
      return a;

    static unsigned int seed = 0xabcd1234;
    const double y = 1.0*rand_r(&seed)/RAND_MAX;

            // find x such that y=erf(x). do so
            // using a Newton method to find
            // the zero of F(x)=erf(x)-y. start
            // at x=0
    double x = 0;
    unsigned int iteration = 0;
    while (true)
      {
        const double residual = 0.5+erf(x/std::sqrt(2.)/sigma)/2-y;
        if (std::fabs(residual) < 1e-7)
          break;
        const double F_prime = 1./std::sqrt(2*3.1415926536)/sigma *
            std::exp(-x*x/sigma/sigma/2);
        x += -residual / F_prime;

          // make sure that we don't
          // recurse endlessly
        ++iteration;
        Assert (iteration < 20, ExcInternalError());
      }
    return x+a;
  }

  namespace System
  {
    std::string getHostname ()
    {
      const unsigned int N=1024;
      char hostname[N];
      gethostname (&(hostname[0]), N-1);
      return hostname;
    }

    std::string getTime ()
    {
      std::time_t  time1= std::time (0);
      std::tm     *time = std::localtime(&time1); 

      std::ostringstream o;
      o << time->tm_hour << ":"
        << (time->tm_min < 10 ? "0" : "") << time->tm_min << ":"
        << (time->tm_sec < 10 ? "0" : "") << time->tm_sec;
      return o.str();
    }

    bool
    dirExists( std::string const & dirName )
    {
        DIR* dirStream = opendir( dirName.c_str() );
        if ( 0 == dirStream )
            return false;
        closedir( dirStream );
        return true;
    }
  }
}