//---------------------------------------------------------------------------
//utilities.h,v 1.5 2005/09/08 23:56:45 wolf Exp Version: Version-5-2-0
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed without copyright
//    and license information. Please refer to the file deal.II/doc/license.html
//    for the  text  and further information on this license.
//
//---------------------------------------------------------------------------
#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <string>
#include <sstream>

/** A namespace for utility functions that are not particularly specific to
 * finite element computing or numerical programs, but nevertheless are needed
 * in various contexts when writing applications.
 *
 * @author Wolfgang Bangerth, 2005
 */
namespace Utilities
{

    /** Convert a number @p i to a string, with as many digits as given to fill
     * with leading zeros.
     */
    std::string intToString (const unsigned int i, const unsigned int digits);

    /** Determine how many digits are needed to represent numbers at most as large
     * as the given number.
     */
    unsigned int neededDigits (const unsigned int max_number);

    /** Given a string, convert it to an integer. Throw an assertion if that is
     * not possible.
     */
    int stringToInt (const std::string &s);


    /** Convert a string to a numeric type. Return false is the conversion is not
     * possible
     */
    template <class T> bool fromString(
            T& t,
            const std::string& s,
            std::ios_base& (*f)(std::ios_base&) )
    {
        std::istringstream iss(s);
        return !(iss >> f >> t).fail();
    }

    template<typename T>
    std::string toString(const T& t)
    {
        std::ostringstream s;
        s << t;

        return s.str();
    }

    /** Given a list of strings, convert it to a list of integers. Throw an
     * assertion if that is not possible.
     */
    std::vector<int> stringToInt (const std::vector<std::string> &s);

    /** Given a string that contains text separated by a @p delimiter, split it
     * into its components; for each component, remove leading and trailing
     * spaces.
     *
     * The default value of the delimiter is a comma, so that the function
     * splits comma separated lists of strings.
     */
    std::vector<std::string> splitStringList (
            const std::string &s,
            const char delimiter = ',' );

    /** Take a text, usually a documentation or something, and try to break it
     * into individual lines of text at most @p width characters wide, by
     * breaking at spaces in the text. If this is not possible, return the
     * shortest lines than are longer than @p width.
     */
    std::vector<std::string> breakTextIntoLines (
            const std::string &original_text,
            const unsigned int width );

    /** Generate a random number from a normalized Gaussian probability
     * distribution centered around @p a and with standard deviation @p sigma.
     */
    double generateNormalRandomNumber (const double a, const double sigma);


    /** A namespace for utility functions that probe system properties.
     */
    namespace System
    {

        /** Return the name of the host this process runs on.
         */
        std::string getHostname ();


        /** Return the present time as HH:MM:SS.
        */
        std::string getTime ();

        /** Return true, if the directory @p dirName exists */
        bool dirExists( std::string const & dirName );
    }
}

#endif
