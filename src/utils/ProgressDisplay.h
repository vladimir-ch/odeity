#ifndef PROGRESS_DISPLAY_H
#define PROGRESS_DISPLAY_H

#include <iostream>           // for ostream, cout, etc
#include <string>             // for string

class ProgressDisplay
{
private:
    std::ostream &     m_os;  // may not be present in all imps
    const std::string  m_s1;  // string is more general, safer than
    const std::string  m_s2;  //  const char *, and efficiency or size are
    const std::string  m_s3;  //  not issues

    unsigned long _count, _expected_count, _next_tic_count;
    unsigned int _tic;

public:
    ProgressDisplay( unsigned long expected_count,
                     std::ostream & os = std::cout,
                     const std::string & s1 = "\n", //leading strings
                     const std::string & s2 = "",
                     const std::string & s3 = "" );

    //  Effects: display appropriate scale
    //  Postconditions: count()==0, expected_count()==expected_count
    void restart( unsigned long expected_count );

    //  Effects: Display appropriate progress tic if needed.
    //  Postconditions: count()== original count() + increment
    //  Returns: count().
    unsigned long  operator+=( unsigned long increment );
    unsigned long  operator++()
    {
        return operator+=( 1 );
    }
    unsigned long  count() const
    {
        return _count;
    }
    unsigned long  expected_count() const
    {
        return _expected_count;
    }

private:
    void display_tic();
}; // ProgressDisplay

#endif // PROGRESS_DISPLAY_H

