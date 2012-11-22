#include "ProgressDisplay.h"

ProgressDisplay :: ProgressDisplay( unsigned long expected_count,
                                    std::ostream & os,
                                    const std::string & s1, //leading strings
                                    const std::string & s2,
                                    const std::string & s3 )
        : m_os(os), m_s1(s1), m_s2(s2), m_s3(s3)
{
    restart(expected_count);
}

//  Effects: display appropriate scale
//  Postconditions: count()==0, expected_count()==expected_count
void ProgressDisplay :: restart( unsigned long expected_count )
{
    _count = _next_tic_count = _tic = 0;
    _expected_count = expected_count;

    m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
    << m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
    << std::endl  // endl implies flush, which ensures display
    << m_s3;

    if ( !_expected_count )
    {
        _expected_count = 1;  // prevent divide by zero
    }
} // restart

//  Effects: Display appropriate progress tic if needed.
//  Postconditions: count()== original count() + increment
//  Returns: count().
unsigned long  ProgressDisplay :: operator+=( unsigned long increment )
{
    if ((_count += increment) >= _next_tic_count )
    {
        display_tic();
    }
    return _count;
}

// use of floating point ensures that both large and small counts
// work correctly.  static_cast<>() is also used several places
// to suppress spurious compiler warnings.
void ProgressDisplay :: display_tic()
{
    unsigned int tics_needed = static_cast<unsigned int>( (static_cast<double>(_count)/_expected_count)*50.0 );

    do
    {
        m_os << '*' << std::flush;
    }
    while (++_tic < tics_needed);

    _next_tic_count = static_cast<unsigned long>((_tic/50.0)*_expected_count);

    if ( _count == _expected_count )
    {
        if ( _tic < 51 )
        {
            m_os << '*';
        }
        m_os << std::endl;
    }
} // display_tic

