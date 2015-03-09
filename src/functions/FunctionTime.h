#ifndef FUNCTION_TIME_H
#define FUNCTION_TIME_H

#include <sundials/sundials_types.h>

class FunctionTime
{
  public:
  
    FunctionTime (const realtype initial_time = 0.0);

    virtual ~FunctionTime();
  
    realtype getTime () const;

    virtual void setTime (const realtype new_time);

    virtual void advanceTime (const realtype delta_t);

  private:

    realtype time;
};

inline realtype
FunctionTime::getTime () const
{
  return time;
}


#endif // FUNCTION_TIME_H

