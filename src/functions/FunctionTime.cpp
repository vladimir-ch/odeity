#include "FunctionTime.h"


FunctionTime::FunctionTime(const realtype initial_time)	:
		time(initial_time)
{}



FunctionTime::~FunctionTime()
{}



void
FunctionTime::setTime (const realtype new_time)
{
  time = new_time;
}



void
FunctionTime::advanceTime (const realtype delta_t)
{
  setTime( time + delta_t );
}

