#include "RungeKutta23.h"

RungeKutta23::RungeKutta23()
  : ExplicitRungeKuttaBase( 2, 4, true ) // 2th order, 4 stages, FSAL true
{
  solverName_ = "Runge-Kutta 3(2) (Bogacki-Shampine, also implemented in Matlab as ode23)";

  c_(0) = 0.0;
  c_(1) = 1.0/2.0;
  c_(2) = 3.0/4.0;
  c_(3) = 1.0;

  b_(0) = 2.0/9.0;
  b_(1) = 1.0/3.0; 
  b_(2) = 4.0/9.0;
  b_(3) = 0.0;

  e_(0) = -5.0/72.0;
  e_(1) =  1.0/12.0; 
  e_(2) =  1.0/9.0;
  e_(3) = -1.0/8.0;

  a_(0,0) = 0.0;     a_(0,1) = 0.0;     a_(0,2) = 0.0;     a_(0,3) = 0.0;
  a_(1,0) = 1.0/2.0; a_(1,1) = 0.0;     a_(1,2) = 0.0;     a_(1,3) = 0.0;
  a_(2,0) = 0.0;     a_(2,1) = 3.0/4.0; a_(2,2) = 0.0;     a_(2,3) = 0.0;
  a_(3,0) = 2.0/9.0; a_(3,1) = 1.0/3.0; a_(3,2) = 4.0/9.0; a_(3,3) = 0.0;
}

OdeIntegratorBase * createRungeKutta23Solver()
{
    return new RungeKutta23();
}
