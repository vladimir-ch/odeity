#include "RungeKuttaMerson45.h"

RungeKuttaMerson45::RungeKuttaMerson45()
  : ExplicitRungeKuttaBase( 4, 5, false ) // 4th order, 5 stages, FSAL false
{
  solverName_ = "Runge-Kutta-Merson5(4)5 (opravdu k vyvoji reseni pouzivame reseni 5. radu?)";

  c_(0) = 0.0;
  c_(1) = 1.0/3.0;
  c_(2) = 1.0/3.0;
  c_(3) = 1.0/2.0;
  c_(4) = 1.0;

  b_(0) = 1.0/6.0;
  b_(1) = 0.0; 
  b_(2) = 0.0;
  b_(3) = 2.0/3.0;
  b_(4) = 1.0/6.0;

  e_(0) =  1.0/15.0;
  e_(1) =  0.0; 
  e_(2) = -3.0/10.0;
  e_(3) =  4.0/15.0;
  e_(4) = -1.0/30.0;

  a_(0,0) = 0.0;     a_(0,1) = 0.0;     a_(0,2) = 0.0;   a_(0,3) = 0.0; a_(0,4) = 0.0;
  a_(1,0) = 1.0/3.0; a_(1,1) = 0.0;     a_(1,2) = 0.0;   a_(1,3) = 0.0; a_(1,4) = 0.0;
  a_(2,0) = 1.0/6.0; a_(2,1) = 1.0/6.0; a_(2,2) = 0.0;   a_(2,3) = 0.0; a_(2,4) = 0.0;
  a_(3,0) = 0.125;   a_(3,1) = 0.0;     a_(3,2) = 0.375; a_(3,3) = 0.0; a_(3,4) = 0.0;
  a_(4,0) = 0.5;     a_(4,1) = 0.0;     a_(4,2) = -1.5;  a_(4,3) = 2.0; a_(4,4) = 0.0;
}


OdeIntegratorBase * createRungeKuttaMerson45Solver()
{
    return new RungeKuttaMerson45();
}


