/* Copyright (c) 2010 Vladimir Chalupecky
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


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


