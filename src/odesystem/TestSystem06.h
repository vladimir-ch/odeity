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

#ifndef TEST_SYSTEM_06_H
#define TEST_SYSTEM_06_H

#include "ExplicitOde.h"

// Matlab ODE Suite - chm6ode
// initial condition y(*) = 0.0
class TestSystem06 : public ExplicitOde {
public:
  TestSystem06( )
  : ExplicitOde( 4 )
  {}

  virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
  {
    realtype K = std::pow(M_E, 20.7 - 1500.0/y(0) );
    ydot(0) = 1.3 * (y(2)-y(0)) + 10400.0 * K * y(1);
    ydot(1) = 1880.0 * (y(3) - y(1)*(1+K));
    ydot(2) = 1752.0 - 269.0 * y(2) + 267.0 * y(0);
    ydot(3) = 0.1 + 320.0 * y(1) - 321.0 * y(3);

    return 0;
  }

};

#endif // TEST_SYSTEM_06_H



