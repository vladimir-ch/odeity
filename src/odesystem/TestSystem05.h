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

#ifndef TEST_SYSTEM_05_H
#define TEST_SYSTEM_05_H

#include "ExplicitOde.h"

// Stiff DETEST A3
// initial condition y(*) = 1.0
class TestSystem05 : public ExplicitOde {
public:
  TestSystem05( )
  : ExplicitOde( 4 )
  {}

  virtual int rhs( const realtype &t,  const Vector<realtype> &y, Vector<realtype> &ydot )
  {
    ydot(0) = -1.0e4 * y(0) + 1.0e2 * y(1) - 1.0e1 * y(2) + y(3);
    ydot(1) = -1.0e3 * y(1) + 1.0e1 * y(2) - 1.0e1 * y(3);
    ydot(2) = -y(2) + 1.0e1 * y(3);
    ydot(3) = -1.0e-1 * y(3);

    return 0;
  }

};

#endif // TEST_SYSTEM_04_H


