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


#include "ExplicitOde.h"
#include "../utils/Exceptions.h"

ExplicitOde::ExplicitOde( bool hasJacobian )
  :
    hasJacobian_( hasJacobian ),
    useJacobian_( false )
{}



void ExplicitOde::setUseJacobian( bool useJac )
{
    useJacobian_ = useJac && hasJacobian_;
}



bool ExplicitOde::useJacobian() const
{
  return useJacobian_;
}



bool ExplicitOde::hasJacobian() const
{
  return hasJacobian_;
}



int ExplicitOde::jacobian( const Vector<realtype>& v, Vector<realtype>& Jv,
    const realtype t, const Vector<realtype>& y, const Vector<realtype>& fy, const Vector<realtype>& tmp )
{
  Assert( false, ExcPureFunctionCalled() );

  return 0;
}

