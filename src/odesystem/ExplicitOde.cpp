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


