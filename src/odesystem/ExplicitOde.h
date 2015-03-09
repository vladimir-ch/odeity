#ifndef EXPLICITODE_H
#define EXPLICITODE_H

#include "OdeSystemBase.h"

#include <sundials/sundials_types.h>

// forward declaration of Vector
template <typename T> class Vector;

class ExplicitOde : public OdeSystemBase
{
  public:
    ExplicitOde( bool hasJacobian );

    virtual int rhs( const realtype &t, const Vector<realtype> &y, Vector<realtype> &ydot ) = 0;
    bool useJacobian() const;
    void setUseJacobian( bool useJac );
    bool hasJacobian() const;
    virtual int jacobian( const Vector<realtype>& v, Vector<realtype>& Jv,
        const realtype t, const Vector<realtype>& y, const Vector<realtype>& fy,
        const Vector<realtype>& tmp );

  private:
    bool hasJacobian_;
    bool useJacobian_;
};


#endif
