#include <mms/mms.hpp>

#include <assert.h>

namespace mms
{
  // first partial derivatives for first order SFADs
  double DX(const SFAD_Type &f) { return f.dx(_X);}
  double DY(const SFAD_Type &f) { return f.dx(_Y);}
  double DZ(const SFAD_Type &f) { return f.dx(_Z);}
  double DT(const SFAD_Type &f) { return f.dx(_T);}

  double D(const SFAD_Type &f, const int i) { 
    assert(_X <= i && i <= _T);
    return f.dx(i);
  }

  SFAD_Type D(const SFAD2_Type &f, const int i) { 
    assert(_X <= i && i <= _T);
    return f.dx(i);
  }

  // first order partial derivatives for second order SFADs
  SFAD_Type DX2(const SFAD2_Type &f) {return f.dx(_X);}
  SFAD_Type DY2(const SFAD2_Type &f) {return f.dx(_Y);}
  SFAD_Type DZ2(const SFAD2_Type &f) {return f.dx(_Z);}
  SFAD_Type DT2(const SFAD2_Type &f) {return f.dx(_T);}

  SFAD_Type D2(const SFAD2_Type &f, const int i) {
    assert(_X <= i && i <= _T);
    return f.dx(i);
  }
} // namespace mms
