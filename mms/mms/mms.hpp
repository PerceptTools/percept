#ifndef MMS_HPP
#define MMS_HPP

#include <Sacado.hpp>

namespace mms
{

  typedef Sacado::Fad::DFad<double> FAD_Type;
  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FAD2_Type;

  typedef Sacado::Fad::SFad<double,4> SFAD_Type;
  typedef Sacado::Fad::SFad<Sacado::Fad::SFad<double,4>,4> SFAD2_Type;
  
  typedef SFAD_Type (scalar_SFAD_function)(const SFAD_Type & x,
					 const SFAD_Type & y,
					 const SFAD_Type & z,
					 const SFAD_Type & t
					 );

  typedef SFAD2_Type (scalar_SFAD2_function)(const SFAD2_Type & x,
					   const SFAD2_Type & y,
					   const SFAD2_Type & z,
					   const SFAD2_Type & t
					   );

  enum PARTIAL_DERIVATIVES {_X=0, _Y, _Z, _T};
  
  // first partial derivatives for first order FADs
  double DX(const SFAD_Type &f);
  double DY(const SFAD_Type &f);
  double DZ(const SFAD_Type &f);
  double DT(const SFAD_Type &f);
  
  double   D(const  SFAD_Type &f, const int i);

  SFAD_Type D(const SFAD2_Type &f, const int i);

  // first order partial derivatives for second order FADs
  SFAD_Type DX2(const SFAD2_Type &f);
  SFAD_Type DY2(const SFAD2_Type &f);
  SFAD_Type DZ2(const SFAD2_Type &f);
  SFAD_Type DT2(const SFAD2_Type &f);

  SFAD_Type D2(const SFAD2_Type &f, const int i);

  // NOTE I would like to overload using DX like this:
  //   FAD_Type DX(const FAD2_Type &f) { return f.dx(0);}
  // but that code does not compile

/*   FAD_Type DX2(const FAD2_Type &f) { return f.dx(0);} */
/*   FAD_Type DY2(const FAD2_Type &f) { return f.dx(1);} */
/*   FAD_Type DZ2(const FAD2_Type &f) { return f.dx(2);} */
/*   FAD_Type DT2(const FAD2_Type &f) { return f.dx(3);} */
  
/*   // second order partial derivatives for second order FADs */
/*   double DXX(const FAD2_Type &f) { return DX(DX2(f));} */
/*   double DXY(const FAD2_Type &f) { return DX(DY2(f));} */
/*   double DXZ(const FAD2_Type &f) { return DX(DZ2(f));} */

/*   double DYX(const FAD2_Type &f) { return DY(DX2(f));} */
/*   double DYY(const FAD2_Type &f) { return DY(DY2(f));} */
/*   double DYZ(const FAD2_Type &f) { return DY(DZ2(f));} */

/*   double DZX(const FAD2_Type &f) { return DZ(DX2(f));} */
/*   double DZY(const FAD2_Type &f) { return DZ(DY2(f));} */
/*   double DZZ(const FAD2_Type &f) { return DZ(DZ2(f));} */
  
/*   // examples of scalar differential operators */
/*   double LAPLACE2D(const FAD2_Type &f) { return DXX(f)+DYY(f);} */
/*   double LAPLACE3D(const FAD2_Type &f) { return DXX(f)+DYY(f)+DZZ(f);} */
  
  // API for scalar source terms  
  class ScalarMMSBase
  {
  public:
    
    ScalarMMSBase() {}
    
    virtual ~ScalarMMSBase() {}
    
    // scalar source term interface
    virtual double source(const FAD_Type & x,
			  const FAD_Type & y,
			  const FAD_Type & z,
			  const FAD_Type & time) = 0;
  };
  
} // namespace mms

#endif // MMS_HPP
