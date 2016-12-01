/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <mms/mms.hpp>
#include <stk_util/diag/UserPlugin.hpp>

namespace mms {
namespace regression_tests {

//=============================================================================
//=============================================================================
//=============================================================================

TEST(exactScalarSolution, euler3D)
{
  const double xb=1.234, yb=2.345, zb=2.3, tb = 23.4;

  SFAD_Type x(4, 0, xb), y(4, 1, yb), z(4, 2, zb), t(4, 3, tb);

  const std::string file_name = "./mms_user.so";
  const std::string register_function_name = "register_euler_density";
  
  sierra::Plugin::Registry::rootInstance().
    registerDL(file_name.c_str(), register_function_name.c_str());

  const std::string function_name = "euler_density";
  scalar_SFAD_function * densitySub 
    = sierra::Plugin::UserSubroutine<scalar_SFAD_function>::getFunction(function_name.c_str());

  SFAD_Type rho = (*densitySub) (x,y,z,t);

  // check values
  EXPECT_DOUBLE_EQ(rho.val(), 1.0+(xb-tb)*yb*zb);

  // check derivatives
  EXPECT_DOUBLE_EQ(rho.dx(_X), yb*zb);  
  EXPECT_DOUBLE_EQ(rho.dx(_Y), (xb-tb)*zb);  
  EXPECT_DOUBLE_EQ(rho.dx(_Z), (xb-tb)*yb);  
  EXPECT_DOUBLE_EQ(rho.dx(_T), -yb*zb);  
}
  
}// namespace regression_tests
}// namespace mms

