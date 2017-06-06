// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <percept/Percept.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/FunctionOperator.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/Util.hpp>
#include <percept/norm/Norm.hpp>
#include <percept/norm/H1Norm.hpp>
#include <percept/math/Math.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/fixtures/Fixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <Teuchos_ScalarTraits.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

namespace percept {
namespace heavy_tests {

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================

/// This test uses a back door to the function that passes in the element to avoid the lookup of the element when the
///  StringFunction contains references to FieldFunctions
void TEST_norm_string_function_turbo_timings(TurboOption turboOpt, const size_t nxyz = 4, bool only_turbo=false)
{
  EXCEPTWATCH;

  dw().m(LOG_NORM) << "TEST.norm.string_function " << stk::diag::dendl;

  /// create a meta data/bulk data empty pair
  PerceptMesh eMesh(3u);

  if (1)
  {
    // Need a symmetric mesh around the origin for some of the tests below to work correctly (i.e. have analytic solutions)
    const size_t num_x = nxyz;
    const size_t num_y = nxyz;
    const size_t num_z = nxyz;
    std::string config_mesh =
      Ioss::Utils::to_string(num_x) + "x" +
      Ioss::Utils::to_string(num_y) + "x" +
      Ioss::Utils::to_string(num_z) + "|bbox:-0.5,-0.5,-0.5,0.5,0.5,0.5";
	
    eMesh.new_mesh(GMeshSpec(config_mesh));

    eMesh.commit();
  }

  stk::mesh::MetaData& metaData = *eMesh.get_fem_meta_data();
  stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *coords_field = metaData.get_field(stk::topology::NODE_RANK, "coordinates");

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", coords_field, &bulkData,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// the function to be integrated:  sqrt(Integral[x^2, dxdydz]) =?= sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)
  StringFunction sfx("x", Name("sfx"), Dimensions(3), Dimensions(1) );

  ff_coords.add_alias("mc");
  //StringFunction sfcm("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", Name("sfcm"), Dimensions(3), Dimensions(1));
  StringFunction sfx_mc("mc[0]", Name("sfx_mc"), Dimensions(3), Dimensions(1) );

  /// the function to be integrated:  sqrt(Integral[x^2, dxdydz]) =?= sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)

  /// A place to hold the result.
  /// This is a "writable" function (we may want to make this explicit - StringFunctions are not writable; FieldFunctions are
  /// since we interpolate values to them from other functions).
  ConstantFunction sfx_res(0.0, "sfx_res");
  ConstantFunction sfx_res_turbo(0.0, "sfx_res_turbo");
  ConstantFunction sfx_res_slow(0.0, "sfx_res_slow");
  ConstantFunction sfx_res_fast(0.0, "sfx_res_fast");

#define COL_SEP "|"
#define EXPR_CELL_WIDTH (80)

#define TIME_IT2(expr_none,expr_turbo,msg,topt)                         \
  {                                                                     \
    double TURBO_NONE_time    = 0;                                      \
    double TURBO_ON_time = 0;                                           \
    TIME_IT(expr_none,TURBO_NONE_time);                                 \
    TIME_IT(expr_turbo,TURBO_ON_time);                                  \
    if (1) std::cout << msg << #topt << " for expression= " << QUOTE(expr_none) << " timings= " << std::endl; \
    if (1) std::cout << "TURBO_NONE_time= " << TURBO_NONE_time << " "   \
                     << ( turboOpt==TURBO_ELEMENT?"TURBO_ELEMENT_time":"TURBO_BUCKET_time") <<"= " << TURBO_ON_time \
                     << " ratio= " << TURBO_NONE_time/TURBO_ON_time << std::endl; \
  }

  int numIter = 1;
  for (int iter = 0; iter < numIter; iter++)
  {
    /// Create the operator that will do the work
    /// get the l2 norm
    Norm<2> l2Norm      (bulkData, &metaData.universal_part(), TURBO_NONE);
    Norm<2> l2Norm_turbo(bulkData, &metaData.universal_part(), turboOpt);

    TIME_IT2((only_turbo ? (void)0 : l2Norm(sfx, sfx_res)); , l2Norm_turbo(sfx, sfx_res_turbo);, "Should be the same turboOpt= ", turboOpt );
    TIME_IT2((only_turbo ? (void)0 : l2Norm(sfx_mc, sfx_res_slow)); , l2Norm_turbo(sfx, sfx_res_fast); , "StringFunction with ref to FF, slow vs. fast" ,turboOpt );

    double sfx_expect = std::sqrt(0.25/3.);
    if (!only_turbo) EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    EXPECT_DOUBLE_EQ_APPROX_TOL( sfx_expect, sfx_res_turbo.getValue(), 1.e-8);
    if (!only_turbo) EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_slow.getValue());
    EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_fast.getValue());

    /// the function to be integrated:  (Integral[ abs(x), dxdydz]) =?= (2 * |x|^2/2 @ [0, 0.5]) ==> .25)
    Norm<1> l1Norm(bulkData, &metaData.universal_part(), TURBO_NONE);
    l1Norm(sfx, sfx_res);

    Norm<1> l1Norm_turbo(bulkData, &metaData.universal_part(), TURBO_NONE);
    l1Norm_turbo(sfx, sfx_res_turbo);

    sfx_expect = 0.25;
    if (!only_turbo) EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());

    /// the function to be integrated:  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest1.py)
    StringFunction sfxyz("x*y*z", Name("sfxyz"), Dimensions(3), Dimensions(1) );

    TIME_IT2((only_turbo ? (void)0 : l2Norm(sfxyz, sfx_res)); ,  l2Norm_turbo(sfxyz, sfx_res_turbo); ,  "should be the same", turboOpt );

    sfx_expect = 0.0240562612162344;
    if (!only_turbo) EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
    EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());

    /// the function to be integrated (but over a rotated domain):  sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest2.py)
    /// now rotate the mesh
    Math::Matrix rmz = Math::rotationMatrix(2, 30);
    Math::Matrix rm = rmz;
    eMesh.transform_mesh(rm);

    TIME_IT2((only_turbo ? (void)0 : l2Norm(sfxyz, sfx_res)); , l2Norm_turbo(sfxyz, sfx_res_turbo); , "should be the same", turboOpt );

    sfx_expect = 0.0178406008037016;
    // NOTE: we need extra quadrature accuracy to reproduce this result (cubDegree==4 in IntegratedOp almost gets it right)
    //   for now, we are satisfied with 3 digits
    //EXPECT_DOUBLE_EQ_APPROX(sfx_res.getValue(), sfx_expect);
    if (!only_turbo && std::fabs(sfx_res.getValue()-sfx_expect) > 0.01*sfx_expect)
    {
      EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res.getValue());
      EXPECT_TRUE(false);
    }
    if (!only_turbo && std::fabs(sfx_res_turbo.getValue()-sfx_expect) > 0.01*sfx_expect)
    {
      EXPECT_DOUBLE_EQ_APPROX( sfx_expect, sfx_res_turbo.getValue());
      EXPECT_TRUE(false);
    }
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

TEST(heavy_norm, string_function_turbo_timings)
{
  EXCEPTWATCH;
  TEST_norm_string_function_turbo_timings(TURBO_ELEMENT);
}
TEST(heavy_norm, string_function_turbo_timings_bucket)
{
  EXCEPTWATCH;
  TEST_norm_string_function_turbo_timings(TURBO_BUCKET);
}
TEST(heavy_norm, string_function_turbo_timings_bucket_prof)
{
  EXCEPTWATCH;
  //TEST_norm_string_function_turbo_timings(TURBO_BUCKET, 50, true);
  TEST_norm_string_function_turbo_timings(TURBO_BUCKET, 10, true);
}


}
}

