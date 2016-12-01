// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/util/Loops.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>

#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>
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
namespace unit_tests {

#define EXTRA_PRINT 0

      static int print_infoLevel = 0;

//=============================================================================
//=============================================================================
//=============================================================================

TEST(function, fieldFunction_demo_1_0_0)
{
  EXCEPTWATCH;

  // start_demo_fieldFunction_1
  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
  eMesh.commit();
  eMesh.print_info("fieldFunction_demo_1_0_0",  print_infoLevel);

  // the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  // create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3);

  // here we could evaluate this field function
  double x=0.123, y=0.234, z=0.345, time=0.0;
  eval_vec3_print(x, y, z, time, ff_coords);
  // end_demo

}

TEST(function, fieldFunction_read_print)
{
  EXCEPTWATCH;
  // just reads a mesh file and prints some info about the meta data

  bool print_info = false;

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";

  PerceptMesh eMesh(3u);
  eMesh.new_mesh_read_only(GMeshSpec(config_mesh));

  stk::mesh::MetaData& metaData = *eMesh.get_fem_meta_data();

  const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();

  unsigned nparts = parts.size();
  if (print_info) std::cout << "Number of parts = " << nparts << std::endl;

  // here's where we can add parts
  // ...
  // ... then we would have to commit the metaData

  const stk::mesh::FieldVector & fields =  metaData.get_fields();
  unsigned nfields = fields.size();
  if (print_info)
  {
    std::cout << "Number of fields = " << fields.size() << std::endl;
    for (unsigned ifld = 0; ifld < nfields; ifld++)
    {
      stk::mesh::FieldBase *field = fields[ifld];
      if (print_info) std::cout << "Field[" << ifld << "]= " << field->name() << " rank= " << field->field_array_rank() << std::endl;
      if (print_info) std::cout << *field << std::endl;
      unsigned nfr = field->restrictions().size();
      if (print_info) std::cout << " number of field restrictions= " << nfr << std::endl;
      for (unsigned ifr = 0; ifr < nfr; ifr++)
      {
        const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
	stk::mesh::Selector frselector = fr.selector();
        if (print_info) std::cout << " field restriction " << ifr << " stride[0] = " << fr.num_scalars_per_entity() <<
                         " type= " << field->entity_rank() << " selector= " << frselector << std::endl;
      }
    }
  }
}


//=============================================================================
//=============================================================================
//=============================================================================

#define EXPR_COORD_MAG (sqrt(x*x + y*y + z*z))

class CheckCoordMag : public GenericFunction
{
public:
  bool m_error;
  std::string m_name;
  CheckCoordMag(std::string name="") : m_error(false), m_name(name) {}
  virtual void operator()(MDArray& domain, MDArray& codomain, double time_value_optional=0.0)
  {
    double x = domain(0);
    double y = domain(1);
    double z = domain(2);
    double v = EXPR_COORD_MAG;
    double cmag_field_node = codomain(0);
    //EXPECT_DOUBLE_EQ(v, cmag_field_node);
    if (fabs(v-cmag_field_node) > 1.e-6)
    {
      std::cout << "CheckCoordMag:: " << m_name <<
        " v= " << v << " x= " << x << " y= " << y << " z= "<< z << " cmag_field_node= " << cmag_field_node << std::endl;
      Util::pause(true, "cmag_field_node");
      ASSERT_NEAR(v, cmag_field_node, 1.e-9);
      m_error = true;
    }
  }

};

//=============================================================================
//=============================================================================
//=============================================================================
TEST(function, fieldFunction_demo_1)
{
  EXCEPTWATCH;


  {
    stk::io::util::Gmesh_STKmesh_Fixture gms(MPI_COMM_WORLD, "3x3x3|bbox:0,0,0,1,1,1");
    std::cout << "gms= " << &gms << std::endl;
  }

  std::cout << "gms= end"  << std::endl;

  // start_demo_fieldFunction_1
  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
  eMesh.commit();

  // the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  // create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3);

  // here we could evaluate this field function
  double x=0.123, y=0.234, z=0.345, time=0.0;
  eval_vec3_print(x, y, z, time, ff_coords);
  // end_demo

}

TEST(function, fieldFunction_demo_2)
{
  EXCEPTWATCH;

  // start_demo_fieldFunction_2
  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1")); // create a 3x3x3 hex mesh in the unit cube

  // add a new field
  // NOTE: we have to create the fields here before committing the mesh
  int vectorDimension = 0;  // signifies a scalar field
  eMesh.add_field("coords_mag_field", stk::topology::NODE_RANK, vectorDimension);
  eMesh.commit();

  // the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  // get the new field created by PerceptMesh
  stk::mesh::FieldBase* coords_mag_field = eMesh.get_field(stk::topology::NODE_RANK, "coords_mag_field");

  // create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh,  3, 3);
  eval_vec3_print(0.1,0.1,0.1,0.0, ff_coords);

  // create a StringFunction to define the magnitude of the coordinates
  StringFunction coords_mag_sf( "sqrt(x*x + y*y + z*z)" , Name("coords_mag_sf"), 3, 1);
  double x=0.123, y=0.234, z=0.345;
  double vv = std::sqrt(x*x + y*y + z*z);            // evaluate the expression in C++
  double v1 = eval(x, y, z, 0, coords_mag_sf);       // evaluate the analytic expression
  ASSERT_NEAR(vv, v1, 1.e-9);                          // the two results should be the same

  // Interpolate the analytic field defined by "coords_mag_sf" to the field we created to hold the coordinate magnitude field
  // 1. create a field function to represent the new coordinate magnitude field, and interpolate the string function to its nodes
  FieldFunction coords_mag_field_function("coords_mag_field_function", coords_mag_field, eMesh, 3, 1);

  coords_mag_field_function.interpolateFrom(coords_mag_sf);

  // We can now write the model with the new coordinates magnitude field to an Exodus file
  eMesh.save_as("./cube_hex8_withCoordMag_out.e");
  // end_demo

  // start_demo_fieldFunction_3

  // tell Percept that we want to refer to the ff_coords FieldFunction by a simple alias "mc"
  ff_coords.add_alias("mc");

  // define a new StringFunction that does the same thing as coords_mag_sf, evaluates the coordinate magnitudes
  StringFunction sfcm("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", Name("sfcm"), 3, 1);
  // end_demo

}

TEST(function, fieldFunction_readMesh_createField_interpolateFrom)
{
  EXCEPTWATCH;
  // more i/o; adding a field; writing the resulting mesh; creating a FieldFunction, invoking interpolateFrom

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";

  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec(config_mesh));
  int vectorDimension = 0;  // signifies a scalar field
  eMesh.add_field("coords_mag_field", stk::topology::NODE_RANK, vectorDimension);
  eMesh.commit();

  unsigned p_rank = eMesh.get_bulk_data()->parallel_rank();
  //unsigned p_size = eMesh.get_bulk_data()->parallel_size();
  Util::setRank(p_rank);

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  /// get the new field created by readModelCreateOptionalFields()
  stk::mesh::FieldBase* coords_mag_field = eMesh.get_field(stk::topology::NODE_RANK, "coords_mag_field");
  VERIFY_OP_ON(coords_mag_field, !=, 0, "TEST::function::fieldFunction_readMesh_createField_interpolateFrom: null coords_mag_field");

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3, FieldFunction::SIMPLE_SEARCH );

  /// here we could evaluate this field function
  if (0)
  {
    if (p_rank == 0)
    {
      std::cout << "TEST::function::fieldFunction_readMesh_createField_interpolateFrom eval ff_coords=" << std::endl;
      eval_vec3_print(0.1, 0.2, 0.3, 0.0, ff_coords);
    }
  }

  StringFunction coords_mag_sf( EXPAND_AND_QUOTE(EXPR_COORD_MAG) , Name("coords_mag_sf"), 3, 1);

  /// create a field function to represent the new coordinate magnitude field, and interpolate the string function to its nodes
  FieldFunction coords_mag_field_function("coords_mag_field_function", coords_mag_field, eMesh, 3, 3, FieldFunction::SIMPLE_SEARCH );
  coords_mag_field_function.interpolateFrom(coords_mag_sf);

  /// check that the coordinates mag field is set correctly
  {
    EXCEPTWATCH;
    CheckCoordMag checkCoordMag;
    //if (!p_rank) std::cout << "checkCoordMag..." << std::endl;
    nodalOpLoop(*eMesh.get_bulk_data(), checkCoordMag, coords_mag_field);
    //if (!p_rank) std::cout << "checkCoordMag...done" << std::endl;
    EXPECT_FALSE(checkCoordMag.m_error);
  }

  try {
    ff_coords.add_alias("mc");
    StringFunction sfcm("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", Name("sfcm"), Dimensions(3), Dimensions(1));

    double tol1 = 1.e-12;

    {
      MDArray vv = eval_vec3(0.1, 0.2, 0.3, 0.0, ff_coords);

      ASSERT_NEAR(vv(0), 0.1, tol1);
      ASSERT_NEAR(vv(1), 0.2, tol1);
      ASSERT_NEAR(vv(2), 0.3, tol1);
    }

    {
      double vv = eval(0.1, 0.2, 0.3, 0.0, sfcm);
      double v_expect = std::sqrt(0.1*0.1 + 0.2*0.2 + 0.3*0.3);
      ASSERT_NEAR(vv, v_expect, tol1);
    }

    coords_mag_field_function.interpolateFrom(sfcm);
    CheckCoordMag checkCoordMag1(std::string(EXPAND_AND_QUOTE(__FILE__))+": "+toString(__LINE__));
    if (!p_rank) std::cout << "checkCoordMag1..." << std::endl;
    nodalOpLoop(*eMesh.get_bulk_data(), checkCoordMag1, coords_mag_field);
    if (!p_rank) std::cout << "checkCoordMag1...done" << std::endl;
    EXPECT_FALSE(checkCoordMag1.m_error);
  }
  catch ( const std::exception * X ) {
    std::cout << "  unexpected exception: " << X->what() << std::endl;
    exit(123);
  }
  catch ( const std::exception & X ) {
    std::cout << "  unexpected exception: " << X.what() << std::endl;
    exit(124);
  }
  catch( ... ) {
    std::cout << "  ... exception" << std::endl;
    exit(125);
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

enum {NPTS = 4};
static double testpoints[NPTS][4] = {
  {0.1234,     0.5678,    0.9,    0.812   },
  {0.1234e-3,  0.5678e-5, 0.97,   0.01    },
  {0.101,      0.02,      0.1020, 0.0122  },
  {0.003,      0.89,      0.01,   0.5     }
};

TEST(function, fieldFunction_multiplePoints)
{
  EXCEPTWATCH;
  std::cout << "TEST::function::fieldFunction_multiplePoints" <<  std::endl;
  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";

  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec(config_mesh));
  int vectorDimension = 0;  // signifies a scalar field
  eMesh.add_field("coords_mag_field", stk::topology::NODE_RANK, vectorDimension);
  eMesh.commit();

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  FieldFunction ff_coords("ff_coords", f_coords, eMesh,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );
  MDArray val1 = eval_vec3(0.2, 0.3, 0.4, 0.0, ff_coords);
  std::cout << "val1= \n" << val1 << std::endl;

  MDArray points(NPTS, 3);
  MDArray output(NPTS, 3);
  MDArray output_expect(NPTS, 3);

  //StringFunction sf1("x+y*z");
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    double x = testpoints[ipts][0];
    double y = testpoints[ipts][1];
    double z = testpoints[ipts][2];
    double t = testpoints[ipts][3];
    points(ipts, 0) = x;
    points(ipts, 1) = y;
    points(ipts, 2) = z;
    //points(ipts, 3) = t;

    //std::cout << "field_op: ipts= " << ipts << std::endl;

    MDArray vec = eval_vec3(x, y, z, t, ff_coords);
    EXPECT_NEAR(vec(0), x, fabs(1.e-5*x));
    EXPECT_NEAR(vec(1), y, fabs(1.e-5*y));
    EXPECT_NEAR(vec(2), z, fabs(1.e-5*z));
    output_expect(ipts, 0) = x;
    output_expect(ipts, 1) = y;
    output_expect(ipts, 2) = z;
  }
  std::cout << "field_op: NPTS= " << NPTS << std::endl;
  //         ff_coords.setDomainDimensions(Dimensions(NPTS,3));
  //         ff_coords.setCodomainDimensions(Dimensions(NPTS,3));
  ff_coords.setDomainDimensions(Dimensions(3));
  ff_coords.setCodomainDimensions(Dimensions(3));
  ff_coords(points, output, 0.0);
  for (unsigned ipts = 0; ipts < NPTS; ipts++)
  {
    EXPECT_NEAR(output(ipts, 0), output_expect(ipts, 0), 1.e-5*(fabs(output_expect(ipts,0))) );
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

TEST(function, fieldFunction_point_eval_verify)
{
  EXCEPTWATCH;
  /// test evaluation of field function at a point

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";

  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec(config_mesh));

  eMesh.commit();
  // no need for this in create mode: eMesh.readBulkData();

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// here we evaluate this field function
  MDArray val1 = eval_vec3(0.2, 0.3, 0.4, 0.0, ff_coords);
  //std::cout << "eval = \n" << val1 << std::endl;

  stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();


  bool didCatch = false;
  try {
    // evaluate a point that is known to be outside the domain
    MDArray val10 = eval_vec3(1.2, 1.3, 1.4, 0.0, ff_coords);
  }
  catch ( const std::exception & X ) {
    std::cout << "  expected to catch this exception: " << X.what() << std::endl;
    didCatch = true;
  }
  catch( ... ) {
    std::cout << "  P:" << bulkData.parallel_rank()
              << " Caught unknown exception"
              << std::endl ;
    std::cout.flush();
    didCatch = false;
  }
  EXPECT_TRUE(didCatch);

  //double value = eval(1.2, 2.3, 3.4, 0.0, ff_coords);
  MDArray pts(3);
  MDArray output_pts(3);
  pts(0) = 0.2; pts(1) = 0.3; pts(2) = 0.4; //pts(3) = 0.0;
  ff_coords(pts, output_pts);
  ASSERT_NEAR(pts(0), output_pts(0), 1.e-9);
  ASSERT_NEAR(pts(1), output_pts(1), 1.e-9);
  ASSERT_NEAR(pts(2), output_pts(2), 1.e-9);
}

//=============================================================================
//=============================================================================
//=============================================================================

TEST(function, fieldFunction_point_eval_deriv_verify)
{
  EXCEPTWATCH;
  /// test evaluation of field function at a point

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";

  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec(config_mesh));

  int vectorDimension = 0;  // signifies a scalar field
  stk::mesh::FieldBase *f_test = eMesh.add_field("test", stk::topology::NODE_RANK, vectorDimension);
  eMesh.commit();
  // no need for this in create mode: eMesh.readBulkData();

  StringFunction sf1("x + y + z + t");

  FieldFunction ff1("ff1", f_test, eMesh, 3, 1);
  ff1.interpolateFrom(sf1);

  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  /// create a field function from the existing coordinates field
  FieldFunction ff_coords("ff_coords", f_coords, eMesh,
                          Dimensions(3), Dimensions(3), FieldFunction::SIMPLE_SEARCH );

  /// here we evaluate this field function
  MDArray val1 = eval_vec3(0.2, 0.3, 0.4, 0.0, ff_coords);
  //std::cout << "eval = \n" << val1 << std::endl;

  //double value = eval(1.2, 2.3, 3.4, 0.0, ff_coords);
  MDArray pts(3);
  MDArray output_pts(3);
  pts(0) = 0.2; pts(1) = 0.3; pts(2) = 0.4; //pts(3) = 0.0;
  ff_coords(pts, output_pts);
  ASSERT_NEAR(pts(0), output_pts(0), 1.e-9);
  ASSERT_NEAR(pts(1), output_pts(1), 1.e-9);
  ASSERT_NEAR(pts(2), output_pts(2), 1.e-9);

#if 0
  Teuchos::RCP<Function> deriv_ff = ff1.gradient();
  MDArray outp1(3);
  deriv_ff->operator()(pts, outp1);
  ASSERT_NEAR(outp1(0), 1.0, 1.e-9);
  ASSERT_NEAR(outp1(1), 1.0, 1.e-9);
  ASSERT_NEAR(outp1(2), 1.0, 1.e-9);
#endif
}

//=============================================================================
//=============================================================================
//=============================================================================

TEST(function, fieldFunction_point_eval_timing)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  /// test evaluation of field function at a point

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";

  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec(config_mesh));

  eMesh.commit();
  // no need for this in create mode: eMesh.readBulkData();

  //unsigned p_rank = eMesh.get_bulk_data()->parallel_rank();
  unsigned p_size = eMesh.get_bulk_data()->parallel_size();
  // FIXME
  if (p_size > 1) return;
  /// the coordinates field is always created by the PerceptMesh read operation, here we just get the field
  stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

  for (unsigned iSearchType = 0; iSearchType < 2; iSearchType++)
  {
    /// create a field function from the existing coordinates field
    FieldFunction::SearchType search_type = (iSearchType == 0 ? FieldFunction::SIMPLE_SEARCH : FieldFunction::STK_SEARCH);
    //std::cout <<  "P[" << Util::get_rank() <<  "] search_type = " << search_type << " = " << typeid(search_type).name() << std::endl;
    FieldFunction ff_coords("ff_coords", f_coords, eMesh,
                            Dimensions(3), Dimensions(3), search_type
                            );

    // The first point that is evaluated fires the setup of the stk::search data structure (oct-tree, bih-tree)
    double t1st =  stk::wall_time();
    MDArray val11 = eval_vec3(0.2, 0.3, 0.4, 0.0, ff_coords);
    val11 = eval_vec3(0.2, 0.3, 0.4, 0.0, ff_coords);
    t1st = stk::wall_dtime(t1st);

    // timings
    unsigned numIter = 10000;
    //unsigned numIter = 1000;
    MDArray pts(3);
    MDArray output_pts(3);

    // ensure the same set of random data each run
    Teuchos::ScalarTraits<double>::seedrandom(12345);

    double tstart =  stk::wall_time();
    for (unsigned iter = 0; iter < numIter; iter++)
    {
      double rnd = Teuchos::ScalarTraits<double>::random();
      pts(0) = (rnd+1.0)/2.0;
      rnd = Teuchos::ScalarTraits<double>::random();
      pts(1) = (rnd+1.0)/2.0;
      rnd = Teuchos::ScalarTraits<double>::random();
      pts(2) = (rnd+1.0)/2.0;

      // FIXME
      //!! pts(0) = 0.2; pts(1) = 0.3; pts(2)= 0.4;
      // FIXME
      ff_coords(pts, output_pts, 0.0);
#if 0
      EXPECT_DOUBLE_EQ(pts(0), output_pts(0));
      EXPECT_DOUBLE_EQ(pts(1), output_pts(1));
      EXPECT_DOUBLE_EQ(pts(2), output_pts(2));
#endif
    }

    double total_time = stk::wall_dtime(tstart);
    if (1 | EXTRA_PRINT) std::cout
                            << "TEST::function::fieldFunction_point_eval_timing: "
                            << " for search_type= " << (iSearchType==0?"SIMPLE_SEARCH":"STK_SEARCH")<< "\n"
                            << "    time for 1st eval=  " << t1st << "\n"
                            << "    for " << numIter << " iterations, evaluating field(x,y,z) time = " << total_time  << "\n"
                            << "    average per point lookup and eval time = " << (total_time/((double)numIter)) << std::endl;
  }
  //std::cout << "P[" << Util::get_rank() <<  "] TEST::function::fieldFunction_point_eval_timing done " << std::endl;
}

static void test_sync(PerceptMesh& eMesh, bool sync_shared, bool sync_aura)
{
  stk::mesh::FieldBase* pressure_field =  eMesh.get_field(stk::topology::NODE_RANK, "pressure");
  unsigned p_rank = eMesh.get_bulk_data()->parallel_rank();
  unsigned p_size = eMesh.get_bulk_data()->parallel_size();
  (void)p_size;

  const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( stk::topology::NODE_RANK );

  enum Type{ Owned, Shared, Ghost };
  std::string types[] = {"Owned", "Shared", "Ghost" };
  if (!p_rank) std::cout << "\n\n=========================\nsrk_sync: sync_shared= " << sync_shared << " sync_aura= " << sync_aura << std::endl;

  std::ostringstream out;
  out << "\n\n=========================\ntest_sync: sync_shared= " << sync_shared << " sync_aura= " << sync_aura << "\n";
  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
    {
      {
        stk::mesh::Bucket & bucket = **k ;

        const unsigned num_elements_in_bucket = bucket.size();

        for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
          {
            stk::mesh::Entity entity = bucket[iEntity];
            double * const p = eMesh.field_data( pressure_field , entity );
            stk::mesh::EntityId id=eMesh.identifier(entity);

            int type=Owned;
            if (bucket.owned())
              {
                //p[0] = (p_rank?200+id:100+id);
                p[0] = (p_rank+1)*100+id;
              }
            else if (bucket.shared())
              {
                //p[0] = -(p_rank?double(100+id):double(200+id));
                p[0] = -double((eMesh.parallel_owner_rank(entity)+1)*100+id);
                type=Shared;
              }
            else
              {
                //p[0] = (p_rank?2000+id:1000+id);
                p[0] = ((p_rank+1)*1000 + id);
                type=Ghost;
                //std::cout << "P["<<p_rank<<"] ghost= " << p[0] << std::endl;
              }
            out << "P["<<p_rank<<"] id= " << eMesh.identifier(entity) << " p= " << p[0] << " type= " << types[type] << std::endl;

          }
      }
    }
  std::cout << out.str() << std::endl;
  eMesh.save_as("./field_before_sync.e");

  {
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(pressure_field);

    // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
    if (sync_aura) stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);

    // the shared part (just the shared boundary)
    if (sync_shared) stk::mesh::communicate_field_data(*eMesh.get_bulk_data()->ghostings()[0], fields);
  }
  std::ostringstream out1;
  out1 << "test_sync: sync_shared= " << sync_shared << " sync_aura= " << sync_aura << "\n";

  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
    {
      {
        stk::mesh::Bucket & bucket = **k ;

        const unsigned num_elements_in_bucket = bucket.size();

        for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
          {
            stk::mesh::Entity entity = bucket[iEntity];
            stk::mesh::EntityId id=eMesh.identifier(entity);
            double * const p = eMesh.field_data( pressure_field , entity );
            int type=Owned;
            //double p_e = (p_rank?200+id:100+id);
            double p_e = (p_rank+1)*100+id;
            if (bucket.owned())
              {
                ASSERT_NEAR(p[0], p_e, 1.e-6);
              }
            else if (bucket.shared())
              {
                //p_e = (p_rank?100+id:200+id);
                p_e = ((eMesh.parallel_owner_rank(entity)+1)*100+id);
                type = Shared;
                if (sync_shared)
                  {
                    if (std::fabs(p[0]-p_e) > 1.e-6)
                      {
                        out1 << "P[" << p_rank << "] ERROR: p[0] = " << p[0] << " p_e= " << p_e << std::endl;
                      }
                    ASSERT_NEAR(p[0], p_e, 1.e-6);
                  }
              }
            else
              {
                //p_e = (p_rank?100+id:200+id);
                p_e = ((eMesh.parallel_owner_rank(entity)+1)*100+id);
                type = Ghost;
                if (sync_aura)
                  ASSERT_NEAR(p[0], p_e, 1.e-6);
              }
            out1 << "P["<<p_rank<<"] after id= " << eMesh.identifier(entity) << " p= " << p[0] << " type= " << types[type] << std::endl;

          }
      }
    }
  std::cout << out1.str() << std::endl;

  eMesh.save_as("./field_sync.e");
}

typedef stk::mesh::Field<int>  PressureFieldType ;

static void test_sync_1(stk::mesh::BulkData& eMesh, PressureFieldType& pressure_field, bool sync_shared, bool sync_aura)
{
  unsigned p_rank = eMesh.parallel_rank();
  unsigned p_size = eMesh.parallel_size();
  (void)p_size;

  const stk::mesh::BucketVector & buckets = eMesh.buckets( stk::topology::NODE_RANK );

  enum Type{ Owned, Shared, Ghost };
  std::string types[] = {"Owned", "Shared", "Ghost" };
  if (!p_rank) std::cout << "\n\n=========================\nsrk_sync: sync_shared= " << sync_shared << " sync_aura= " << sync_aura << std::endl;

  std::ostringstream out;
  out << "\n\n=========================\ntest_sync_1: sync_shared= " << sync_shared << " sync_aura= " << sync_aura << "\n";
  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
    {
      {
        stk::mesh::Bucket & bucket = **k ;

        const unsigned num_elements_in_bucket = bucket.size();

        for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
          {
            stk::mesh::Entity entity = bucket[iEntity];
            int * const p = stk::mesh::field_data<PressureFieldType>( pressure_field , entity );
            stk::mesh::EntityId id=eMesh.identifier(entity);

            int type=Owned;
            if (bucket.owned())
              {
                p[0] = (p_rank+1)*100+id;
              }
            else if (bucket.shared())
              {
                p[0] = -((eMesh.parallel_owner_rank(entity)+1)*100+id);
                type=Shared;
              }
            else
              {
                p[0] = ((p_rank+1)*1000 + id);
                type=Ghost;
                //std::cout << "P["<<p_rank<<"] ghost= " << p[0] << std::endl;
              }
            out << "P["<<p_rank<<"] id= " << eMesh.identifier(entity) << " p= " << p[0] << " type= " << types[type] << std::endl;

          }
      }
    }
  std::cout << out.str() << std::endl;

  {
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(&pressure_field);

    // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
    if (sync_aura) stk::mesh::communicate_field_data(eMesh.aura_ghosting(), fields);

    // the shared part (just the shared boundary)
    if (sync_shared) stk::mesh::communicate_field_data(*eMesh.ghostings()[0], fields);
  }
  std::ostringstream out1;
  out1 << "test_sync_1: sync_shared= " << sync_shared << " sync_aura= " << sync_aura << "\n";

  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
    {
      {
        stk::mesh::Bucket & bucket = **k ;

        const unsigned num_elements_in_bucket = bucket.size();

        for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
          {
            stk::mesh::Entity entity = bucket[iEntity];
            stk::mesh::EntityId id=eMesh.identifier(entity);
            int * const p = stk::mesh::field_data<PressureFieldType>( pressure_field , entity );
            int type=Owned;
            int p_e = (p_rank+1)*100+id;
            if (bucket.owned())
              {
                ASSERT_EQ(p[0], p_e);
              }
            else if (bucket.shared())
              {
                p_e = ((eMesh.parallel_owner_rank(entity)+1)*100+id);
                type = Shared;
                if (sync_shared)
                  {
                    if (p[0] != p_e)
                      {
                        out1 << "P[" << p_rank << "] ERROR: p[0] = " << p[0] << " p_e= " << p_e << std::endl;
                      }
                    ASSERT_EQ(p[0], p_e);
                  }
              }
            else
              {
                p_e = ((eMesh.parallel_owner_rank(entity)+1)*100+id);
                type = Ghost;
                if (sync_aura)
                  ASSERT_EQ(p[0], p_e);
              }
            out1 << "P["<<p_rank<<"] after id= " << eMesh.identifier(entity) << " p= " << p[0] << " type= " << types[type] << std::endl;

          }
      }
    }
  std::cout << out1.str() << std::endl;

}

TEST(function, fieldFunction_field_sync)
{
  EXCEPTWATCH;

  //PerceptMesh eMesh(3u);
  //eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1")); // create a 3x3x3 hex mesh in the unit cube

  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  bool sidesets_on = false;
  percept::QuadFixture<double> fixture( pm , 2 , 2, sidesets_on);
  fixture.set_bounding_box(0,1,0,1);

  percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, false);
  // add a new field
  int vectorDimension = 0;  // signifies a scalar field
  stk::mesh::FieldBase* pressure_field =  eMesh.add_field("pressure", stk::topology::NODE_RANK, vectorDimension);
  PressureFieldType& p_field = fixture.meta_data.declare_field<PressureFieldType>(stk::topology::NODE_RANK, "p");
  stk::mesh::put_field( p_field , fixture.meta_data.universal_part());

  (void)pressure_field;
  (void)p_field;
  fixture.meta_data.commit();
  fixture.generate_mesh();

  unsigned p_rank = eMesh.get_bulk_data()->parallel_rank();
  (void)p_rank;
  unsigned p_size = eMesh.get_bulk_data()->parallel_size();
  if (p_size <= 4)
    {
      bool exercise_field_sync_bug = false;
      test_sync(eMesh, false, false);
      test_sync(eMesh, false, true);
      if (exercise_field_sync_bug || p_size <=2)
        {
          test_sync(eMesh, true, false);
          test_sync(eMesh, true, true);
        }

       test_sync_1(*eMesh.get_bulk_data(), p_field, false, false);
       test_sync_1(*eMesh.get_bulk_data(), p_field, false, true);
       if (exercise_field_sync_bug || p_size <=2)
         {
           test_sync_1(*eMesh.get_bulk_data(), p_field, true, false);
           test_sync_1(*eMesh.get_bulk_data(), p_field, true, true);
         }

    }


}


#if 0
int main()
{
  StringFunction sf1("x + y + z + t");
  sf1(xyz, out);

  StringFunction sf2("x - y");
  sfx("x");
  sfy("y");
  sfxy("x-y");
  sfxy1== sfx-sfy;

  StringFunction sf21dif = sf2-sf1;

  mesh::MetaData m;
  //...  setup field, etc.
  FieldFunction ff1("ff1", part1, field_1);  // can be nodal or elemental

  ff1.interpolateFrom(sf1);

  FieldFunction ff2(ff1);  // copy
  FieldFunction ff3 = ff1;  // copy

  Function ff21diff = ff2-ff1;  // should be zero
  ZeroFunction zero_func;

  //assert( ff21diff == zero_func );
  random_probe_assert(ff21diff, zero_func, 100);

  DifferenceFunction dfsf21(sf1,sf2);
  //assert( sf21dif == dfsf21 );
  random_probe_assert( sf21dif, dfsf21, 100);

  // check copy
  ff2.interpolateFrom(sf2);
  ff2.assertMostlyEqual(ff1);
}
#endif

} // namespace unit_tests
} // namespace percept
