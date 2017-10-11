// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef RegressionTestLocalRefiner_hpp
#define RegressionTestLocalRefiner_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <gtest/gtest.h>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <gtest/gtest.h>
#include <boost/lexical_cast.hpp>
#include <stk_io/IossBridge.hpp>

#include <percept/Percept.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>

#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <percept/RunEnvironment.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

#include <adapt/IEdgeBasedAdapterPredicate.hpp>
#include <adapt/IElementBasedAdapterPredicate.hpp>
#include <adapt/PredicateBasedElementAdapter.hpp>
#include <adapt/PredicateBasedEdgeAdapter.hpp>

#include <percept/function/ElementOp.hpp>

  namespace percept
  {
    namespace regression_tests
    {

      //=============================================================================
      //=============================================================================
      //=============================================================================

      static void normalize(double input_normal[3], double normal[3])
      {
        double sum = std::sqrt(input_normal[0]*input_normal[0]+
                               input_normal[1]*input_normal[1]+
                               input_normal[2]*input_normal[2]);
        normal[0] = input_normal[0] / sum;
        normal[1] = input_normal[1] / sum;
        normal[2] = input_normal[2] / sum;
      }

      static void normalize(double input_output_normal[3])
      {
        normalize(input_output_normal, input_output_normal);
      }

      static double distance(double c0[3], double c1[3])
      {
        return std::sqrt((c0[0]-c1[0])*(c0[0]-c1[0]) + (c0[1]-c1[1])*(c0[1]-c1[1]) + (c0[2]-c1[2])*(c0[2]-c1[2]) );
      }

      static void difference(double v01[3], double c0[3], double c1[3])
      {
        v01[0] = c0[0] - c1[0];
        v01[1] = c0[1] - c1[1];
        v01[2] = c0[2] - c1[2];
      }
      static double dot(double c0[3], double c1[3])
      {
        return c0[0]*c1[0] + c0[1]*c1[1] + c0[2]*c1[2];
      }

      static double plane_dot_product(double plane_point[3], double plane_normal[3], double point[3])
      {
        double normal[3]={0,0,0};
        normalize(plane_normal, normal);
        double dot = 0.0;
        for (int i = 0; i < 3; i++)
          {
            dot += (point[i] - plane_point[i])*normal[i];
          }
        return dot;
      }


      class SetElementRefineFieldValue : public percept::ElementOp
      {
        RefineFieldType_type m_val;
      public:
        SetElementRefineFieldValue(percept::PerceptMesh& eMesh, int val=0) : m_val(val) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          f_data[0] = m_val;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}

      };

      class SetRefineField : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
      public:
        SetRefineField(percept::PerceptMesh& eMesh) : m_eMesh(eMesh) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          double plane_point[3] = {2,0,0};
          double plane_normal[3] = {1, .5, -.5};

          const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK );
          //stk::mesh::Entity const *elem_nodes_e = m_eMesh.get_bulk_data()->begin_entities(element, stk::topology::NODE_RANK );

          unsigned num_node = elem_nodes.size();
          //unsigned num_node = m_eMesh.get_bulk_data()->num_connectivity(element, stk::topology::NODE_RANK );
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

          bool found = false;
          for (unsigned inode=0; inode < num_node-1; inode++)
            {
              stk::mesh::Entity node_i = elem_nodes[ inode ].entity();
              //mesh::Entity node_i = elem_nodes_e[ inode ];
              double *coord_data_i = m_eMesh.field_data(coordField, node_i);

              for (unsigned jnode=inode+1; jnode < num_node; jnode++)
                {
                  stk::mesh::Entity node_j = elem_nodes[ jnode ].entity();
                  double *coord_data_j = m_eMesh.field_data(coordField, node_j);

                  double dot_0 = plane_dot_product(plane_point, plane_normal, coord_data_i);
                  double dot_1 = plane_dot_product(plane_point, plane_normal, coord_data_j);

                  // if edge crosses the plane...
                  if (dot_0*dot_1 < 0)
                    {
                      found=true;
                      break;
                    }
                }
            }
          if (found)
            f_data[0] = 1;
          else
            f_data[0] = 0;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}
      };

      class SetUnrefineField : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
      public:
        SetUnrefineField(percept::PerceptMesh& eMesh) : m_eMesh(eMesh) {}
        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

          bool found = true;
          for (unsigned inode=0; inode < num_node; inode++)
            {
              stk::mesh::Entity node = elem_nodes[ inode ].entity();
              double *coord_data = m_eMesh.field_data(coordField, node);

              if (coord_data[0] < 0.0 || coord_data[1] < 0.0) // || coord_data[2] > 1.1)
                {
                  found=false;
                  break;
                }
            }
          if (found)
            f_data[0] = -1;
          else
            f_data[0] = 0;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}
      };

      struct MyEdgeBasedRefinePredicate : public IEdgeBasedAdapterPredicate {

        MyEdgeBasedRefinePredicate(PerceptMesh& eMesh, stk::mesh::Selector * selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
          IEdgeBasedAdapterPredicate(eMesh, selector, field, tolerance) {}

        /// Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
                       double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
        {
          int mark = 0;
          if (0 == m_selector || (*m_selector)(m_eMesh.bucket(element)))
            {
              double plane_point[3] = {2,0,0};
              double plane_normal[3] = {1, 0, 0};

              double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
              double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

              // if edge crosses the plane...
              if (dot_0 * dot_1 < 0)
                {
                  mark |= DO_REFINE;
                }

              if (coord0[1] < 0 && coord1[1] < 0)
                {
                  mark |= DO_UNREFINE;
                }
            }

          return mark;
        }
      };

      struct PlaneShock
      {
        double plane_point_init[3]; // = {2 + shock_displacement,0,0};
        double plane_normal[3]; // = {1, 0, 0};
        double plane_point[3]; // = {2 + shock_displacement,0,0};
        double shock_width;
        
        PlaneShock()
        {
          plane_point_init[0]=0;
          plane_point_init[1]=0;
          plane_point_init[2]=0;
          plane_point[0]=0;
          plane_point[1]=0;
          plane_point[2]=0;
          plane_normal[0]=1;
          plane_normal[1]=0;
          plane_normal[2]=0;
          shock_width = 0.0;
        }

        void setCurrentPlanePoint(double shock_displacement)
        {
          normalize(plane_normal);
          plane_point[0] = plane_point_init[0] + shock_displacement*plane_normal[0];
          plane_point[1] = plane_point_init[1] + shock_displacement*plane_normal[1];
          plane_point[2] = plane_point_init[2] + shock_displacement*plane_normal[2];
        }

      };

      static double shock_width = 1./5.0;

      static double shock_function(double x)
      {
        // normalize by width
        double width = shock_width; // 1./5.0;
        x /= width;
        return std::tanh(x);
      }


      static double shock_diff(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh,
                               stk::mesh::Entity node0, stk::mesh::Entity node1, double *coord0, double *coord1, PlaneShock& shock, double shock_displacement)
      {
        shock.setCurrentPlanePoint(shock_displacement);
        double *plane_point = shock.plane_point;
        double *plane_normal = shock.plane_normal;

        double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
        double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

        double v01[3] = {0,0,0};
        difference(v01, coord1, coord0);
        normalize(v01);
        normalize(plane_normal);
        double v01dotn = std::abs(dot(v01, plane_normal));

        double d01p = std::abs(dot_0)+std::abs(dot_1);
        dot_0 = shock_function(dot_0);
        dot_1 = shock_function(dot_1);

        if (nodal_refine_field)
          {
            double *fd0 = eMesh.field_data(nodal_refine_field, node0);
            double *fd1 = eMesh.field_data(nodal_refine_field, node1);
            fd0[0] = dot_0;
            fd1[0] = dot_1;
          }

        double d01 = distance(coord0, coord1);

        return (1 + 0*d01 + 0*d01p + 0*v01dotn)*std::abs(dot_0 - dot_1);
      }

      static int shock_diff1(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh,
                               stk::mesh::Entity node0, stk::mesh::Entity node1, double *coord0, double *coord1, PlaneShock& shock, double shock_displacement)
      {
        shock.setCurrentPlanePoint(shock_displacement);
        double *plane_point = shock.plane_point;
        double *plane_normal = shock.plane_normal;

        if (distance(coord0, coord1) < 1./200.)
          return DO_UNREFINE;

        double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
        double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

        if (nodal_refine_field)
          {
            double *fd0 = eMesh.field_data(nodal_refine_field, node0);
            double *fd1 = eMesh.field_data(nodal_refine_field, node1);
            fd0[0] = dot_0 < 0? -1 : 1;
            fd1[0] = dot_1 < 0? -1 : 1;
          }

        if (dot_0*dot_1 < 1.e-6)
          return DO_REFINE;
        else
          return DO_UNREFINE;

      }

      // This can be used as an edge or element-based predicate

      struct ShockBasedRefinePredicate : public IEdgeBasedAdapterPredicate, IElementBasedAdapterPredicate {

        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase * m_nodal_refine_field;
        PlaneShock m_shock;
        double m_shock_displacement;
        double m_shock_diff_criterion;

        ShockBasedRefinePredicate(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh, stk::mesh::Selector* selector, stk::mesh::FieldBase *field, double tolerance,
                                  PlaneShock shock, double shock_displacement=0, double shock_diff_criterion=0.4) :
          IEdgeBasedAdapterPredicate(eMesh, selector, field, tolerance),
          IElementBasedAdapterPredicate(eMesh, selector, field, tolerance),
          m_eMesh(eMesh),m_nodal_refine_field(nodal_refine_field),
          m_shock(shock), m_shock_displacement(shock_displacement),
          m_shock_diff_criterion(shock_diff_criterion) {}



        /// Element-based: Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity element)
        {
          if (0 == m_eb_selector || (*m_eb_selector)(m_eMesh.bucket(element)))
            {

              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
              int spatialDimension = m_eMesh.get_spatial_dim();
              CellTopology cell_topo(cell_topo_data);
              const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);

              CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

              unsigned numSubDimNeededEntities = 0;
              numSubDimNeededEntities = cell_topo_data->edge_count;

              unsigned ref_count=0, unref_count=0;
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = m_eMesh.field_data( *coordField , node0 );
                  double * const coord1 = m_eMesh.field_data( *coordField , node1 );
                  double  dcoord0[3] = {coord0[0],coord0[1], (spatialDimension==2?0:coord0[2])};
                  double  dcoord1[3] = {coord1[0],coord1[1], (spatialDimension==2?0:coord1[2])};

                  int markInfo = this->operator()(element, iSubDimOrd, node0, node1, dcoord0, dcoord1, 0);
                  bool do_ref = markInfo & DO_REFINE;
                  if (do_ref)
                    {
                      ++ref_count;
                    }
                  bool do_unref = markInfo & DO_UNREFINE;
                  if (!do_ref && do_unref)
                    {
                      ++unref_count;
                    }
                }
              if (ref_count == numSubDimNeededEntities)
                return DO_REFINE;
              if (unref_count == numSubDimNeededEntities)
                return DO_UNREFINE;
            }

          return DO_NOTHING;
        }

        /// Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
                        double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
        {
          int mark=0;
          if (0 == m_selector || (*m_selector)(m_eMesh.bucket(element)))
            {

              // refine check
              double d01 = shock_diff(m_nodal_refine_field, m_eMesh, node0, node1, coord0, coord1, m_shock, m_shock_displacement);
              if ( d01 > m_shock_diff_criterion)  // 0.05, 0.2
                {
                  mark |= DO_REFINE;
                }

              // unrefine check
              d01 = shock_diff(0, m_eMesh, node0, node1, coord0, coord1, m_shock, m_shock_displacement);
              if ( d01 <= m_shock_diff_criterion/2.0 )  // 0.05, 0.2
                {
                  mark |= DO_UNREFINE;
                }
            }
          return mark;
        }

        //double *fdata = eMesh.field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        //return m_selector(entity) && fdata[0] > 0;
      };

      // This can be used as an edge or element-based predicate

      struct ShockBasedRefinePredicate1 : public IEdgeBasedAdapterPredicate, IElementBasedAdapterPredicate, percept::ElementOp,  GenericFunction {

        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase * m_nodal_refine_field;
        PlaneShock m_shock;
        double m_shock_displacement;
        double m_shock_diff_criterion;

        ShockBasedRefinePredicate1(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh, stk::mesh::Selector* selector, stk::mesh::FieldBase *field, double tolerance,
                                  PlaneShock shock, double shock_displacement=0, double shock_diff_criterion=0.4) :
          IEdgeBasedAdapterPredicate(eMesh, selector, field, tolerance),
          IElementBasedAdapterPredicate(eMesh, selector, field, tolerance),
          m_eMesh(eMesh),m_nodal_refine_field(nodal_refine_field),
          m_shock(shock), m_shock_displacement(shock_displacement),
          m_shock_diff_criterion(shock_diff_criterion) {}


        /// Element-based: Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity element)
        {
          if (0 == m_eb_selector || (*m_eb_selector)(m_eMesh.bucket(element)))
            {

              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
              int spatialDimension = m_eMesh.get_spatial_dim();
              CellTopology cell_topo(cell_topo_data);
              const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);

              CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

              unsigned numSubDimNeededEntities = 0;
              numSubDimNeededEntities = cell_topo_data->edge_count;

              unsigned ref_count=0, unref_count=0;
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = m_eMesh.field_data( *coordField , node0 );
                  double * const coord1 = m_eMesh.field_data( *coordField , node1 );
                  double  dcoord0[3] = {coord0[0],coord0[1], (spatialDimension==2?0:coord0[2])};
                  double  dcoord1[3] = {coord1[0],coord1[1], (spatialDimension==2?0:coord1[2])};

                  int markInfo = this->operator()(element, iSubDimOrd, node0, node1, dcoord0, dcoord1, 0);
                  bool do_ref = markInfo & DO_REFINE;
                  if (do_ref)
                    {
                      ++ref_count;
                    }
                  bool do_unref = markInfo & DO_UNREFINE;
                  if (!do_ref && do_unref)
                    {
                      ++unref_count;
                    }
                }
              if (ref_count == numSubDimNeededEntities)
                return DO_REFINE;
              if (unref_count == numSubDimNeededEntities)
                return DO_UNREFINE;
            }

          return DO_NOTHING;
        }

        /// Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
                        double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
        {
          int mark=0;
          if (0 == m_selector || (*m_selector)(m_eMesh.bucket(element)))
            {

              // refine check
              mark = shock_diff1(m_nodal_refine_field, m_eMesh, node0, node1, coord0, coord1, m_shock, m_shock_displacement);
            }
          return mark;
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK );
          //stk::mesh::Entity const *elem_nodes_e = m_eMesh.get_bulk_data()->begin_entities(element, stk::topology::NODE_RANK );

          unsigned num_node = elem_nodes.size();
          //unsigned num_node = m_eMesh.get_bulk_data()->num_connectivity(element, stk::topology::NODE_RANK );
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

          bool found = false;
          for (unsigned inode=0; inode < num_node-1; inode++)
            {
              stk::mesh::Entity node_i = elem_nodes[ inode ].entity();
              //mesh::Entity node_i = elem_nodes_e[ inode ];
              double *coord_data_i = m_eMesh.field_data(coordField, node_i);

              for (unsigned jnode=inode+1; jnode < num_node; jnode++)
                {
                  stk::mesh::Entity node_j = elem_nodes[ jnode ].entity();
                  double *coord_data_j = m_eMesh.field_data(coordField, node_j);

                  int mark = shock_diff1(m_nodal_refine_field, m_eMesh, node_i, node_j, coord_data_i, coord_data_j, m_shock, m_shock_displacement);
                  if (mark & DO_REFINE)
                    {
                      found=true;
                      break;
                    }
                }
            }
          if (found)
            f_data[0] = 1;
          else
            f_data[0] = -1;

          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}

        // nodal op
        virtual void operator()(MDArray& domain, MDArray& codomain, double time_value_optional=0.0)
        {
          double x = domain(0);
          double y = domain(1);
          int spatialDimension = m_eMesh.get_spatial_dim();
          double z = (spatialDimension==3?domain(2):0);
          m_shock.setCurrentPlanePoint(m_shock_displacement);
          double *plane_point = m_shock.plane_point;
          double *plane_normal = m_shock.plane_normal;

          double coord0[] = {x,y,z};
          double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
          codomain(0) = dot_0 < 0? -1 : 1;
        }

      };

      class SetRefineField1 : public percept::ElementOp
      {
        IEdgeAdapter& m_breaker;
      public:
        SetRefineField1(IEdgeAdapter& breaker) : m_breaker(breaker) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          bool isParent = m_breaker.getMesh().isParentElement(element, false);
          bool hasFamilyTree = m_breaker.getMesh().hasFamilyTree(element);
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          int mark = m_breaker.markUnrefine(element);
          f_data[0] = 0;
          if (hasFamilyTree && !isParent)
            {
              if (mark & DO_UNREFINE)
                f_data[0] = -2;
              else
                f_data[0] = 2;
            }

          return false;
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}
      };

      class SetElementFieldBoundingBox : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase *m_field;
      public:
        std::vector<double> m_bbox;
        SetElementFieldBoundingBox(percept::PerceptMesh& eMesh, const std::vector<double>& bbox, stk::mesh::FieldBase *field) : m_eMesh(eMesh), m_field(field), m_bbox(bbox) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          stk::mesh::FieldBase* coord_field = m_eMesh.get_coordinates_field();

          const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
          //int spatialDimension = m_eMesh.get_spatial_dim();
          CellTopology cell_topo(cell_topo_data);

          const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );

          unsigned num_node = elem_nodes.size();
          double c[] = {0,0,0};
          for (unsigned inode=0; inode < num_node; inode++)
            {
              stk::mesh::Entity node = elem_nodes[ inode ].entity();
              double *c_data = m_eMesh.field_data(coord_field, node);
              c[0] += c_data[0]/double(num_node);
              c[1] += c_data[1]/double(num_node);
              c[2] += c_data[2]/double(num_node);
            }
          std::vector<double> cc(&c[0], &c[0]+3);
          if (0) std::cout << "SetElementFieldBoundingBox element " << m_eMesh.identifier(element) << " topo= " << cell_topo.getName() 
                    << " c= " << cc << " bbox= " << m_bbox
                    <<  std::endl;

          if (
              (m_bbox[0] < c[0] && c[0] < m_bbox[1]) &&
              (m_bbox[2] < c[1] && c[1] < m_bbox[3]) &&
              (m_bbox[4] < c[2] && c[2] < m_bbox[5])
              )
            {
              if (0) std::cout << "SetElementFieldBoundingBox marking element " << m_eMesh.identifier(element) << " for refine" << std::endl;
              f_data[0] = 1;
            }
          else
            {
              f_data[0] = -1;
            }
          return false;  // don't terminate the loop
        }
        virtual void init_elementOp()
        {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }
        virtual void fini_elementOp() {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }
      };

      class SetElementFieldSelector : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase *m_field;
        stk::mesh::Selector m_selector;
      public:
        SetElementFieldSelector(percept::PerceptMesh& eMesh, stk::mesh::Selector sel, stk::mesh::FieldBase *field) : m_eMesh(eMesh), m_field(field), m_selector(sel) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          if (m_selector(m_eMesh.bucket(element)))
            {
              f_data[0] = 1;
            }
          else
            {
              f_data[0] = -1;
            }
          return false;  // don't terminate the loop
        }
        virtual void init_elementOp()
        {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }
        virtual void fini_elementOp() {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }
      };


      struct PointSourceRefinePredicate : public IElementBasedAdapterPredicate
      {

        percept::PerceptMesh& m_eMesh;
        std::vector<double> m_pointSourcePosition;

        PointSourceRefinePredicate(percept::PerceptMesh& eMesh, stk::mesh::Selector* selector, stk::mesh::FieldBase *field, double tolerance,
                                   std::vector<double>& pointSourcePosition) :
        IElementBasedAdapterPredicate(eMesh, selector, field, tolerance),
          m_eMesh(eMesh), m_pointSourcePosition(pointSourcePosition)
        {
        }

        /// Element-based: Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity element)
        {
          if (0 == m_eb_selector || (*m_eb_selector)(m_eMesh.bucket(element)))
            {
              const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
              //int spatialDimension = m_eMesh.get_spatial_dim();
              CellTopology cell_topo(cell_topo_data);
              const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
              VERIFY_OP_ON(elem_nodes.size(), ==, 3, "only for tris");

              CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

              unsigned numSubDimNeededEntities = 0;
              numSubDimNeededEntities = cell_topo_data->edge_count;

              //int count=0;
              unsigned ref_count=0, unref_count=0;
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = m_eMesh.field_data( *coordField , node0 );
                  double * const coord1 = m_eMesh.field_data( *coordField , node1 );
                  double cross =
                    (coord1[0] - coord0[0])*(m_pointSourcePosition[1] - coord0[1])
                    - (coord1[1] - coord0[1])*(m_pointSourcePosition[0] - coord0[0]);

                  if (cross > 0)
                    ++ref_count;
                  else
                    ++unref_count;
                }
              if (ref_count == numSubDimNeededEntities)
                return DO_REFINE;
              else
                return DO_UNREFINE;

            }

          return DO_NOTHING;
        }

      };
    }

    struct EdgeBasedInterpError : public IEdgeBasedAdapterPredicate {

    	  EdgeBasedInterpError(PerceptMesh& eMesh, stk::mesh::Selector * selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
          IEdgeBasedAdapterPredicate(eMesh, selector, field, tolerance) {}

        /// Return DO_NOTHING, DO_REFINE, DO_UNREFINE or sum of these
        int operator()(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
                       double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
        {
          int mark = DO_NOTHING;

          if (0 == m_selector || (*m_selector)(m_eMesh.bucket(element)))
            {
              double edge_vec[2];
              edge_vec[0] = coord0[0] - coord1[0];
              edge_vec[1] = coord0[1] - coord1[1];

              double midpoint[2];
              midpoint[0] = 0.5*(coord0[0] + coord1[0]);
              midpoint[1] = 0.5*(coord0[1] + coord1[1]);

              double Hess[2][2];
              memset(&Hess[0][0], 0, sizeof(double)*4);

              // assume u(x,y) = x*x*x
              Hess[0][0]= 6.0*midpoint[0]; // 6*x

              double error = 0;

              for (int i=0; i<2; i++) {
            	  for (int j=0; j<2; j++) {
            		  error += Hess[i][j] * edge_vec[i] * edge_vec[j];
            	  }
              }

              if (error > 0.05) {
            	  mark |= DO_REFINE;
              }
            }

          return mark;
        }
      };

  }

#endif
