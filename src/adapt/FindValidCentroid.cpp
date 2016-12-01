// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/FindValidCentroid.hpp>

#if defined(STK_BUILT_IN_SIERRA) && !STK_PERCEPT_LITE
#include <percept/mesh/geometry/volume/sierra_only/FiniteVolumeMesh.hpp>
#endif

namespace percept {

  double FindValidCentroid::getVolumes(std::vector<double>& volumes, stk::mesh::Entity element)
  {
    volumes.resize(0);
    const CellTopologyData *elem_cell_topo_data = m_eMesh.get_cell_topology(element);
    double element_volume = m_eMesh.volume(element, m_eMesh.get_coordinates_field(), elem_cell_topo_data);
    if (m_eMesh.numChildren(element) == 0)
      return element_volume;
    std::vector<stk::mesh::Entity> children;
    m_eMesh.getChildren(element, children, false, false);
    volumes.resize(children.size());
    for (unsigned ii=0; ii < children.size(); ++ii)
      {
        const CellTopologyData *cell_topo_data = m_eMesh.get_cell_topology(children[ii]);
        volumes[ii] = m_eMesh.volume(children[ii], m_eMesh.get_coordinates_field(), cell_topo_data);
        if (m_debug) std::cout << "ii= " << ii << " child= " << m_eMesh.print_entity_compact(children[ii]) << " volume= " << volumes[ii] << std::endl;
        if (0)
          {
            volumes[ii] = std::numeric_limits<double>::max();
            VolumeUtil jacA;
            shards::CellTopology cell_topo(cell_topo_data);
            double volScale = jacA.getJacobianToVolumeScale(cell_topo);
            double jacobian = 0.0;
            jacA(jacobian, m_eMesh, children[ii], m_eMesh.get_coordinates_field(), cell_topo_data);
            double cellVol = jacobian*volScale;
            for (int in=0; in < jacA.m_num_nodes; ++in)
              {
                if (m_debug) std::cout << "ii= " << ii << " in= " << in << " v[in]= " << jacA.m_detJ[in] << " cellVol= " << cellVol << std::endl;
                volumes[ii] = std::min(volumes[ii], jacA.m_detJ[in]*volScale);
              }
          }
#if defined(STK_BUILT_IN_SIERRA) && !STK_PERCEPT_LITE
        if (m_use_finite_volume)
          {
            volumes[ii] = std::numeric_limits<double>::max();
            double sc_volume[8];
            FiniteVolumeMesh3D fvm(*m_eMesh.get_bulk_data());

            shards::CellTopology cell_topo(cell_topo_data);
            unsigned numNodes = cell_topo.getNodeCount();
            fvm.elementVolume(children[ii], sc_volume);
            double cellVol = 0.0;
            for (unsigned in=0; in < numNodes; ++in)
              {
                cellVol += sc_volume[in];
              }
            for (unsigned in=0; in < numNodes; ++in)
              {
                if (m_debug) std::cout << "ii= " << ii << " in= " << in << " v[in]= " << sc_volume[in] << " cellVol= " << cellVol << std::endl;
                volumes[ii] = std::min(volumes[ii], sc_volume[in]);
              }
          }
#endif
      }
    return element_volume;
  }

  double FindValidCentroid::metric(std::vector<double>& volumes, bool& foundBad)
  {
    double met = 0.0;
    foundBad = false;
    for (unsigned ii=0; ii < volumes.size(); ++ii)
      {
        if (volumes[ii] <= 0.0)
          {
            met += -volumes[ii];
            foundBad = true;
            break;
          }
      }
    if (foundBad)
      return met;

    met = 0.0;
    for (unsigned ii=0; ii < volumes.size(); ++ii)
      {
        met += 1.0/volumes[ii];
      }
    return met;
  }

  // return if changed
  bool FindValidCentroid::findCentroid(stk::mesh::Entity element, double *c_p, std::vector<stk::mesh::Entity>& nodes, stk::mesh::Entity c_node)
  {
    // only for coordinate field
    int fieldDim = 3;
    stk::mesh::FieldBase *field = m_eMesh.get_coordinates_field();
    unsigned nsz = nodes.size();
    double *c_node_p = (double *)stk::mesh::field_data(*field, c_node);
    double c_p_save[3] = {c_p[0], c_p[1], c_p[2]};
    if (0)
      for (int isp = 0; isp < fieldDim; isp++)
        {
          c_node_p[isp] = c_p[isp];
        }
    if (m_debug) std::cout << "element= " << m_eMesh.print_entity_compact(element) << std::endl;


    double c_node_p_save[3] = {c_node_p[0], c_node_p[1], c_node_p[2]};
    (void)c_node_p_save;

    // check volumes
    stk::topology topo = m_eMesh.topology(element);
    if (topo == stk::topology::PYRAMID_5 || topo == stk::topology::WEDGE_6 || topo == stk::topology::HEX_8)
      {
        std::vector<double> volumes;
        bool foundBad = false;
        double element_volume = getVolumes(volumes, element);
        (void)element_volume;
        double met = metric(volumes, foundBad);
        if (m_debug) std::cout << "element= " << m_eMesh.print_entity_compact(element) << " initial met= " << met << " foundBad= " << foundBad << " element_volume= " << element_volume << std::endl;
        if (!foundBad)
          {
            if (0)
              for (int isp = 0; isp < fieldDim; isp++)
                {
                  c_p[isp] = c_node_p[isp];
                }
            return false;
          }
        if (m_debug) std::cout << "element= " << m_eMesh.print_entity_compact(element) << " initial met= " << met << " foundBad= " << foundBad << " element_volume= " << element_volume << std::endl;

        double metMin = std::numeric_limits<double>::max();
        double c_p_min[3] = {0,0,0};
        bool foundGood = false;
        for (int ix = 1; ix < ndiv; ++ix)
          {
            double xi = double(ix)/double(ndiv);
            for (int iy = 1; iy < ndiv; ++iy)
              {
                double eta = double(iy)/double(ndiv);
                for (int iz = 1; iz < ndiv; ++iz)
                  {
                    double zeta = double(iz)/double(ndiv);

                    double hexBases[] = {
                      (1-xi)*(1-eta)*(1-zeta),
                      (  xi)*(1-eta)*(1-zeta),
                      (  xi)*(  eta)*(1-zeta),
                      (1-xi)*(  eta)*(1-zeta),
                      (1-xi)*(1-eta)*(  zeta),
                      (  xi)*(1-eta)*(  zeta),
                      (  xi)*(  eta)*(  zeta),
                      (1-xi)*(  eta)*(  zeta)
                    };
                    double pyrBases[] = {
                      (1-xi)*(1-eta)*(1-zeta),
                      (  xi)*(1-eta)*(1-zeta),
                      (  xi)*(  eta)*(1-zeta),
                      (1-xi)*(  eta)*(1-zeta),
                      zeta
                    };
                    double wedgeBases[] = {
                      (1-xi-eta)*(1-zeta),
                      (  xi    )*(1-zeta),
                      (     eta)*(1-zeta),
                      (1-xi-eta)*(  zeta),
                      (  xi    )*(  zeta),
                      (     eta)*(  zeta)
                    };

                    double *bases = 0;
                    switch(topo.value()) {
                    case stk::topology::PYRAMID_5:
                      bases = &pyrBases[0];
                      break;
                    case stk::topology::WEDGE_6:
                      bases = &wedgeBases[0];
                      break;
                    case stk::topology::HEX_8:
                      bases = &hexBases[0];
                      break;
                    default:
                      VERIFY_MSG("bad topo");
                    }

                    for (int isp = 0; isp < fieldDim; isp++)
                      {
                        c_node_p[isp] = 0.0;
                      }

                    for (unsigned ipts=0; ipts < nsz; ipts++)
                      {
                        stk::mesh::Entity node = nodes[ipts];
                        double *  field_data = m_eMesh.field_data_inlined(field, node);
                        if (field_data)
                          {
                            for (int isp = 0; isp < fieldDim; isp++)
                              {
                                c_node_p[isp] += field_data[isp]*bases[ipts];
                              }
                          }
                      }

                    foundBad = false;
                    element_volume = getVolumes(volumes, element);
                    met = metric(volumes, foundBad);
                    if (m_debug) std::cout << "ix,y,z= " << ix << " " << iy << " " << iz << " xi= " << xi << " " << eta << " " << zeta << " met= " << met << " foundBad= " << foundBad << std::endl;
                    if (!foundBad)
                      {
                        foundGood = true;
                        if (met < metMin)
                          {
                            metMin = met;
                            for (int isp = 0; isp < fieldDim; isp++)
                              {
                                c_p_min[isp] = c_node_p[isp];
                              }
                          }
                      }

                  }
              }
          }

        if (foundGood)
          {
            for (int isp = 0; isp < fieldDim; isp++)
              {
                c_node_p[isp] = c_p_min[isp];
                c_p[isp] = c_p_min[isp];
              }
            if (m_debug)
              {
                foundBad = false;
                element_volume = getVolumes(volumes, element);
                met = metric(volumes, foundBad);
                std::cout << "foundGood c_node= " << m_eMesh.print_entity_compact(c_node) << " element= " << m_eMesh.id(element) << " met= " << met << " foundBad= " << foundBad << std::endl;
                VERIFY_OP_ON(foundBad, ==, false, "bad foundBad 2");
              }
          }
        else
          {
            std::ostringstream str;
            if (0 && !m_debug)
              {
                for (int isp = 0; isp < fieldDim; isp++)
                  {
                    c_node_p[isp] = c_node_p_save[isp];
                    c_p[isp] = c_p_save[isp];
                  }
                m_debug = true;
                findCentroid(element, c_p, nodes, c_node);
              }
            m_debug = true;
            element_volume = getVolumes(volumes, element);
            str << "negative volumes in refined mesh - couldn't find centroid, m_use_finite_volume= " << m_use_finite_volume << " parent topo= " << topo << " element= " << m_eMesh.print_entity_compact(element);

            if (1)
              {
                std::set<stk::mesh::Entity> ll;
                ll.insert(element);
                std::vector<stk::mesh::Entity> children;
                m_eMesh.getChildren(element, children, false, false);
                for (unsigned ii=0; ii < children.size(); ++ii)
                  {
                    ll.insert(children[ii]);
                  }
                std::string ff = "bad-fvc-"+toString(m_eMesh.get_rank())+".vtk";
                m_eMesh.dump_vtk(ff, false, &ll);
              }
            VERIFY_MSG(str.str());
          }

      }
    return true;
  }

}
