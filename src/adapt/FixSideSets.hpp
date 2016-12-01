// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef FixSideSets_hpp
#define FixSideSets_hpp

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <percept/Percept.hpp>

#include <adapt/Refiner.hpp>
#include <adapt/RefinerUtil.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/MeshUtil.hpp>
#include <adapt/AdaptedMeshVerifier.hpp>

#define DEBUG_FSS 0

namespace percept {


  class FixSideSets {
  public:
    PerceptMesh& m_eMesh;
    stk::mesh::PartVector& m_excludeParts;
    SidePartMap& m_side_part_map;
    std::string m_geomFile;
    bool m_avoidFixSideSetChecks;

    FixSideSets(PerceptMesh& eMesh, stk::mesh::PartVector& excludeParts, SidePartMap& side_part_map, const std::string& geomFile, bool avoidFixSideSetChecks)
      : m_eMesh(eMesh), m_excludeParts(excludeParts), m_side_part_map(side_part_map), m_geomFile(geomFile), m_avoidFixSideSetChecks(avoidFixSideSetChecks) {}


    bool is_perm_bad(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord, stk::mesh::Permutation& perm);
    bool is_positive_perm(stk::mesh::Entity element, stk::mesh::Entity side, unsigned side_ord);
    bool has_default_perm(stk::mesh::Entity side);
    void fix_permutation(SetOfEntities& side_set);
    bool connect(stk::mesh::Entity side, bool& valid_side_part_map, SetOfEntities* avoid_elems, bool onlyPosPerm=false);

    // if the element (element) has a side that matches  the given side (side), connect them but first delete old connections
    bool connectSidesForced(stk::mesh::Entity element, stk::mesh::Entity side, bool& valid_side_part_map, bool onlyPosPerm = false);
    void disconnect_entity(stk::mesh::Entity entity);

    template<class SetOfEntities>
    static void print_set(PerceptMesh& eMesh, SetOfEntities& side_set, const std::string& msg = "sides_to_remove.size")
    {
      if (DEBUG_FSS)
        {
          std::cout << "print_set: " << msg << " side_set.size= " << side_set.size();
          for (typename SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
            {
              std::cout << " " << eMesh.id(*it_side); // << " E: " << eMesh.print_entity_compact(*it_side)
              //<< " nc: " << eMesh.numChildren(*it_side);
              ;
            }
          std::cout << std::endl;
        }
    }

    void delete_unattached_sides(SetOfEntities& side_set, SetOfEntities *avoid_sides);
    bool bucket_acceptable(stk::mesh::Bucket& bucket, stk::mesh::EntityRank rank);
    void build_side_set(SetOfEntities& side_set);
    void disconnect_sides(SetOfEntities& side_set);
    void reconnect_sides(SetOfEntities& side_set, SetOfEntities *avoid_elems, bool onlyPosPerm);
    void check_connect(SetOfEntities& side_set, SetOfEntities *avoid_elems);
    void end_begin();

    // fast reconnector

    void
    fix_side_sets_2(bool allow_not_found, SetOfEntities *avoid_elems, SetOfEntities *avoid_sides);
  };


}
#endif
