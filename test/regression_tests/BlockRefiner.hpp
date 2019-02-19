// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_BlockRefiner_hpp
#define adapt_BlockRefiner_hpp

#include <percept/PerceptMesh.hpp>
#include <percept/FieldTypes.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/ElementRefinePredicate.hpp>
#include <adapt/AdaptHelperFunctions.hpp>

  namespace percept {

    class BlockRefiner {
      PerceptMesh m_eMesh;
      BlockNamesType m_block_names;
      std::string m_input_file, m_ouput_file;
      bool m_verify;
      bool m_plot;
    public:
      //BlockRefiner(PerceptMesh& eMesh) : m_eMesh(eMesh) {}
      BlockRefiner(const std::string& input, const std::string& output, bool verify=true, bool plot=false) : m_input_file(input), m_ouput_file(output), m_verify(verify), m_plot(plot) {}

      void set_refine_field(stk::mesh::Selector& block_selector)
      {
        PerceptMesh& eMesh = m_eMesh;
        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            RefineFieldType_type val = static_cast<RefineFieldType_type>(0);
            if (block_selector(bucket))
              val = static_cast<RefineFieldType_type>(1);

            if (1)
              {
                const unsigned num_elements_in_bucket = bucket.size();
                for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                  {
                    stk::mesh::Entity element = bucket[iElement];
                    RefineFieldType_type *fdata = stk::mesh::field_data( *eMesh.m_refine_field , element );
                    fdata[0] = static_cast<RefineFieldType_type>(val);
                  }
              }
          }
      }

      void refine(const std::string& block_names_str, int nrefs=1)
      {
        PerceptMesh& eMesh = m_eMesh;
        eMesh.open(m_input_file);

        BlockNamesType block_names(percept::EntityRankEnd+1u);

        if (1)
          {
            std::cout << "block_names: original block_name option = " << block_names_str << std::endl;
            block_names = RefinerUtil::getBlockNames(block_names_str, (unsigned)eMesh.get_rank(), eMesh);
            std::cout << "block_names: after getBlockNames block_name option = " << block_names << std::endl;
            if (1)
              {
                eMesh.commit();
                std::string input_geometry = ""; // FIXME
                block_names = RefinerUtil::correctBlockNamesForPartPartConsistency(eMesh, block_names, input_geometry);

                eMesh.close();
                eMesh.open(m_input_file);
              }
          }
        m_block_names = block_names;

        eMesh.register_and_set_refine_fields();

        Teuchos::RCP<UniformRefinerPatternBase> localBreakPattern = make_local_break_pattern(eMesh);

        eMesh.commit();

        std::set<stk::mesh::Part *> pl;
        if (1)
          {
            for (unsigned ii=0; ii < block_names[eMesh.element_rank()].size(); ++ii)
              {
                std::string bn = block_names[eMesh.element_rank()][ii];
                VERIFY_OP_ON ((bn[0] == '+' || bn[0] == '-'), ==, true, "bad block name: "+bn);
                std::string bname = bn.substr(1);
                stk::mesh::Part *part = eMesh.get_fem_meta_data()->get_part(bname);
                VERIFY_OP_ON(part, !=, 0, "couldn't find part: "+bname);

                if (bn[0] == '+')
                  {
                    pl.insert(part);
                  }
              }
            for (unsigned ii=0; ii < block_names[eMesh.element_rank()].size(); ++ii)
              {
                std::string bn = block_names[eMesh.element_rank()][ii];
                std::string bname = bn.substr(1);
                stk::mesh::Part *part = eMesh.get_fem_meta_data()->get_part(bname);

                if (bn[0] == '-')
                  {
                    if (pl.find(part) != pl.end())
                      pl.erase(part);
                  }
              }
          }
        stk::mesh::PartVector pv(pl.begin(), pl.end());
        if (!eMesh.get_rank()) std::cout << "block_names= " << block_names << "\nfrom parts = " << eMesh.print_part_vector_string(pv) << std::endl;

        stk::mesh::Selector block_selector = stk::mesh::selectUnion( pv );

        eMesh.output_active_children_only(true);

        std::cout << "block_selector= " << block_selector << std::endl;
        bool debug = false;
        AdaptedMeshVerifier adaptedMeshVerifier(debug);
        if (m_verify && !adaptedMeshVerifier.isValid(eMesh, true))
          VERIFY_MSG("Invalid initial mesh");

        stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());
        ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);
        TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, *localBreakPattern, 0);

        breaker.setRemoveOldElements(false);
        breaker.setAlwaysInitializeNodeRegistry(false);

        int iplot=0;
        if (m_plot)
          {
            char buf[1000];
            sprintf(buf, "%04d", iplot);
            if (iplot == 0)
              eMesh.save_as(m_ouput_file+"-anim.e");
            else
              eMesh.save_as(m_ouput_file+"-anim.e-s"+std::string(buf));
            ++iplot;
          }

        for (int ir = 0; ir < nrefs; ++ir)
          {
            if (!eMesh.get_rank())
              std::cout << "Refinement pass # " << ir << " start..." << std::endl;
            set_refine_field(block_selector);
            breaker.refine();
            if (m_plot)
              {
                char buf[1000];
                sprintf(buf, "%04d", iplot);
                if (iplot == 0)
                  eMesh.save_as(m_ouput_file+"-anim.e");
                else
                  eMesh.save_as(m_ouput_file+"-anim.e-s"+std::string(buf));
                ++iplot;
              }
          }

        if (m_verify && !adaptedMeshVerifier.isValid(eMesh, false))
          VERIFY_MSG("Invalid adapted mesh");

        eMesh.save_as(m_ouput_file);
      }


    };

  }

#endif
