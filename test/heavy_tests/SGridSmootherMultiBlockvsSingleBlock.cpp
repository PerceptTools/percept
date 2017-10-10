// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/Percept.hpp>
#if defined( STK_PERCEPT_HAS_GEOMETRY )

#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherAlgebraic.hpp>
#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_total_element_metric.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_update_coordinates.hpp>
#include <percept/mesh/mod/smoother/gradient_functors.hpp>
#include <percept/mesh/mod/smoother/get_alpha_0_refmesh.hpp>
#include <percept/mesh/mod/smoother/get_edge_len_avg.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>
#include <percept/structured/StructuredGridRefiner.hpp>

#include <iostream>
#include <cstdlib>
#include <array>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

#define DO_TESTS 1
#if DO_TESTS

#include <percept/PerceptUtils.hpp>
#include <percept/Util.hpp>

#include <percept/math/DenseMatrix.hpp>

namespace std {
}

namespace percept
{

struct zero_a4d
{
    StructuredGrid::MTField::Array4D a4d;

    zero_a4d(StructuredGrid::MTField::Array4D a4d_in){a4d=a4d_in;}

    void zero_out_fields()
    {
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,(unsigned)a4d.size()),*this);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const
    {
        a4d.data()[index]=0;
    }

};

void is_field_data_equivalent(std::string fld_name,unsigned noElemsPerDim, unsigned spatialDim, std::shared_ptr<BlockStructuredGrid> bsg_mb, std::shared_ptr<BlockStructuredGrid> bsg_sb)
{
    StructuredGrid::MTField *coord_field_mb = bsg_mb->m_fields["coordinates_saved"].get();

    if(!coord_field_mb)
    {
        std::cout << "no coordinates_saved field, aborting the diff\n";
        return;
    }

    Array4D sb_ela_fld =  *(bsg_sb->m_fields[fld_name].get()->m_block_fields[0]);

    double data[3];
    for (unsigned iblock=0; iblock < bsg_mb->m_sblocks.size(); ++iblock)
    {
        Array4D coord_field_iblock = *(coord_field_mb->m_block_fields[iblock]);
        Array4D::HostMirror host_coord_field_iblock = Kokkos::create_mirror_view(coord_field_iblock);
        Kokkos::deep_copy(host_coord_field_iblock,coord_field_iblock);
        const SGridSizes sgridsizes = bsg_mb->m_sblocks[iblock]->m_sizes;
        const std::array<unsigned int, 3> loop_orderings = {{bsg_mb->m_sblocks[iblock]->m_loop_ordering[0],bsg_mb->m_sblocks[iblock]->m_loop_ordering[1],bsg_mb->m_sblocks[iblock]->m_loop_ordering[2]}};

        Array4D mb_ela_fld =  *(bsg_mb->m_fields[fld_name].get()->m_block_fields[iblock]);

        unsigned nNodes = (1+sgridsizes.node_max[0])*(1+sgridsizes.node_max[1])*(1+sgridsizes.node_max[2]);
        for(unsigned iNode=0;iNode<nNodes;iNode++) {

            unsigned indx[3];
            const int L0 = loop_orderings[0], L1 =
                    loop_orderings[1], L2 = loop_orderings[2];
            const unsigned int sgsizes[3] = { 1 + sgridsizes.node_max[L0]
                                                                      - sgridsizes.node_min[L0], 1
                                                                      + sgridsizes.node_max[L1]
                                                                                            - sgridsizes.node_min[L1], 1
                                                                                            + sgridsizes.node_max[L2]
                                                                                                                  - sgridsizes.node_min[L2] };
            indx[L2] = sgridsizes.node_min[L2]
                                           + (iNode / (sgsizes[L0] * sgsizes[L1]));
            indx[L1] = sgridsizes.node_min[L1]
                                           + ((iNode / sgsizes[L0]) % sgsizes[L1]);
            indx[L0] = sgridsizes.node_min[L0] + (iNode % sgsizes[L0]);

            unsigned i=indx[0];
            unsigned j=indx[1];
            unsigned k=indx[2];

            for (int d=0; d<3; d++)
                data[d]=host_coord_field_iblock(i,j,k,d);


            int org_ijk[3];

            org_ijk[0] = (int)((data[0]*noElemsPerDim) + 0.5);
            org_ijk[1] = (int)((data[1]*noElemsPerDim) + 0.5);
            org_ijk[2] = (int)((data[2]*noElemsPerDim) + 0.5);

            for(unsigned iDim=0;iDim<spatialDim;iDim++)
                EXPECT_NEAR(mb_ela_fld(i,j,k,iDim),sb_ela_fld(org_ijk[0],org_ijk[1],org_ijk[2],iDim),1e-12);
        }//foreach node
    }
}

void build_perturbed_cube_multiblock(PerceptMesh& eMesh, const unsigned nele,PerceptMesh& eMesh_2 , std::string filename , bool add_all_fields, bool refine = false, unsigned noRefs = 1)
{
    const unsigned nxyz = nele+1;
    std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
    std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);

    unsigned noElems = nele;
    if(refine)
    {
        unsigned print_level = 0;
        for(unsigned iRef=0;iRef<noRefs;iRef++)
        {
            StructuredGridRefiner sgr(bsg, print_level);
            sgr.do_refine();
            noElems *= 2; //double the number of elements in each direction
            for(unsigned i=0;i<3;i++)
                sizes[i] = noElems + 1;
            bsg = sgr.m_output;
        }
    }

    eMesh.set_block_structured_grid(bsg);

    if(!add_all_fields) {
        bsg->register_field("coordinates", eMesh.get_spatial_dim());
        bsg->register_field("coordinates_N", eMesh.get_spatial_dim());
        bsg->register_field("coordinates_NM1", eMesh.get_spatial_dim());
        bsg->register_field("cg_g",eMesh.get_spatial_dim());
    }
    else{
        bsg->register_field("coordinates", eMesh.get_spatial_dim());
        bsg->register_field("coordinates_saved", eMesh.get_spatial_dim());
        eMesh.copy_field("coordinates_saved","coordinates");
        eMesh.add_coordinate_state_fields();
    }

    typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
    VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
    typename StructuredGrid::MTField *coord_field_N = bsg->m_fields["coordinates_N"].get();
    VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
    typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
    VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

    // save state of original mesh
    // field, dst, src:
    eMesh.copy_field(coord_field_NM1, coord_field);

    Array4D coord_field_0 = *(coord_field->m_block_fields[0]);//hard coded to single block case
    Array4D::HostMirror host_coord_field = Kokkos::create_mirror_view(coord_field_0);
    Kokkos::deep_copy(host_coord_field,coord_field_0);

    // randomly perturb interior nodes
    const double neleCoeff = 1.0;//increase in order to observe smaller difference between multi and single block smootherlog files
    const double epsilon = 1.0/((double)nele*neleCoeff);
    const double tol = 1e-6;
    std::srand(1);
    double data[3];

    unsigned noManipedCoord1block = 0;
    for (unsigned k=0; k<sizes[0]; k++) {
        for (unsigned j=0; j<sizes[1]; j++) {
            for (unsigned i=0; i<sizes[2]; i++) {
                for (int d=0; d<3; d++) {
                    data[d]=host_coord_field(i,j,k,d);
                }

                const double x = data[0], y = data[1], z = data[2];
                if (x > tol && x< 1-tol &&
                        y > tol && y< 1-tol &&
                        z > tol && z< 1-tol) {

                    for (int d=0; d<3; d++) {
                        double rand_coeff = (double)std::rand()/(double)RAND_MAX;
                        data[d] = data[d] +epsilon*( rand_coeff - 0.5);
                    }

                    noManipedCoord1block++;
                }//tolerance
                for (int d=0; d<3; d++) {
                    host_coord_field(i,j,k,d)=data[d];
                }
            }
        }
    }

    Kokkos::deep_copy(coord_field_0, host_coord_field);

    // save state of projected mesh
    // field, dst, src:
    eMesh.copy_field(coord_field_N, coord_field);

    std::string dbtype("cgns_structured");
    eMesh_2.open(filename, dbtype);
    eMesh_2.readBulkData();
    std::shared_ptr<BlockStructuredGrid> bsg_2 = eMesh_2.get_block_structured_grid();

    if(refine)
    {
        unsigned print_level = 0;
        for(unsigned iRef=0;iRef<noRefs;iRef++)
        {
            StructuredGridRefiner sgr(eMesh_2.get_block_structured_grid(), print_level);
            sgr.do_refine();
            bsg_2 = sgr.m_output;
            eMesh_2.set_block_structured_grid(bsg_2);
        }
    }

    if(!add_all_fields) {
        bsg_2->register_field("coordinates", eMesh_2.get_spatial_dim());
        bsg_2->register_field("coordinates_N", eMesh_2.get_spatial_dim());
        bsg_2->register_field("coordinates_NM1", eMesh_2.get_spatial_dim());
        bsg_2->register_field("cg_g",eMesh_2.get_spatial_dim());
    }
    else{
        bsg_2->register_field("coordinates", eMesh_2.get_spatial_dim());
        bsg_2->register_field("coordinates_saved", eMesh.get_spatial_dim());
        eMesh_2.copy_field("coordinates_saved","coordinates");
        eMesh_2.add_coordinate_state_fields();
    }

    StructuredGrid::MTField *coord_field_2 = bsg_2->m_fields["coordinates"].get();
    StructuredGrid::MTField *coord_field_N_2 = bsg_2->m_fields["coordinates_N"].get();
    StructuredGrid::MTField *coord_field_NM1_2 = bsg_2->m_fields["coordinates_NM1"].get();
    eMesh_2.copy_field(coord_field_NM1_2, coord_field_2);

    for (unsigned iblock=0; iblock < bsg_2->m_sblocks.size(); ++iblock)
    {
        Array4D coord_field_iblock = *(coord_field_2->m_block_fields[iblock]);
        Array4D::HostMirror host_coord_field_iblock = Kokkos::create_mirror_view(coord_field_iblock);
        Kokkos::deep_copy(host_coord_field_iblock,coord_field_iblock);
        const SGridSizes sgridsizes = bsg_2->m_sblocks[iblock]->m_sizes;
        const std::array<unsigned int, 3> loop_orderings = {{bsg_2->m_sblocks[iblock]->m_loop_ordering[0],bsg_2->m_sblocks[iblock]->m_loop_ordering[1],bsg_2->m_sblocks[iblock]->m_loop_ordering[2]}};

        unsigned nNodes = (1+sgridsizes.node_max[0])*(1+sgridsizes.node_max[1])*(1+sgridsizes.node_max[2]);
        for(unsigned iNode=0;iNode<nNodes;iNode++) {

            unsigned indx[3];
            const int L0 = loop_orderings[0], L1 =
                    loop_orderings[1], L2 = loop_orderings[2];
            const unsigned int sgsizes[3] = { 1 + sgridsizes.node_max[L0]
                                                                      - sgridsizes.node_min[L0], 1
                                                                      + sgridsizes.node_max[L1]
                                                                                            - sgridsizes.node_min[L1], 1
                                                                                            + sgridsizes.node_max[L2]
                                                                                                                  - sgridsizes.node_min[L2] };
            indx[L2] = sgridsizes.node_min[L2]
                                           + (iNode / (sgsizes[L0] * sgsizes[L1]));
            indx[L1] = sgridsizes.node_min[L1]
                                           + ((iNode / sgsizes[L0]) % sgsizes[L1]);
            indx[L0] = sgridsizes.node_min[L0] + (iNode % sgsizes[L0]);

            unsigned i=indx[0];
            unsigned j=indx[1];
            unsigned k=indx[2];

            for (int d=0; d<3; d++) {
                data[d]=host_coord_field_iblock(i,j,k,d);
            }
            const double x = data[0], y = data[1], z = data[2];
            if (x > tol && x< 1-tol &&
                    y > tol && y< 1-tol &&
                    z > tol && z< 1-tol) {

                int org_ijk[3];

                org_ijk[0] = (int)((data[0]*noElems) + 0.5);
                org_ijk[1] = (int)((data[1]*noElems) + 0.5);
                org_ijk[2] = (int)((data[2]*noElems) + 0.5);

                for (int d=0; d<3; d++)
                    data[d] = host_coord_field(org_ijk[0],org_ijk[1],org_ijk[2],d);

            }//tolerance
            for (int d=0; d<3; d++) {
                host_coord_field_iblock(i,j,k,d)=data[d];
            }
        }//foreach node
        Kokkos::deep_copy(coord_field_iblock, host_coord_field_iblock); //sync back block field
    }
    eMesh_2.copy_field(coord_field_N_2, coord_field_2);
}

TEST(heavy_perceptMeshSmoother_mbvsb, total_metric_and_gradient)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD ;
    MPI_Barrier( MPI_COMM_WORLD );

    const unsigned p_size = stk::parallel_machine_size( pm );
    if (p_size > 1) return;

    Double mtot_1block = 0;
    unsigned n_invalid_1block = 0;


    PerceptMesh multiBlock(3u);
    std::string fileName = "cube_4blocks.cgns";

    bool doRef = true;
    unsigned noRefs = 4;//mesh resolution doesn't seem to have an effect on the dot product of the gradient fields for this mesh or for the metrics
    PerceptMesh singleBlock(3u);
    bool add_all_fields = true;
    unsigned noElems = 5;
    build_perturbed_cube_multiblock(singleBlock, noElems, multiBlock, fileName, add_all_fields, doRef, noRefs);
    if(doRef)
        for(unsigned iRef=0;iRef<noRefs;iRef++)
            noElems *= 2;
    SGridGenericAlgorithm_total_element_metric< HexMeshSmootherMetric >  ga_tem_1block(&singleBlock, mtot_1block, n_invalid_1block);

    ga_tem_1block.m_metric.m_untangling = true;
    ga_tem_1block.m_metric.m_use_ref_mesh = true;
    ga_tem_1block.run(0);
    mtot_1block = ga_tem_1block.mtot;
    n_invalid_1block = ga_tem_1block.n_invalid;

    std::cout << std::setprecision(15) << "mtot and n_invalid from single block case UNTANGLE metric " <<  mtot_1block << "  " << n_invalid_1block <<std::endl;

    mtot_1block = ga_tem_1block.mtot = 0.0 ;
    n_invalid_1block = ga_tem_1block.n_invalid = 0.0;
    ga_tem_1block.m_metric.m_untangling = false;
    ga_tem_1block.run(0);
    ga_tem_1block.run(0);
    mtot_1block = ga_tem_1block.mtot;
    n_invalid_1block = ga_tem_1block.n_invalid;

    std::cout << std::setprecision(15) << "mtot and n_invalid from single block case SMOOTHING metric " <<  mtot_1block << "  " << n_invalid_1block <<std::endl;

    ///////////////////////////////multiblock metric
    Double mtot_multiBlock = 0;
    unsigned n_invalid_multiBlock = 0;
    SGridGenericAlgorithm_total_element_metric< HexMeshSmootherMetric >  ga_tem_multiBlock(&multiBlock, mtot_multiBlock, n_invalid_multiBlock);

    ga_tem_multiBlock.m_metric.m_untangling = true;
    ga_tem_multiBlock.m_metric.m_use_ref_mesh = true;
    for(unsigned iBlock=0;iBlock<multiBlock.get_block_structured_grid()->m_sblocks.size();iBlock++) {
        ga_tem_multiBlock.run(iBlock);

        mtot_multiBlock += ga_tem_multiBlock.mtot;
        n_invalid_multiBlock += ga_tem_multiBlock.n_invalid;
    }
    std::cout << std::setprecision(15) << "mtot and n_invalid from multi block case UNTANGLE metric " <<  mtot_multiBlock << "  " << n_invalid_multiBlock <<std::endl;
    mtot_multiBlock = ga_tem_multiBlock.mtot = 0;
    n_invalid_multiBlock = ga_tem_multiBlock.n_invalid = 0;
    ga_tem_multiBlock.m_metric.m_untangling = false;
    for(unsigned iBlock=0;iBlock<multiBlock.get_block_structured_grid()->m_sblocks.size();iBlock++) {
        ga_tem_multiBlock.run(iBlock);

        mtot_multiBlock += ga_tem_multiBlock.mtot;
        n_invalid_multiBlock += ga_tem_multiBlock.n_invalid;
    }
    std::cout << std::setprecision(15) << "mtot and n_invalid from multi block case SMOOTHING metric " <<  mtot_multiBlock << "  " << n_invalid_multiBlock <<std::endl;

    ///////////////////////////////multiblock metric

    EXPECT_NEAR(mtot_multiBlock,mtot_1block,1e-10);
    EXPECT_EQ(n_invalid_multiBlock,n_invalid_1block);

    /////////////////single block gradients
    StructuredGrid::MTField::Array4D grad_fld_data = *singleBlock.get_block_structured_grid()->m_fields["cg_g"].get()->m_block_fields[0];
    StructuredGrid::MTField * grad_field           =  singleBlock.get_block_structured_grid()->m_fields["cg_g"].get();

    sGrid_GenericAlgorithm_get_gradient_1 gg_1_sgrid (&singleBlock, NULL);

    // grad for smoothing
    gg_1_sgrid.m_metric.m_untangling=false;
    gg_1_sgrid.run();
    double g_1block_smooth = std::sqrt(singleBlock.nodal_field_dot(grad_field, grad_field));
    std::cout << std::setprecision(12) << "           single block grad smoothing = " << g_1block_smooth << std::endl;

    // grad for untangling
    gg_1_sgrid.m_metric.m_untangling=true;
    zero_a4d za_1(grad_fld_data);
    za_1.zero_out_fields();

    gg_1_sgrid.run();
    double g_1block_untangle = std::sqrt(singleBlock.nodal_field_dot(grad_field, grad_field));
    std::cout << std::setprecision(12) << "           single block grad untangling = " << g_1block_untangle << std::endl;

    za_1.zero_out_fields();
    /////////////////single block gradients

    /////////////////MULTI block gradients
    StructuredGrid::MTField * grad_field_2           =  multiBlock.get_block_structured_grid()->m_fields["cg_g"].get();
    sGrid_GenericAlgorithm_get_gradient_1 gg_2_sgrid (&multiBlock, NULL);
    // grad for smoothing
    gg_2_sgrid.m_metric.m_untangling=false;
    gg_2_sgrid.run();
    double g_multi_smooth = std::sqrt(multiBlock.nodal_field_dot(grad_field_2, grad_field_2));
    std::cout << std::setprecision(12) << "           multiblock grad smoothing = " << g_multi_smooth << std::endl;

    // grad for untangling
    gg_2_sgrid.m_metric.m_untangling=true;

    for(unsigned iBlock=0;iBlock<multiBlock.get_block_structured_grid()->m_sblocks.size();iBlock++){//zero out all blocks
        StructuredGrid::MTField::Array4D grad_fld_data_2 = *multiBlock.get_block_structured_grid()->m_fields["cg_g"].get()->m_block_fields[iBlock];
        zero_a4d za_2(grad_fld_data_2);
        za_2.zero_out_fields();
    }

    gg_2_sgrid.run();
    double g_multi_untangle = std::sqrt(multiBlock.nodal_field_dot(grad_field_2, grad_field_2));
    std::cout << std::setprecision(12) << "           multiblock grad untangling = " << g_multi_untangle << std::endl;

    for(unsigned iBlock=0;iBlock<multiBlock.get_block_structured_grid()->m_sblocks.size();iBlock++){//zero out all blocks
        StructuredGrid::MTField::Array4D grad_fld_data_2 = *multiBlock.get_block_structured_grid()->m_fields["cg_g"].get()->m_block_fields[iBlock];
        zero_a4d za_2(grad_fld_data_2);
        za_2.zero_out_fields();
    }
    /////////////////MULTI block gradients
    EXPECT_NEAR(g_multi_smooth,g_1block_smooth,1e-8);
    EXPECT_NEAR(g_multi_untangle,g_1block_untangle,1e-10);

    if(add_all_fields)
    {
        unsigned spatialDim = 3;
        is_field_data_equivalent("cg_g", noElems, spatialDim, multiBlock.get_block_structured_grid(), singleBlock.get_block_structured_grid() );
    }
}

TEST(heavy_perceptMeshSmoother_mbvsb, dot_product)
{
    PerceptMesh multiBlock(3u);
    std::string fileName = "cube_4blocks.cgns";
    std::string dbtype("cgns_structured");
    multiBlock.open(fileName, dbtype);
    multiBlock.readBulkData();
    std::shared_ptr<BlockStructuredGrid> bsg_mb = multiBlock.get_block_structured_grid();

    PerceptMesh singleBlock(3u);
    std::cout << "about to build unperturbed cube\n";
    unsigned nele = 5;
    const unsigned nxyz = nele+1;
    std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
    std::shared_ptr<BlockStructuredGrid> bsg_sb = BlockStructuredGrid::fixture_1(singleBlock.parallel(), sizes);
    singleBlock.set_block_structured_grid(bsg_sb);


    bsg_mb->register_field("coordinates", multiBlock.get_spatial_dim());

    StructuredGrid::MTField * multiblock_field           = bsg_mb->m_fields["coordinates"].get();
    StructuredGrid::MTField * singlelblock_field         = bsg_sb->m_fields["coordinates"].get();

    double mbDP2 = (multiBlock.nodal_field_dot(multiblock_field, multiblock_field));
    double sbDP2 = (singleBlock.nodal_field_dot(singlelblock_field, singlelblock_field));
    std::cout << "dot product on multiblock mesh " << mbDP2 << std::endl;
    std::cout << "dot product on single block mesh " << sbDP2 << std::endl;

    EXPECT_NEAR(mbDP2,sbDP2,1e-12);
}

TEST(heavy_perceptMeshSmoother_mbvsb, smooth)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD ;
    MPI_Barrier( MPI_COMM_WORLD );

    const unsigned p_size = stk::parallel_machine_size( pm );
    if (p_size > 1) return;

    int innerIter = 1001;
    const double smooth_tol = 1.e-3;

    PerceptMesh multiBlock(3u);
    std::string fileName = "cube_4blocks.cgns";

    PerceptMesh singleBlock(3u);
    build_perturbed_cube_multiblock(singleBlock, 5, multiBlock, fileName, true);
    singleBlock.setProperty("in_filename","singleBlock");

    std::cout << "Dot of mb cf " << multiBlock.nodal_field_dot("coordinates","coordinates") << std::endl;
    std::cout << "Dot of sb cf " << singleBlock.nodal_field_dot("coordinates","coordinates") << std::endl;

    {
        stk::diag::Timer root_timer(singleBlock.getProperty("in_filename"), rootTimerStructured());
        stk::diag::TimeBlock root_block(root_timer);
        SGridBoundarySelector boundarySelector (&singleBlock);
        ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&singleBlock, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);
        pmmpsi.m_do_animation = 0;
        pmmpsi.run();
    }
    multiBlock.setProperty("in_filename","multiBlock");
    {
        stk::diag::Timer root_timer(multiBlock.getProperty("in_filename"), rootTimerStructured());
        stk::diag::TimeBlock root_block(root_timer);
        SGridBoundarySelector boundarySelector (&multiBlock);
        ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&multiBlock, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);
        pmmpsi.m_do_animation = 0;
        pmmpsi.run();
    }

    std::cout << "Dot of mb cf " << multiBlock.nodal_field_dot("coordinates","coordinates") << std::endl;
    std::cout << "Dot of sb cf " << singleBlock.nodal_field_dot("coordinates","coordinates") << std::endl;
}

TEST(performance, refine_smooth_cgns)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD ;
    MPI_Barrier( MPI_COMM_WORLD );

    const unsigned p_size = stk::parallel_machine_size( pm );
    if (p_size > 1) return;

    int innerIter = 1001;
    const double smooth_tol = 1.e-3;

    PerceptMesh multiBlock(3u);
    std::string fileName = "cube_4blocks.cgns";

    PerceptMesh singleBlock(3u);
    bool refine = true;

    unsigned nRefs = 2;
    build_perturbed_cube_multiblock(singleBlock, 5, multiBlock, fileName, true,refine,nRefs);
    singleBlock.setProperty("in_filename","singleBlock");

    double  multiBlock_dot =  multiBlock.nodal_field_dot("coordinates","coordinates");
    double singleBlock_dot = singleBlock.nodal_field_dot("coordinates","coordinates");
    EXPECT_NEAR(multiBlock_dot, singleBlock_dot, 1e-12);

    {
        stk::diag::Timer root_timer(singleBlock.getProperty("in_filename"), rootTimerStructured());
        stk::diag::TimeBlock root_block(root_timer);
        SGridBoundarySelector boundarySelector (&singleBlock);
        ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&singleBlock, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);
        pmmpsi.m_do_animation = 0;
        pmmpsi.run();
    }
    multiBlock.setProperty("in_filename","multiBlock");
    {
        stk::diag::Timer root_timer(multiBlock.getProperty("in_filename"), rootTimerStructured());
        stk::diag::TimeBlock root_block(root_timer);
        SGridBoundarySelector boundarySelector (&multiBlock);
        ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&multiBlock, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);
        pmmpsi.m_do_animation = 0;
        pmmpsi.run();
    }

     multiBlock_dot =  multiBlock.nodal_field_dot("coordinates","coordinates");
    singleBlock_dot = singleBlock.nodal_field_dot("coordinates","coordinates");
    EXPECT_NEAR(multiBlock_dot, singleBlock_dot, 2e-1);
}

void grade_mesh_along_x_axis_and_warp_y(PerceptMesh * eMesh, std::string out_name = "graded_cube_mesh", bool dump = false)
{
    //Hard coded to grade a cube mesh on intervals [0,1] [0,1] [0,1]
    if(!eMesh)
    {
        std::cout << "invalid mesh pointer, aborting grading" << std::endl;
        return;
    }
    if(!eMesh->get_block_structured_grid())
    {
        std::cout << "invalid bsg pointer, aborting grading" << std::endl;
    }

    std::shared_ptr<BlockStructuredGrid> bsg = eMesh->get_block_structured_grid();

    for(unsigned iBlock = 0; iBlock < bsg->m_sblocks.size(); iBlock++)
    {
        Array4D field_data = *bsg->m_fields["coordinates_NM1"].get()->m_block_fields[iBlock];
        Array4D::HostMirror field_data_mir = Kokkos::create_mirror_view(field_data);
        Kokkos::deep_copy(field_data_mir,field_data);

        SGridSizes sizes = bsg->m_sblocks[iBlock]->m_sizes;
        Kokkos::Array<unsigned, 3> lo;
        lo[0] = bsg->m_sblocks[iBlock]->m_loop_ordering[0];
        lo[1] = bsg->m_sblocks[iBlock]->m_loop_ordering[1];
        lo[2] = bsg->m_sblocks[iBlock]->m_loop_ordering[2];

        for(unsigned iNode = 0; iNode < sizes.node_size[0]*sizes.node_size[1]*sizes.node_size[2];iNode++ )
        {
            Kokkos::Array<unsigned, 3> ijk;
            sgrid_multi_dim_indices_from_index_node(sizes,lo,iNode,ijk);

            Double x = field_data_mir(ijk[0],ijk[1],ijk[2],0);
            Double y = field_data_mir(ijk[0],ijk[1],ijk[2],1);
            Double z = field_data_mir(ijk[0],ijk[1],ijk[2],2);
            Double y_mut=0.0;
            if(x>0.5) //make sure deformation isn't symmetric across block boundaries
                y_mut = 0.2;
            if(x<0.5)
                y_mut = -0.3;


            if(device_safe_abs_flt( Double(1.0-x) )<1e-14 || device_safe_abs_flt( Double(0.0-x) )<1e-14) //only deform nodes between x=0 and x=1
                continue;
            x = x/2.0;
            z = z*z;
            y += y_mut;
            field_data_mir(ijk[0],ijk[1],ijk[2],0) = x;
            field_data_mir(ijk[0],ijk[1],ijk[2],1) = y;
            field_data_mir(ijk[0],ijk[1],ijk[2],2) = z;
        }
        Kokkos::deep_copy(field_data,field_data_mir);
    }
    eMesh->copy_field("coordinates","coordinates_NM1");

    if(dump)
        bsg->dump_vtk(out_name);
}

TEST(heavy_perceptMeshSmoother_mbvsb, edge_len_avg)
{   //verifies that both multi block and single block ela achieve the same results
    PerceptMesh multiBlock(3u);
    PerceptMesh singleBlock(3u);

    bool do_print = true;

    std::string fileName = "cube_4blocks.cgns";
    std::string dbtype("cgns_structured");
    multiBlock.open(fileName, dbtype);
    multiBlock.readBulkData();
    std::shared_ptr<BlockStructuredGrid> bsg_mb = multiBlock.get_block_structured_grid();

    unsigned nele = 5;
    const unsigned nxyz = nele+1;
    std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
    std::shared_ptr<BlockStructuredGrid> bsg_sb = BlockStructuredGrid::fixture_1(singleBlock.parallel(), sizes);
    singleBlock.set_block_structured_grid(bsg_sb);

    bsg_mb->register_field("coordinates_saved", multiBlock.get_spatial_dim());
    bsg_mb->register_field("coordinates_NM1", multiBlock.get_spatial_dim());
    bsg_mb->register_field("coordinates", multiBlock.get_spatial_dim());
    bsg_mb->register_field("cg_edge_length", 1);
    bsg_mb->register_field("num_adj_elems", 1);

    multiBlock.copy_field("coordinates_NM1","coordinates");
    multiBlock.copy_field("coordinates_saved","coordinates");
    multiBlock.nodal_field_set_value("cg_edge_length",0.0);

    bsg_sb->register_field("coordinates_saved", multiBlock.get_spatial_dim());
    bsg_sb->register_field("coordinates_NM1", multiBlock.get_spatial_dim());
    bsg_sb->register_field("cg_edge_length", 1);
    bsg_sb->register_field("num_adj_elems", 1);

    singleBlock.copy_field("coordinates_NM1","coordinates");
    singleBlock.copy_field("coordinates_saved","coordinates");
    singleBlock.nodal_field_set_value("cg_edge_length",0.0);

    sGrid_GenericAlgorithm_get_edge_lengths sb_ela(&singleBlock);
    sGrid_GenericAlgorithm_get_edge_lengths mb_ela(&multiBlock);

    sb_ela.calc_edge_lengths(do_print);

    mb_ela.calc_edge_lengths(do_print);

    double uniform_coord_field_dp_mb = (multiBlock.nodal_field_dot("cg_edge_length", "cg_edge_length"));
    double uniform_coord_field_dp_sb = (singleBlock.nodal_field_dot("cg_edge_length", "cg_edge_length"));

    std::cout << "dot product on ela_field multiblock mesh after fixup " << uniform_coord_field_dp_mb << std::endl;
    std::cout << "dot product on ela_field single block mesh after fixup " << uniform_coord_field_dp_sb << std::endl;

    EXPECT_NEAR(uniform_coord_field_dp_mb,uniform_coord_field_dp_sb,1e-12);

    grade_mesh_along_x_axis_and_warp_y(&multiBlock,"multi_grade");
    grade_mesh_along_x_axis_and_warp_y(&singleBlock,"single_grade");

    singleBlock.nodal_field_set_value("cg_edge_length",0.0);
    multiBlock.nodal_field_set_value("cg_edge_length",0.0);

    sb_ela.adj_elem_field_filled = false;
    sb_ela.adj_elem_field_fixed = false;
    sb_ela.calc_edge_lengths(do_print);

    mb_ela.adj_elem_field_filled = false;
    mb_ela.adj_elem_field_fixed = false;
    mb_ela.calc_edge_lengths(do_print);

    double non_uniform_coord_field_dp_mb = (multiBlock.nodal_field_dot("cg_edge_length", "cg_edge_length"));
    double non_uniform_coord_field_dp_sb = (singleBlock.nodal_field_dot("cg_edge_length", "cg_edge_length"));

    std::cout << "dot product on ela_field multiblock mesh after fixup " << non_uniform_coord_field_dp_mb << std::endl;
    std::cout << "dot product on ela_field singleblock mesh after fixup " << non_uniform_coord_field_dp_sb << std::endl;

    EXPECT_NEAR(non_uniform_coord_field_dp_mb,non_uniform_coord_field_dp_sb,1e-12);


    StructuredGrid::MTField *coord_field_mb = bsg_mb->m_fields["coordinates_saved"].get();
    Array4D sb_ela_fld =  *(bsg_sb->m_fields["cg_edge_length"].get()->m_block_fields[0]);

    double data[3];
    for (unsigned iblock=0; iblock < bsg_mb->m_sblocks.size(); ++iblock)
    {
        Array4D coord_field_iblock = *(coord_field_mb->m_block_fields[iblock]);
        Array4D::HostMirror host_coord_field_iblock = Kokkos::create_mirror_view(coord_field_iblock);
        Kokkos::deep_copy(host_coord_field_iblock,coord_field_iblock);
        const SGridSizes sgridsizes = bsg_mb->m_sblocks[iblock]->m_sizes;
        const std::array<unsigned int, 3> loop_orderings = {{bsg_mb->m_sblocks[iblock]->m_loop_ordering[0],bsg_mb->m_sblocks[iblock]->m_loop_ordering[1],bsg_mb->m_sblocks[iblock]->m_loop_ordering[2]}};

        Array4D mb_ela_fld =  *(bsg_mb->m_fields["cg_edge_length"].get()->m_block_fields[iblock]);

        unsigned nNodes = (1+sgridsizes.node_max[0])*(1+sgridsizes.node_max[1])*(1+sgridsizes.node_max[2]);
        for(unsigned iNode=0;iNode<nNodes;iNode++) {

            unsigned indx[3];
            const int L0 = loop_orderings[0], L1 =
                    loop_orderings[1], L2 = loop_orderings[2];
            const unsigned int sgsizes[3] = { 1 + sgridsizes.node_max[L0]
                                                                      - sgridsizes.node_min[L0], 1
                                                                      + sgridsizes.node_max[L1]
                                                                                            - sgridsizes.node_min[L1], 1
                                                                                            + sgridsizes.node_max[L2]
                                                                                                                  - sgridsizes.node_min[L2] };
            indx[L2] = sgridsizes.node_min[L2]
                                           + (iNode / (sgsizes[L0] * sgsizes[L1]));
            indx[L1] = sgridsizes.node_min[L1]
                                           + ((iNode / sgsizes[L0]) % sgsizes[L1]);
            indx[L0] = sgridsizes.node_min[L0] + (iNode % sgsizes[L0]);

            unsigned i=indx[0];
            unsigned j=indx[1];
            unsigned k=indx[2];

            for (int d=0; d<3; d++)
                data[d]=host_coord_field_iblock(i,j,k,d);


            int org_ijk[3];

            org_ijk[0] = (int)((data[0]*nele) + 0.5);
            org_ijk[1] = (int)((data[1]*nele) + 0.5);
            org_ijk[2] = (int)((data[2]*nele) + 0.5);

            EXPECT_NEAR(mb_ela_fld(i,j,k,0),sb_ela_fld(org_ijk[0],org_ijk[1],org_ijk[2],0),1e-12);
        }//foreach node
    }
}

TEST(heavy_perceptMeshSmoother_mbvsb, edge_len_avg_LARGE)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD ;
    MPI_Barrier( MPI_COMM_WORLD );

    const unsigned p_size = stk::parallel_machine_size( pm );
    if (p_size > 1) return;

    PerceptMesh multiBlock(3u);
    std::string fileName = "cube_4blocks.cgns";

    bool doRef = true;
    unsigned noRefs = 3;
    bool add_all_fields = true;
    PerceptMesh singleBlock(3u);
    unsigned nele = 5;
    build_perturbed_cube_multiblock(singleBlock, nele, multiBlock, fileName, add_all_fields, doRef, noRefs);
    bool do_print = true;
    for(unsigned iRef=0;iRef<noRefs;iRef++)
        nele *=2;

    sGrid_GenericAlgorithm_get_edge_lengths sb_ela(&singleBlock);
    sGrid_GenericAlgorithm_get_edge_lengths mb_ela(&multiBlock);

    sb_ela.calc_edge_lengths(do_print);
    mb_ela.calc_edge_lengths(do_print);

    std::shared_ptr<BlockStructuredGrid> bsg_mb = multiBlock.get_block_structured_grid();
    std::shared_ptr<BlockStructuredGrid> bsg_sb = singleBlock.get_block_structured_grid();

    StructuredGrid::MTField *coord_field_mb = bsg_mb->m_fields["coordinates_saved"].get();
    Array4D sb_ela_fld =  *(bsg_sb->m_fields["cg_edge_length"].get()->m_block_fields[0]);

    double data[3];
    for (unsigned iblock=0; iblock < bsg_mb->m_sblocks.size(); ++iblock)
    {
        Array4D coord_field_iblock = *(coord_field_mb->m_block_fields[iblock]);
        Array4D::HostMirror host_coord_field_iblock = Kokkos::create_mirror_view(coord_field_iblock);
        Kokkos::deep_copy(host_coord_field_iblock,coord_field_iblock);
        const SGridSizes sgridsizes = bsg_mb->m_sblocks[iblock]->m_sizes;
        const std::array<unsigned int, 3> loop_orderings = {{bsg_mb->m_sblocks[iblock]->m_loop_ordering[0],bsg_mb->m_sblocks[iblock]->m_loop_ordering[1],bsg_mb->m_sblocks[iblock]->m_loop_ordering[2]}};

        Array4D mb_ela_fld =  *(bsg_mb->m_fields["cg_edge_length"].get()->m_block_fields[iblock]);

        unsigned nNodes = (1+sgridsizes.node_max[0])*(1+sgridsizes.node_max[1])*(1+sgridsizes.node_max[2]);
        for(unsigned iNode=0;iNode<nNodes;iNode++) {

            unsigned indx[3];
            const int L0 = loop_orderings[0], L1 =
                    loop_orderings[1], L2 = loop_orderings[2];
            const unsigned int sgsizes[3] = { 1 + sgridsizes.node_max[L0]
                                                                      - sgridsizes.node_min[L0], 1
                                                                      + sgridsizes.node_max[L1]
                                                                                            - sgridsizes.node_min[L1], 1
                                                                                            + sgridsizes.node_max[L2]
                                                                                                                  - sgridsizes.node_min[L2] };
            indx[L2] = sgridsizes.node_min[L2]
                                           + (iNode / (sgsizes[L0] * sgsizes[L1]));
            indx[L1] = sgridsizes.node_min[L1]
                                           + ((iNode / sgsizes[L0]) % sgsizes[L1]);
            indx[L0] = sgridsizes.node_min[L0] + (iNode % sgsizes[L0]);

            unsigned i=indx[0];
            unsigned j=indx[1];
            unsigned k=indx[2];

            for (int d=0; d<3; d++)
                data[d]=host_coord_field_iblock(i,j,k,d);

            int org_ijk[3];

            org_ijk[0] = (int)((data[0]*nele) + 0.5);
            org_ijk[1] = (int)((data[1]*nele) + 0.5);
            org_ijk[2] = (int)((data[2]*nele) + 0.5);

            EXPECT_NEAR(mb_ela_fld(i,j,k,0),sb_ela_fld(org_ijk[0],org_ijk[1],org_ijk[2],0),1e-12);
        }//foreach node
    }//foreach block
}

TEST(heavy_perceptMeshSmoother_mbvsb, bsg_field_additions)
{
    PerceptMesh multiBlock(3u);
    PerceptMesh singleBlock(3u);

    std::string fileName = "cube_4blocks.cgns";
    std::string dbtype("cgns_structured");
    multiBlock.open(fileName, dbtype);
    multiBlock.readBulkData();
    std::shared_ptr<BlockStructuredGrid> bsg_mb = multiBlock.get_block_structured_grid();

    unsigned nele = 5;
    const unsigned nxyz = nele+1;
    std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
    std::shared_ptr<BlockStructuredGrid> bsg_sb = BlockStructuredGrid::fixture_1(singleBlock.parallel(), sizes);
    singleBlock.set_block_structured_grid(bsg_sb);

    bsg_sb->register_field("alpha_field",singleBlock.get_spatial_dim());
    bsg_mb->register_field("alpha_field",multiBlock.get_spatial_dim());

    bsg_sb->register_field("beta_field",singleBlock.get_spatial_dim());
    bsg_mb->register_field("beta_field",multiBlock.get_spatial_dim());

    bsg_sb->register_field("gamma_field",singleBlock.get_spatial_dim());
    bsg_mb->register_field("gamma_field",multiBlock.get_spatial_dim());

    StructuredGrid::MTField* alpha_mb = bsg_mb->m_fields["alpha_field"].get();
    StructuredGrid::MTField* beta_mb = bsg_mb->m_fields["beta_field"].get();
    StructuredGrid::MTField* gamma_mb = bsg_mb->m_fields["gamma_field"].get();

    StructuredGrid::MTField* alpha_sb = bsg_sb->m_fields["alpha_field"].get();
    StructuredGrid::MTField* beta_sb = bsg_sb->m_fields["beta_field"].get();
    StructuredGrid::MTField* gamma_sb = bsg_sb->m_fields["gamma_field"].get();

    singleBlock.nodal_field_set_value("beta_field",0.0);
    multiBlock.nodal_field_set_value("beta_field",0.0);

    singleBlock.nodal_field_set_value("gamma_field",0.0);
    multiBlock.nodal_field_set_value("gamma_field",0.0);

    singleBlock.nodal_field_set_value("alpha_field",1.0);
    multiBlock.nodal_field_set_value("alpha_field",1.0);

    bsg_sb->nodal_field_axpby(2.0,alpha_sb,0.0,beta_sb);
    bsg_mb->nodal_field_axpby(2.0,alpha_mb,0.0,beta_mb);

    double beta_fld_sqr_mb = (multiBlock.nodal_field_dot("beta_field", "beta_field"));
    double beta_fld_sqr_sb = (singleBlock.nodal_field_dot("beta_field", "beta_field"));

    EXPECT_NEAR(beta_fld_sqr_mb,beta_fld_sqr_sb,1e-12);

    bsg_sb->nodal_field_axpbypgz(2.0,alpha_sb,1.0,beta_sb,0.0,gamma_sb);
    bsg_mb->nodal_field_axpbypgz(2.0,alpha_mb,1.0,beta_mb,0.0,gamma_mb);

    double gamma_fld_sqr_mb = (multiBlock.nodal_field_dot("gamma_field", "gamma_field"));
    double gamma_fld_sqr_sb = (singleBlock.nodal_field_dot("gamma_field", "gamma_field"));

    EXPECT_NEAR(gamma_fld_sqr_mb,gamma_fld_sqr_sb,1e-12);
}

TEST(heavy_perceptMeshSmoother_mbvsb, get_alpha_0)
{
    PerceptMesh multiBlock(3u);
    PerceptMesh singleBlock(3u);

    std::string fileName = "cube_4blocks.cgns";
    std::string dbtype("cgns_structured");
    multiBlock.open(fileName, dbtype);
    multiBlock.readBulkData();
    std::shared_ptr<BlockStructuredGrid> bsg_mb = multiBlock.get_block_structured_grid();

    unsigned nele = 5;
    const unsigned nxyz = nele+1;
    std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
    std::shared_ptr<BlockStructuredGrid> bsg_sb = BlockStructuredGrid::fixture_1(singleBlock.parallel(), sizes);
    singleBlock.set_block_structured_grid(bsg_sb);

    bsg_sb->register_field("cg_s",singleBlock.get_spatial_dim());
    bsg_sb->register_field("cg_edge_length",1);

    bsg_mb->register_field("cg_s",multiBlock.get_spatial_dim());
    bsg_mb->register_field("cg_edge_length",1);

    singleBlock.nodal_field_set_value("cg_s",1.0);
    singleBlock.nodal_field_set_value("cg_edge_length",1.0);

    multiBlock.nodal_field_set_value("cg_s",1.0);
    multiBlock.nodal_field_set_value("cg_edge_length",1.0);

    std::cout << "mb cg sqr " << (multiBlock.nodal_field_dot("cg_s", "cg_s"))  <<std::endl;
    std::cout << "mb ela sqr "<<(multiBlock.nodal_field_dot("cg_edge_length", "cg_edge_length")) <<std::endl;

    std::cout << "sb cg sqr " << (singleBlock.nodal_field_dot("cg_s", "cg_s"))  <<std::endl;
    std::cout << "sb ela sqr "<<(singleBlock.nodal_field_dot("cg_edge_length", "cg_edge_length")) <<std::endl;

    SGrid_GenericAlgorithm_get_alpha_0 ga_sb(&singleBlock);
    SGrid_GenericAlgorithm_get_alpha_0 ga_mb(&multiBlock);

    ga_sb.calc_alpha();

    ga_mb.calc_alpha();

    EXPECT_NEAR(ga_mb.alpha,ga_sb.alpha,1e-12);

    bool do_print = true;
    //calc alpha on deformed mesh
    bsg_mb->register_field("coordinates_NM1", multiBlock.get_spatial_dim());
    bsg_mb->register_field("coordinates", multiBlock.get_spatial_dim());
    bsg_mb->register_field("num_adj_elems", 1);

    multiBlock.copy_field("coordinates_NM1","coordinates");
    multiBlock.nodal_field_set_value("cg_edge_length",0.0);
    multiBlock.nodal_field_set_value("num_adj_elems",0.0);

    bsg_sb->register_field("coordinates_NM1", multiBlock.get_spatial_dim());
    bsg_sb->register_field("num_adj_elems", 1);

    singleBlock.copy_field("coordinates_NM1","coordinates");
    singleBlock.nodal_field_set_value("cg_edge_length",0.0);
    singleBlock.nodal_field_set_value("num_adj_elems",0.0);

    grade_mesh_along_x_axis_and_warp_y(&multiBlock,"multi_grade");
    grade_mesh_along_x_axis_and_warp_y(&singleBlock,"single_grade");

    sGrid_GenericAlgorithm_get_edge_lengths sb_ela(&singleBlock);
    sGrid_GenericAlgorithm_get_edge_lengths mb_ela(&multiBlock);

    sb_ela.calc_edge_lengths(do_print);
    mb_ela.calc_edge_lengths(do_print);

    //now that you have a real edge length field, test get_alpha against it
    ga_sb.reset_alpha();
    ga_mb.reset_alpha();

    ga_sb.calc_alpha();

    ga_mb.calc_alpha();

    std::cout << "alpha of sb case is ... " << ga_sb.alpha << std::endl;
    std::cout << "alpha of mb case is ... " << ga_mb.alpha << std::endl;
}

TEST(heavy_perceptMeshSmoother_mbvsb, update_coords)
{
    int innerIter = 1001;
    const double smooth_tol = 1.e-3;

    PerceptMesh multiBlock(3u);
    std::string fileName = "cube_4blocks.cgns";
    std::string dbtype("cgns_structured");
    multiBlock.open(fileName, dbtype);
    multiBlock.readBulkData();

    multiBlock.add_coordinate_state_fields();
    std::shared_ptr<BlockStructuredGrid> bsg_mb = multiBlock.get_block_structured_grid();
    bsg_mb->register_field("coordinates", multiBlock.get_spatial_dim());



    PerceptMesh singleBlock(3u);
    unsigned nele = 5;
    const unsigned nxyz = nele+1;
    std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
    std::shared_ptr<BlockStructuredGrid> bsg_sb = BlockStructuredGrid::fixture_1(singleBlock.parallel(), sizes);
    singleBlock.set_block_structured_grid(bsg_sb);
    singleBlock.add_coordinate_state_fields();

    singleBlock.setProperty("in_filename","singleBlock");
    {
        stk::diag::Timer root_timer(singleBlock.getProperty("in_filename"), rootTimerStructured());
        stk::diag::TimeBlock root_block(root_timer);
        SGridBoundarySelector boundarySelector (&singleBlock);
        ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&singleBlock, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);

        std::cout << "dot of sb coords before update " << singleBlock.nodal_field_dot("coordinates","coordinates") << std::endl;
        pmmpsi.get_edge_lengths(&singleBlock);

        singleBlock.nodal_field_set_value("cg_s",0.5);
        pmmpsi.m_coord_updater.alpha = 0.5;
        pmmpsi.m_coord_updater.run(true);

        std::cout << "dot of sb coords after update " << singleBlock.nodal_field_dot("coordinates","coordinates") << std::endl;
        std::cout << "dmax and relative " << pmmpsi.m_coord_updater.dmax << "  ---  " << pmmpsi.m_coord_updater.dmax_relative <<std::endl;
    }

    multiBlock.setProperty("in_filename","multiBlock ");
    {
        stk::diag::Timer root_timer(multiBlock.getProperty("in_filename"), rootTimerStructured());
        stk::diag::TimeBlock root_block(root_timer);
        SGridBoundarySelector boundarySelector (&multiBlock);
        ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&multiBlock, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);

        std::cout << "dot of mb coords before update " << multiBlock.nodal_field_dot("coordinates","coordinates") << std::endl;
        pmmpsi.get_edge_lengths(&multiBlock);

        multiBlock.nodal_field_set_value("cg_s",0.5);
        pmmpsi.m_coord_updater.alpha = 0.5;
        pmmpsi.m_coord_updater.run(true);

        std::cout << "dot of mb coords after update " << multiBlock.nodal_field_dot("coordinates","coordinates") << std::endl;
        std::cout << "dmax and relative " << pmmpsi.m_coord_updater.dmax << "  ---  " << pmmpsi.m_coord_updater.dmax_relative <<std::endl;
    }
}
}//percept
#endif
#endif
