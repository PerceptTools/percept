// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <stk_unit_tests/stk_mesh_fixtures/HexFixture.hpp>

#include <percept/uq/Percept_API_KLSolver.hpp>

#include <percept/rfgen/RFGen_CovarianceFunction.h>
#include <percept/rfgen/RFGen_KLSolver.h>

// Trilinos includes
#include "Teuchos_ParameterList.hpp"
#include "AnasaziTypes.hpp"

namespace percept {
  namespace unit_tests {

void setupPL(Teuchos::ParameterList &pl)
{
  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  //verbosity += Anasazi::Debug;
    
  // params for BlockDavidsonSolMgr
  const int blockSize = 2;
  const int numBlocks = 2;
  const int maxRestarts = 100;
  MagnitudeType tol = 1.0e-6;

  pl.set( "Verbosity", verbosity );
  pl.set( "Which", "LM" ); // which eigenvalues = largest magnitude?
  pl.set( "Block Size", blockSize );
  pl.set( "Num Blocks", numBlocks );
  pl.set( "Maximum Restarts", maxRestarts );
  pl.set( "Convergence Tolerance", tol );
  pl.set( "Use Locking", true );
  pl.set( "Locking Tolerance", tol/10 );
  
}

      TEST(rfgen, kl_solve)
      {
        // build 3D mesh with 10 cells in x direction
        stk::ParallelMachine pm = MPI_COMM_WORLD;
        const unsigned nx = 10 , ny = 1, nz = 1;

        stk::mesh::fixtures::HexFixture fixture( pm , nx , ny, nz);

        const int maxNev = 4;

        stk::mesh::Field<double,stk::mesh::SimpleArrayTag> & phi = fixture.m_meta.declare_field<stk::mesh::Field<double,stk::mesh::SimpleArrayTag> >(stk::topology::ELEMENT_RANK, "phi");
        stk::mesh::put_field_on_mesh(phi, fixture.m_meta.universal_part(), maxNev, nullptr);

        std::vector<double> lambda(maxNev);

        // build Percept KL solver
        Teuchos::RCP<RFGen::API_KLSolver> klSolver = 
          Teuchos::rcp(new Percept_API_KLSolver(fixture.m_bulk_data, phi, lambda));
 
        fixture.m_meta.commit();
        fixture.generate_mesh();

        Teuchos::RCP<RFGen::CovarianceFunction> covarFunc = 
          Teuchos::rcp(new RFGen::ExpMagL2CovarianceFunction(3, std::vector<double>(3,3.0)));

        const bool useMatrixFree = true;

        Teuchos::ParameterList MyPL;
        setupPL(MyPL);

        Teuchos::RCP<RFGen::KLSolver> klsolver = 
          buildKLSolver(klSolver, covarFunc, MyPL, useMatrixFree);
        
        klsolver->solve(maxNev);

        // eigenvalues computed from Matlab svd(inv(B)*A)
        ASSERT_NEAR(lambda[0]*0.1, 0.442721666273175, 1e-6);
        ASSERT_NEAR(lambda[1]*0.1, 0.223675882509994, 1e-6);
        ASSERT_NEAR(lambda[2]*0.1, 0.113573433726618, 1e-6);
        ASSERT_NEAR(lambda[3]*0.1, 0.065712224680075, 1e-6);
     }

  } // namespace unit_tests
} // namespace percept


