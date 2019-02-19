/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "RFunit_OneD_KL_Example.h"

#include "percept/rfgen/RFGen_KLSolver.h"

// Trilinos includes
#include "AnasaziTypes.hpp"

namespace RFGen
{

using shards::Array;
using shards::NaturalOrder;

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

TEST(UQ, KLSolver_MatFree)
{
  /* %TRACE[OFF]% */  /* %TRACE% */

  // build class for 1d KL problem using exponential covariance function
  const int num_intervals = 10;
  std::vector<double> eigenVal_mem, eigenVec_mem;
  Teuchos::RCP<API_KLSolver> oneD_KLSolver = 
    Teuchos::rcp(new OneD_KLSolver(num_intervals, eigenVal_mem, eigenVec_mem));

  Teuchos::ParameterList MyPL;
  setupPL(MyPL);

  Teuchos::RCP<CovarianceFunction> covarFunc = 
    Teuchos::rcp(new ExpMagL2CovarianceFunction(1, std::vector<double>(1,0.3)));

  const bool useMatrixFree = true;

  Teuchos::RCP<KLSolver> klsolver = 
    buildKLSolver(oneD_KLSolver, covarFunc, MyPL, useMatrixFree);

  const int maxNev = 4;
  klsolver->solve(maxNev);

  // eigenvalues computed from Matlab svd(inv(B)*A)
  ASSERT_NEAR(eigenVal_mem[0], 0.442721666273175, 1e-6);
  ASSERT_NEAR(eigenVal_mem[1], 0.223675882509994, 1e-6);
  ASSERT_NEAR(eigenVal_mem[2], 0.113573433726618, 1e-6);
  ASSERT_NEAR(eigenVal_mem[3], 0.065712224680075, 1e-6);

  // checking on eigenvectors - serial only for now
  MPI_Comm mcomm = oneD_KLSolver->getParallelComm();
  int numProcs, myRank;
  EXPECT_TRUE(MPI_SUCCESS == MPI_Comm_size(mcomm, &numProcs));
  EXPECT_TRUE(MPI_SUCCESS == MPI_Comm_rank(mcomm, &myRank));

  if (numProcs==0) {
    ASSERT_NEAR(eigenVec_mem[0*num_intervals+3], 1.123003589165353, 1e-6);
    ASSERT_NEAR(eigenVec_mem[1*num_intervals+7], 1.140986382233420, 1e-6);
    ASSERT_NEAR(eigenVec_mem[2*num_intervals+2], 0.278054557198896, 1e-6);
    ASSERT_NEAR(eigenVec_mem[3*num_intervals+8], 0.485548660245600, 1e-6);
  }

  /*
   Eigenvectors - first 4 from Matlab "eig(A,B)"

   0.717531379713845  -1.170301548289961   1.337578715001846   1.340569886664697
   0.891335134486157  -1.277983237872343   1.066973064685960   0.485548660245600
   1.028376587404176  -1.140986382233418   0.278054557198896  -0.816795427655446
   1.123003589165353  -0.785539957109542  -0.645936949561721  -1.366647957327511
   1.171313335686556  -0.279696587626080  -1.256146002987107  -0.657444649833892
   1.171313335686555   0.279696587626079  -1.256146002987106   0.657444649833894
   1.123003589165352   0.785539957109541  -0.645936949561720   1.366647957327510
   1.028376587404175   1.140986382233420   0.278054557198897   0.816795427655443
   0.891335134486158   1.277983237872343   1.066973064685959  -0.485548660245600
   0.717531379713845   1.170301548289961   1.337578715001845  -1.340569886664697

   Evecs from Encore - note that last has sign change

   7.175316e-01 -1.170302e+00  1.337579e+00 -1.340564e+00
   8.913353e-01 -1.277983e+00  1.066973e+00 -4.855478e-01
   1.028377e+00 -1.140986e+00  2.780544e-01  8.167856e-01
   1.123004e+00 -7.855398e-01 -6.459371e-01  1.366649e+00
   1.171313e+00 -2.796965e-01 -1.256146e+00  6.574479e-01
   1.171313e+00  2.796960e-01 -1.256146e+00 -6.574380e-01
   1.123003e+00  7.855402e-01 -6.459371e-01 -1.366652e+00
   1.028377e+00  1.140986e+00  2.780544e-01 -8.168016e-01
   8.913352e-01  1.277983e+00  1.066973e+00  4.855478e-01
   7.175316e-01  1.170302e+00  1.337579e+00  1.340575e+00

  */
}


TEST(UQ, KLSolver_FullMat)
{
  /* %TRACE[OFF]% */  /* %TRACE% */

  // build class for 1d KL problem using exponential covariance function
  const int num_intervals = 10;
  std::vector<double> eigenVal_mem, eigenVec_mem;
  Teuchos::RCP<API_KLSolver> oneD_KLSolver = 
    Teuchos::rcp(new OneD_KLSolver(num_intervals, eigenVal_mem, eigenVec_mem));

  Teuchos::ParameterList MyPL;
  setupPL(MyPL);

  Teuchos::RCP<CovarianceFunction> covarFunc = 
    Teuchos::rcp(new ExpMagL2CovarianceFunction(1, std::vector<double>(1,0.3)));

  const bool useMatrixFree = false;

  Teuchos::RCP<KLSolver> klsolver = 
    buildKLSolver(oneD_KLSolver, covarFunc, MyPL, useMatrixFree);

  const int maxNev = 4;
  klsolver->solve(maxNev);

  // eigenvalues computed from Matlab svd(inv(B)*A)
  ASSERT_NEAR(eigenVal_mem[0], 0.442721666273175, 1e-6);
  ASSERT_NEAR(eigenVal_mem[1], 0.223675882509994, 1e-6);
  ASSERT_NEAR(eigenVal_mem[2], 0.113573433726618, 1e-6);
  ASSERT_NEAR(eigenVal_mem[3], 0.065712224680075, 1e-6);

  // checking on eigenvectors - serial only for now
  MPI_Comm mcomm = oneD_KLSolver->getParallelComm();
  int numProcs, myRank;
  EXPECT_TRUE(MPI_SUCCESS == MPI_Comm_size(mcomm, &numProcs));
  EXPECT_TRUE(MPI_SUCCESS == MPI_Comm_rank(mcomm, &myRank));

  if (numProcs==0) {
    ASSERT_NEAR(eigenVec_mem[0*num_intervals+3], 1.123003589165353, 1e-6);
    ASSERT_NEAR(eigenVec_mem[1*num_intervals+7], 1.140986382233420, 1e-6);
    ASSERT_NEAR(eigenVec_mem[2*num_intervals+2], 0.278054557198896, 1e-6);
    ASSERT_NEAR(eigenVec_mem[3*num_intervals+8], 0.485548660245600, 1e-6);
  }
}


} 
