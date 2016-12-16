// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT) && ENABLE_SMOOTHER3

/*
 * TpetraLinearSolver.C
 *
 *  Created on: May 2, 2014
 *      Author: tcfishe
 *  Modified on: May 30, 2015
 *      Author: srkenno
 */


//#include <TpetraLinearSystem.h>
#include <percept/mesh/mod/smoother/TpetraLinearSolver.hpp>
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherNewton.hpp>
//#include <DiagWriter.h>

#include <BelosLinearProblem.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosSolverManager.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosTFQMRSolMgr.hpp>

#include <Ifpack2_Factory.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace percept {

  TpetraLinearSolver::TpetraLinearSolver(TpetraBaseLinearSystem & linsys, std::shared_ptr<tftk::linsys::MeshManagerNonFieldBased> & manager)
  :
    linearSystem(linsys),
    meshManager(manager),
    params_(Teuchos::rcp(new Teuchos::ParameterList)),
    paramsPrecond_(Teuchos::rcp(new Teuchos::ParameterList))
  {
    const double tolerance = 1e-3;
    const int max_iterations = 20;
    const int kspace = 500;
    params_->set("Convergence Tolerance", tolerance);
    params_->set("Maximum Iterations", max_iterations);

    if (0)
      {
        params_->set("Verbosity", Belos::Debug + Belos::Warnings + Belos::IterationDetails
                     + Belos::OrthoDetails + Belos::FinalSummary
                     + Belos::TimingDetails + Belos::StatusTestDetails);
      }

    params_->set("Output Frequency", 1);
    Teuchos::RCP<std::ostream> belosOutputStream = Teuchos::rcpFromRef(std::cout);
    params_->set("Output Stream", belosOutputStream);
    params_->set("Num Blocks", kspace);
    params_->set("Maximum Restarts", std::max(1,max_iterations/kspace));
    std::string orthoType = "ICGS";
    params_->set("Orthogonalization",orthoType);
    params_->set("Implicit Residual Scaling", "Norm of Preconditioned Initial Residual");

    paramsPrecond_->set("relaxation: type","Symmetric Gauss-Seidel");
    paramsPrecond_->set("relaxation: sweeps",4);
    //  paramsPrecond_->set("relaxation: type","Jacobi");
    //  paramsPrecond_->set("relaxation: sweeps",1);


  }

  TpetraLinearSolver::~TpetraLinearSolver()
  {
    problem_ = Teuchos::null;
    preconditioner_ = Teuchos::null;
    solver_ = Teuchos::null;
  }

  void
  TpetraLinearSolver::setupLinearSolver(bool doPrecond)
  {
    problem_ = Teuchos::RCP<LinSys::LinearProblem>(new LinSys::LinearProblem(linearSystem.matrix(), linearSystem.solution(),
                                                                             linearSystem.vector()) );
    //  const std::string preconditionerType ("RELAXATION");
    //  preconditioner_ = Teuchos::rcp(new Ifpack2::Relaxation<tftk::linsys::BlockMatrix>(linearSystem.block_matrix()) );

    if (doPrecond)
      {
        Ifpack2::Factory factory;
        const std::string preconditionerType ("RELAXATION");

        preconditioner_ = factory.create(preconditionerType,
                                         Teuchos::rcp_const_cast<const tftk::linsys::RowMatrix>( Teuchos::rcp_dynamic_cast<tftk::linsys::RowMatrix>( linearSystem.matrix()) ),
                                         0);
        preconditioner_->setParameters(*paramsPrecond_);
        preconditioner_->initialize();
        problem_->setRightPrec(preconditioner_);
      }
    solver_ = Teuchos::RCP<LinSys::GmresSolver>(new LinSys::GmresSolver(problem_, params_) );

  }

  void
  TpetraLinearSolver::solve(stk::mesh::FieldBase * deltaQ)
  {

    if (!preconditioner_.is_null())
      preconditioner_->compute();
    problem_->setProblem();
    solver_->solve();

    meshManager->copy_tpetra_to_stk(linearSystem.solution(), deltaQ);
  }

  int
  TpetraLinearSolver::getLinearIterations()
  {
    return solver_->getNumIters();
  }

  double
  TpetraLinearSolver::getLinearProblemResidual()
  {
    tftk::linsys::Vector resid(linearSystem.vector()->getMap());

    linearSystem.matrix()->apply(*linearSystem.solution(), resid);

    tftk::linsys::OneDVector rhs = linearSystem.vector()->get1dViewNonConst();
    tftk::linsys::OneDVector res = resid.get1dViewNonConst ();
    for (int i=0; i<rhs.size(); ++i)
      res[i] -= rhs[i];

    return (resid.norm2()/linearSystem.vector()->norm2());

  }

}
#endif
