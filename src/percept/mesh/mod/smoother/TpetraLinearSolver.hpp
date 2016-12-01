// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*
 * TpetraLinearSolver.h
 *
 *  Created on: May 2, 2014
 *      Author: tcfishe
 *  Modified on: May 30, 2015
 *      Author: srkenno
 */

#ifndef TPETRALINEARSOLVER_HPP_
#define TPETRALINEARSOLVER_HPP_
#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT) && ENABLE_SMOOTHER3

#include <stk_mesh/base/Field.hpp>

#include <TpetraLinearSystemTypes.h>
#include <TpetraBaseLinearSystem.h>

// Forward declare templates
namespace Teuchos {

  template <typename T>
  class ArrayRCP;

  class ParameterList;

}

namespace Belos {

  template <typename Scalar, typename MultiVector>
  class MultiVecTraits;

  template <typename Scalar, typename MultiVector, typename Operator>
  class OperatorTraits;

  template <typename Scalar, typename MultiVector, typename Operator>
  class LinearProblem;

  template <typename Scalar, typename MultiVector, typename Operator>
  class SolverManager;

  template <typename Scalar, typename MultiVector, typename Operator>
  class PseudoBlockGmresSolMgr;

  template <typename Scalar, typename MultiVector, typename Operator>
  class TFQMRSolMgr;

}

namespace Ifpack2 {

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  class Preconditioner;

}

namespace percept {

  //class TpetraLinearSystem;
  //class MatrixFreeOperator;

  struct LinSys {

    typedef tftk::linsys::GlobalOrdinal   GlobalOrdinal; // MUST be signed
    typedef tftk::linsys::LocalOrdinal    LocalOrdinal;  // MUST be signed
    typedef tftk::linsys::Scalar          Scalar;

    typedef Belos::MultiVecTraits<Scalar, tftk::linsys::MultiVector>                                                MultiVectorTraits;
    typedef Belos::OperatorTraits<Scalar, tftk::linsys::MultiVector, tftk::linsys::Operator>                        OperatorTraits;
    typedef Belos::LinearProblem<Scalar, tftk::linsys::MultiVector, tftk::linsys::Operator>                              LinearProblem;
    typedef Belos::SolverManager<Scalar, tftk::linsys::MultiVector, tftk::linsys::Operator>                         SolverManager;
    typedef Belos::PseudoBlockGmresSolMgr<tftk::linsys::Scalar, tftk::linsys::MultiVector, tftk::linsys::Operator>  GmresSolver;
    typedef Belos::TFQMRSolMgr<tftk::linsys::Scalar, tftk::linsys::MultiVector, tftk::linsys::Operator>             TfqmrSolver;
    typedef Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, tftk::linsys::Node> Preconditioner;

  };


  class TpetraLinearSolver {

    typedef tftk::linsys::TpetraCrsLinearSystemImpl<tftk::linsys::MeshManagerNonFieldBased> TpetraBaseLinearSystem;
    TpetraBaseLinearSystem & linearSystem;
    std::shared_ptr<tftk::linsys::MeshManagerNonFieldBased> meshManager;

    const Teuchos::RCP<Teuchos::ParameterList> params_;
    const Teuchos::RCP<Teuchos::ParameterList> paramsPrecond_;
    Teuchos::RCP<LinSys::LinearProblem> problem_;
    Teuchos::RCP<LinSys::GmresSolver> solver_;
    Teuchos::RCP<LinSys::Preconditioner> preconditioner_;

  public:
    TpetraLinearSolver(TpetraBaseLinearSystem & linsys, std::shared_ptr<tftk::linsys::MeshManagerNonFieldBased> & manager);

    ~TpetraLinearSolver();

    void setupLinearSolver(bool doPrecond = true);

    //void setupMatrixFreeLinearSolver(Teuchos::RCP<MatrixFreeOperator> matrixFreeOperator, Teuchos::RCP<tftk::linsys::MultiVector> residual);

    void solve(stk::mesh::FieldBase * deltaQ);

    int getLinearIterations();

    double getLinearProblemResidual();


  };

}


#endif
#endif /* TPETRALINEARSOLVER_HPP_ */
