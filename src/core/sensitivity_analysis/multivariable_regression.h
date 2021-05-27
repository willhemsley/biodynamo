// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & Newcastle University for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#ifndef MULTIVARIABLE_REGRESSION_H_
#define MULTIVARIABLE_REGRESSION_H_

#include <TMatrixD.h>
#include <TVectorD.h>
// #include <TLinearFitter.h>

namespace bdm {

class MultivariableRegression {
 public:
  /// A class for a multivariable regression (Wrapper around TLinearFitter).
  /// param x: random input variables, param y: prediction y(x), param robust:
  /// use a robust estimator or not.
  MultivariableRegression(TMatrixD& x, TVectorD& y, bool robust)
      : y_{y}, x_{x}, robust_{robust} {}

  /// This function fits a multivariable regression as
  /// \f[ \hat{y} = \alpha + \sum_{i=1}^N \beta_i x_i \f]
  /// and determines the parameters alpha and beta.
  void Fit();

  /// Returns the residuals, i.e. the error of the fit
  /// \f[ y - \hat{y} \f]
  TVectorD& GetResiduals();

  /// Returns the parametes alpha and beta computed with
  /// MultivariableRegression::Fit() .
  TVectorD& GetParameters();

 private:
  TVectorD y_;
  TMatrixD x_;
  TVectorD params_;
  TVectorD residuals_;
  bool robust_;
  double ComputeResult_(const double* x);
  void ComputeResiduals_();
};

}  // namespace bdm

#endif  // MULTIVARIABLE_REGRESSION_H_
