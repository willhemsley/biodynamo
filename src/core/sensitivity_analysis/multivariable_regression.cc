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

#include "core/sensitivity_analysis/multivariable_regression.h"
#include <TF1.h>
#include <TFitter.h>
#include <TLinearFitter.h>
#include <TRandom3.h>
#include <string>

namespace bdm {

void MultivariableRegression::Fit() {
  // Get the number of input variables and data points
  auto n_points = x_.GetNrows();
  auto n_variables = x_.GetNcols();

  // Initialize the TLinearFitter - i.e. use ROOT for fitting
  TLinearFitter lf = TLinearFitter(n_variables);

  // Workaround to set "hyp<n_variables>" as formula
  std::string hyperplane = "hyp" + std::to_string(n_variables);
  char* hyperplane_char;
  hyperplane_char = &hyperplane[0];
  lf.SetFormula(hyperplane_char);

  // Assing the data to the fitters
  lf.AssignData(n_points, n_variables, x_.GetMatrixArray(),
                y_.GetMatrixArray());

  // Compute regression parameters
  if (robust_) {
    lf.EvalRobust();
  } else {
    lf.Eval();
  }

  // Save regression parameteres into params_
  lf.GetParameters(params_);
}

TVectorD& MultivariableRegression::GetParameters() { return params_; }

TVectorD& MultivariableRegression::GetResiduals() {
  // Check if ressiduals have been computed, if not compute them
  if (residuals_.GetNoElements() == 0) {
    ComputeResiduals_();
  }
  return residuals_;
}

void MultivariableRegression::ComputeResiduals_() {
  // Resize residuals_ and fill all elements with 0
  residuals_.ResizeTo(y_.GetNoElements());
  residuals_.Zero();
  double* x_ptr_ = x_.GetMatrixArray();
  int n_variables_ = x_.GetNcols();
  // Compute the residuals
  for (int i = 0; i < y_.GetNoElements(); i++) {
    residuals_[i] = y_[i] - ComputeResult_(x_ptr_);
    x_ptr_ += n_variables_;
  }
}

double MultivariableRegression::ComputeResult_(const double* x) {
  double result{0.0};
  result += params_[0];
  for (int i = 1; i < params_.GetNoElements(); i++) {
    result += x[i - 1] * params_[i];
  }
  return result;
}

}  // namespace bdm
