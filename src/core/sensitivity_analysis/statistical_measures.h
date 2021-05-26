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

// This header contains functions to compute standard statistical measures for
// random variables. They can be used to evaluate results of the sensitivity
// analysis.

#ifndef STATISTICAL_MEASSURES_H_
#define STATISTICAL_MEASSURES_H_

#include <vector>

namespace bdm {

/// Compute the average value of the vector x.
/// \f[ \bar{x} = \frac{1}{dim(x)} \sum_{i=1}^{dim(x)} x_i \f]
double ComputeAverage(const std::vector<double>& x);

/// Compute the variance of the vector x having a an average of x_avg.
/// \f[ \mbox{var}(x) = \sum_{i=1}^{dim(x)} (x_i - x_{avg})^2 \f]
double ComputeVariance(const std::vector<double>& x, double x_avg);

/// Compute the covariance of the two vectors x and y having the averages x_avg
/// and y_avg.
/// \f[ \mbox{cov}(x,y) = \sum_{i=1}^{dim(x)} (x_i-x_{avg}) (y_i-y_{avg}) \f]
double ComputeCovariance(const std::vector<double>& x, double x_avg,
                         const std::vector<double>& y, double y_avg);

/// Computes the rank transformation of the vector x. It computes a vector whose
/// i-th entry gives the rank of x[i], i.e. x[i]'s position in a sorted list.
/// E.g. { 0.1, 0.4, 0.2, 0.0 } -> { 2.0, 3.0, 1.0, 0.0 }.
std::vector<double> RankTransfromation(const std::vector<double>& x);

/// Computes the Pearson Correlation Coefficient between two vectors x and y.
/// The Coefficient is helpful to test for linear realationships.
/// \f[ \mbox{CC}(X,Y) = \frac{\mbox{cov}(x,y)}
/// {\sqrt{\mbox{var}(x)\mbox{var}(y)}} \f]
double PearsonCorrelationCoefficient(const std::vector<double>& x,
                                     const std::vector<double>& y);

/// Not yet implemented.
double PartialCorrelationCoefficient(const std::vector<double>& x,
                                     const std::vector<double>& y);

/// Not yet implemented.
double StandardizedRegressionCoefficient(const std::vector<double>& x,
                                         const std::vector<double>& y);

/// Computes the Spearman Rank Correlation Coefficient of two vectors x and y.
/// It is equivalent to rank-transform x and y and computing the Pearson
/// Correlation Coefficient. It should be used for non-linear but monotonic
/// relationships.
double SpearmanRankCorrelationCoefficient(const std::vector<double>& x,
                                          const std::vector<double>& y);

/// Not yet implemented.
double PartialRankCorrelationCoefficient(const std::vector<double>& x,
                                         const std::vector<double>& y);

/// Not yet implemented.
double StandardizedRankRegressionCoefficient(const std::vector<double>& x,
                                             const std::vector<double>& y);

}  // namespace bdm

#endif  // STATISTICAL_MEASSURES_H_
