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

#include "core/sensitivity_analysis/statistical_measures.h"
#include <cmath>
#include <numeric>
#include "core/util/log.h"

namespace bdm {

double ComputeAverage(const std::vector<double>& x) {
  double sum = std::accumulate(x.begin(), x.end(), 0.0);
  return sum / static_cast<double>(x.size());
};

double ComputeVariance(const std::vector<double>& x, const double x_avg) {
  double variance{0.0};
  for (double value : x) {
    variance += pow(value - x_avg, 2);
  }
  variance /= static_cast<double>(x.size());
  return variance;
};

double ComputeCovariance(const std::vector<double>& x, const double x_avg,
                         const std::vector<double>& y, const double y_avg) {
  if (x.size() != y.size()) {
    Log::Fatal("ComputeCovariance", "Received vectors with different sizes ",
               x.size(), " and ", y.size());
  }
  double covariance{0.0};
  for (size_t i = 0; i < x.size(); i++) {
    covariance += (x[i] - x_avg) * (y[i] - y_avg);
  }
  covariance /= static_cast<double>(x.size());
  return covariance;
};

std::vector<double> RankTransfromation(const std::vector<double>& x) {
  std::vector<int> index(x.size());
  std::vector<double> output(x.size());
  std::iota(index.begin(), index.end(), 0);
  std::sort(index.begin(), index.end(),
            [&](const int& a, const int& b) { return (x[a] < x[b]); });
  for (size_t i = 0; i < index.size(); i++) {
    output[index[i]] = i + 1;
  }
  return output;
};

double PearsonCorrelationCoefficient(const std::vector<double>& x,
                                     const std::vector<double>& y) {
  const double x_avg = ComputeAverage(x);
  const double y_avg = ComputeAverage(y);
  const double x_var = ComputeVariance(x, x_avg);
  const double y_var = ComputeVariance(y, y_avg);
  const double xy_cov = ComputeCovariance(x, x_avg, y, y_avg);
  return xy_cov / sqrt(x_var * y_var);
};

double PartialCorrelationCoefficient(const std::vector<double>& x,
                                     const std::vector<double>& y) {
  Log::Warning("PartialCorrelationCoefficient", "Not implemented, returns 0.");
  return 0.0;
};

double StandardizedRegressionCoefficient(const std::vector<double>& x,
                                         const std::vector<double>& y) {
  Log::Warning("StandardizedRegressionCoefficient",
               "Not implemented, returns 0.");
  return 0.0;
};

double SpearmanRankCorrelationCoefficient(const std::vector<double>& x,
                                          const std::vector<double>& y) {
  const auto x_ranked = RankTransfromation(x);
  const auto y_ranked = RankTransfromation(y);
  return PearsonCorrelationCoefficient(x_ranked, y_ranked);
};

double PartialRankCorrelationCoefficient(const std::vector<double>& x,
                                         const std::vector<double>& y) {
  Log::Warning("PartialRankCorrelationCoefficient",
               "Not implemented, returns 0.");
  return 0.0;
}

double StandardizedRankRegressionCoefficient(const std::vector<double>& x,
                                             const std::vector<double>& y) {
  Log::Warning("PartialRankCorrelationCoefficient",
               "Not implemented, returns 0.");
  return 0.0;
}

}  // namespace bdm
