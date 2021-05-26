// -----------------------------------------------------------------------------
//
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include "core/sensitivity_analysis/latin_hypercube.h"
#include "core/sensitivity_analysis/statistical_measures.h"
#include "gtest/gtest.h"

namespace bdm {

void print_vector(std::vector<double>& x) {
  int cntr = 0;
  std::cout << "[ ";
  for (size_t i = 0; i < x.size(); i++) {
    std::cout << x[i] << " ";
    cntr += 1;
    if (cntr % 5 == 0) {
      std::cout << "\n";
    }
  }
  std::cout << " ]" << std::endl;
}

// Random vector generated with numpy
std::vector<double> vec1{1.23190313,  0.047162,   1.81700215,  -0.53846081,
                         -0.64368265, 1.00109012, 0.26790301,  1.7715311,
                         2.60447774,  0.07117996, 0.54367455,  -1.04162298,
                         0.61680193,  0.07535977, -0.79869194, -1.5391589,
                         0.85492167,  0.12051386, 0.70164365,  -1.8281372};
// Random vector generated with numpy
std::vector<double> vec2{-1.32098991, -0.8961675,  -1.50709042, 0.54661656,
                         2.21926183,  1.2490651,   0.52850182,  0.96283572,
                         2.07803291,  -0.92461447, 1.47606221,  0.0488394,
                         0.05216119,  -0.34312988, 0.41495332,  -0.84718729,
                         -0.71119848, 0.8292938,   -0.1783944,  -0.06716452};

TEST(SensitivityAnalysisTest, Average) {
  std::vector<double> vecA{0.1, 0.2, 0.4, 0.5};
  std::vector<double> vecB{1, 2, 3, 4, 6, 7, 8, 9};
  auto vecA_avg = ComputeAverage(vecA);
  auto vecB_avg = ComputeAverage(vecB);
  // Easy averages checked with double precision.
  EXPECT_DOUBLE_EQ(0.3, vecA_avg);
  EXPECT_DOUBLE_EQ(5, vecB_avg);
}

TEST(SensitivityAnalysisTest, Variance) {
  auto vec1_avg = ComputeAverage(vec1);
  EXPECT_FLOAT_EQ(0.26677050872772495, vec1_avg);
  auto vec1_var = ComputeVariance(vec1, vec1_avg);
  EXPECT_FLOAT_EQ(1.222814885932807, vec1_var);
}

TEST(SensitivityAnalysisTest, Covariance) {
  auto vec1_avg = ComputeAverage(vec1);
  EXPECT_FLOAT_EQ(0.26677050872772495, vec1_avg);
  auto vec2_avg = ComputeAverage(vec2);
  EXPECT_FLOAT_EQ(0.180484349082455, vec2_avg);
  auto cov = ComputeCovariance(vec1, vec1_avg, vec2, vec2_avg);
  EXPECT_FLOAT_EQ(0.12851629444985857, cov);
}

TEST(SensitivityAnalysisTest, RankTransfromation) {
  std::vector<double> result1{17, 7, 19, 6, 5, 16, 11, 18, 20, 8,
                              12, 3, 13, 9, 4, 2,  15, 10, 14, 1};
  std::vector<double> result2{2,  4,  1,  14, 20, 17, 13, 16, 19, 3,
                              18, 10, 11, 7,  12, 5,  6,  15, 8,  9};
  auto rank1 = RankTransfromation(vec1);
  auto rank2 = RankTransfromation(vec2);
  for (size_t i = 0; i < vec1.size(); i++) {
    EXPECT_DOUBLE_EQ(result1[i], rank1[i]);
    EXPECT_DOUBLE_EQ(result2[i], rank2[i]);
  }
}

TEST(SensitivityAnalysisTest, PearsonCorrelationCoefficient) {
  float cc = PearsonCorrelationCoefficient(vec1, vec2);
  EXPECT_FLOAT_EQ(0.11143819243035334, cc);
}

TEST(SensitivityAnalysisTest, SpearmanRankCorrelationCoefficient) {
  float srcc = SpearmanRankCorrelationCoefficient(vec1, vec2);
  EXPECT_FLOAT_EQ(0.05714285714285714, srcc);
}

TEST(SensitivityAnalysisTest, LHSUniform) {
  int n_samples{100};
  double min{1.0};
  double max{3.0};
  double analytic_mean{2.0};
  double analytic_var{1. / 3.};
  auto samples = LatinHypercubeSampleUniform(n_samples, min, max);
  auto mean = ComputeAverage(samples);
  auto var = ComputeVariance(samples, mean);
  auto minmax = std::minmax_element(samples.begin(), samples.end());
  // Test sample size
  EXPECT_EQ(n_samples, samples.size());
  // Test mean value of distribution
  EXPECT_NEAR(analytic_mean, mean, 0.01);
  // Test variance of distribution
  EXPECT_NEAR(analytic_var, var, 0.01);
  // Test if samples are sorted or not.
  EXPECT_FALSE(std::is_sorted(samples.begin(), samples.end()));
  // Test bounds for min and max
  EXPECT_TRUE(min < *minmax.first);
  EXPECT_TRUE(max + 5 * analytic_var > *minmax.second);
}

TEST(SensitivityAnalysisTest, LHSNormal) {
  int n_samples{1000};
  double analytic_mean{2.0};
  double analytic_var{1. / 3.};
  auto samples =
      LatinHypercubeSampleNormal(n_samples, analytic_mean, analytic_var);
  auto mean = ComputeAverage(samples);
  auto var = ComputeVariance(samples, mean);
  auto minmax = std::minmax_element(samples.begin(), samples.end());
  // Test sample size
  EXPECT_EQ(n_samples, samples.size());
  // Test mean value of distribution
  EXPECT_NEAR(analytic_mean, mean, 0.01);
  // Test variance of distribution
  EXPECT_NEAR(pow(analytic_var, 2), var, 0.01);
  // Test if samples are sorted or not.
  EXPECT_FALSE(std::is_sorted(samples.begin(), samples.end()));
  // Test bounds for min and max
  EXPECT_TRUE(analytic_mean - 5 * analytic_var < *minmax.first);
  EXPECT_TRUE(analytic_mean + 5 * analytic_var > *minmax.second);
}

}  // namespace bdm
