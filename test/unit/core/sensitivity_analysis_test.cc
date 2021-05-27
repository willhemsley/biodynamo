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

#include <TLinearFitter.h>
#include <TRandom3.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include "core/sensitivity_analysis/latin_hypercube.h"
#include "core/sensitivity_analysis/multivariable_regression.h"
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

TEST(SensitivityAnalysisTest, MultivariableRegression) {
  // Define parameters for test
  int n_samples{10000};
  int n_variables{5};
  double outlier_fraction{0.2};
  double noise_mean{0.0};
  double noise_sigma{0.01};
  time_t seed = time(NULL);
  TRandom3 rng;
  rng.SetSeed(seed);

  // Define parameters for function. We'll construct a hyperplane, basically
  // y = \alpha + \sum_i beta_i x_i. The alpha justifies the plus one below.
  double param_in[n_variables + 1];
  for (int i = 0; i < n_variables + 1; i++) {
    param_in[i] = rng.Uniform();
    std::cout << "parameter " << i << ": " << param_in[i] << std::endl;
  }

  // Produce random inputs and compute the output.
  TVectorD prediction(n_samples);
  TMatrixD inputs(n_samples, n_variables);
  for (int row = 0; row < n_samples; row++) {
    // Add a random noise of at most 5 sigma to our computation
    double noise = rng.Gaus(noise_mean, noise_sigma);
    if (noise_sigma < noise_mean - 5 * noise_sigma) {
      noise_sigma = noise_mean - 5 * noise_sigma;
    } else if (noise_sigma > noise_mean + 5 * noise_sigma) {
      noise_sigma = noise_mean + 5 * noise_sigma;
    } else {
      ;  // all good
    }
    prediction[row] = param_in[0] + noise;
    for (int col = 0; col < n_variables; col++) {
      inputs[row][col] = rng.Uniform();
      prediction[row] += (param_in[col + 1] * inputs[row][col]);
    }
  }

  // Compute MultivariateRegression
  MultivariableRegression m_reg_1(inputs, prediction, false);
  m_reg_1.Fit();
  auto param_fit = m_reg_1.GetParameters();

  // Compare input and output parameters
  for (int i = 0; i < param_fit.GetNoElements(); i++) {
    EXPECT_NEAR(param_in[i], param_fit[i], noise_sigma / 5.0);
  }

  // Check if residuals are in 5.5 sigma range - "Quality check"
  auto residuals = m_reg_1.GetResiduals();
  EXPECT_EQ(prediction.GetNoElements(), residuals.GetNoElements());
  for (int i = 0; i < residuals.GetNoElements(); i++) {
    EXPECT_TRUE(abs(residuals[i]) <= 5.5 * noise_sigma);
  }

  // Modify data to contain outliers
  for (int i = 0; i < prediction.GetNoElements();
       i += static_cast<int>(1.0 / outlier_fraction)) {
    prediction[i] += rng.Uniform();
  }

  // Compute MultivariateRegression on ill-data
  MultivariableRegression m_reg_2(inputs, prediction, false);
  m_reg_2.Fit();
  auto param_fit_regular = m_reg_2.GetParameters();

  // Compute robust MultivariateRegression on ill-data
  MultivariableRegression m_reg_robust(inputs, prediction, true);
  m_reg_robust.Fit();
  auto param_fit_robust = m_reg_robust.GetParameters();

  // Compare input and output parameters
  double err_regular{0.0};
  double err_robust{0.0};
  for (int i = 0; i < param_fit_robust.GetNoElements(); i++) {
    // Expect matching parameters
    EXPECT_NEAR(param_in[i], param_fit_robust[i], noise_sigma / 2.0);
    // Expect better estimates from robust fit over regular fit.
    err_robust += pow(param_in[i] - param_fit_robust[i], 2);
    err_regular += pow(param_in[i] - param_fit_regular[i], 2);
  }
  EXPECT_TRUE(err_robust < err_regular);
}

}  // namespace bdm
