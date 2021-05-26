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

#include "core/sensitivity_analysis/latin_hypercube.h"
#include <algorithm>
#include <random>
#include "Math/DistFuncMathCore.h"
#include "core/util/log.h"
#include "core/util/random.h"

namespace bdm {

std::vector<double>& LatinHypercubeSampleUniform(int n_samples, double min,
                                                 double max) {
  static std::vector<double> samples(n_samples);
  double interval = (max - min) / static_cast<double>(n_samples);
  // For the Latin Hypercube, we split [min, max] into n_samples intervals and
  // pick a random number inside the interval.
  for (size_t i = 0; i < n_samples; i++) {
    UniformRng rng(min + i * interval, min + (i + 1) * interval);
    samples[i] = rng.Sample();
  }
  // We shuffle the vector such that our values occur in a random instead of an
  // increasing order.
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(samples.begin(), samples.end(), g);
  return samples;
};

std::vector<double>& LatinHypercubeSampleNormal(size_t n_samples, double mean,
                                                double sigma,
                                                size_t steps_per_sample) {
  // We use boundaries of 5 sigma to have a reasonable sampling domain.
  const double min = mean - 5 * sigma;
  const double max = mean + 5 * sigma;
  double low = min;
  double high = min;
  const double step =
      (max - min) / static_cast<double>(n_samples * steps_per_sample);
  static std::vector<double> samples;
  double cntr{1.0};
  // For the Latin Hpercube, we need to draw samples from regions of the
  // probability distribution that are equally likely. For us, this likelihood
  // is 1/n_samples. To find the appropriate regions, we use the normal
  // cumulative distribution function. We move from min to max and determine
  // the intervals such that normal_cft(high) = n/n_samples and
  // normal_cft(low) = (n-1)/n_samples, where n is in [1,..,n_samples].
  // Afterwards, we choose a uniform number between low and high.
  while (high < max) {
    high += step;
    if (ROOT::Math::normal_cdf(high, sigma, mean) > (cntr / n_samples)) {
      UniformRng rng(low, high);
      samples.push_back(rng.Sample());
      low = high;
      cntr += 1.0;
    } else {
      ;
    }
  }
  // For the last interval, the normal_cdf(5*sigma,mean,sigma) won't be larger
  // than 1, so we explictly draw a random number from the last interval
  // [low,max].
  if (samples.size() == n_samples - 1) {
    UniformRng rng(low, max);
    samples.push_back(rng.Sample());
  }
  // Catch error if we somehow don't get the right number of samples.
  else if (samples.size() < n_samples - 1) {
    Log::Fatal(
        "LatinHypercubeSampleNormal", "Unable to create sufficient",
        "data points. Consider increasing the parameter steps_per_samples.");
  } else if (samples.size() == n_samples) {
    ;  // all good
  }
  // In this code block, we catch errors that we do not expect.
  else {
    Log::Fatal("LatinHypercubeSampleNormal", "Unknown error, wrong number",
               "off data points. Consider increasing the parameter "
               "steps_per_samples.");
  }
  // We shuffle the vector such that our values occur in a random instead of an
  // increasing order.
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(samples.begin(), samples.end(), g);
  return samples;
}

}  // namespace bdm
