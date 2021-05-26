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

// This header contains functions to sample from probability distributions with
// the Latin Hypercube Methodology. See https://doi.org/10.2307/1268522

#ifndef LATIN_HYPERCUBE_H_
#define LATIN_HYPERCUBE_H_

#include <vector>

namespace bdm {

/// Get n_samples from a Uniform distribution [min, max] sampled in a Latin
/// Hypercube fashion.
std::vector<double>& LatinHypercubeSampleUniform(int n_samples, double min,
                                                 double max);

/// Get n_samples from a Normal distribution characterized by (mean, sigma)
/// sampled in a Latin Hypercube fashion. Note that for large number of samples,
/// the default value for steps_per_samples might not be large enough. As a
/// rule of thumb, steps_per_samples = n_samples / 10 worked in most cases.
std::vector<double>& LatinHypercubeSampleNormal(int n_samples, double mean,
                                                double sigma,
                                                size_t steps_per_sample = 1000);

}  // namespace bdm

#endif  // LATIN_HYPERCUBE_H_
