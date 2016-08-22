//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Spatiocyte package
//
//        Copyright (C) 2006-2009 Keio University
//        Copyright (C) 2010-2014 RIKEN
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Spatiocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Spatiocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Spatiocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// 
//
//

#include <math.h>
#include <curand_kernel.h>
#include <thrust/execution_policy.h>
#include <Common.hpp>
#include <Random.hpp>

struct generate {
  __host__ __device__ generate(const unsigned _a, const unsigned _b):
    a(_a), b(_b) {;} 
  __device__ float operator()(const unsigned n) const {
    curandState s;
    curand_init(n, 0, 0, &s);
    float ranf(curand_uniform(&s));
    ranf *= (b - a + 0.999999);
    ranf += a;
    return (unsigned)truncf(ranf);
  }
  unsigned a, b;
};

Random::Random(const unsigned min, const unsigned max, const unsigned size,
               const unsigned seed):
  max_(max),
  min_(min),
  size_(size),
  cnt_(0),
  seed_(seed),
  data_(size) {
  initialize();
}

void Random::initialize() { 
  cnt_ = 0;
  thrust::counting_iterator<unsigned> begin(seed_);
  thrust::transform(thrust::device, begin, begin+size_, data_.begin(), 
                    generate(min_, max_));
  seed_ += size_;
}

umol_t Random::ran() {
  if(cnt_ >= size_) {
    initialize();
  }
  return data_[cnt_++];
}
