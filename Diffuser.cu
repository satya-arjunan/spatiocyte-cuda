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
// written by Satya Arjunan <satya.arjunan@gmail.com>
//

#include <time.h>
#include <thrust/execution_policy.h>
#include <thrust/random.h>
#include <curand_kernel.h>
#include <Diffuser.hpp>
#include <Compartment.hpp>
#include <Model.hpp>
#include <random>

Diffuser::Diffuser(const double D, Species& species):
  D_(D),
  species_(species),
  compartment_(species_.get_compartment()),
  box_mols_(species_.get_box_mols()),
  species_id_(species_.get_id()),
  vac_id_(species_.get_vac_id()),
  vac_xor_(species_.get_vac_xor()),
  //rng_(0, 11, 1000000*bool(D), time(0)),
  seed_(0),
  offsets_(ADJS*4),
  mols_(0) {
    std::cout << "mols size:" << mols_.size() << std::endl;
  }

void Diffuser::initialize() {
  box_voxels_ = compartment_.get_lattice().get_box_voxels();
  nbit_ = species_.get_compartment().get_model().get_nbit();
  one_nbit_ = pow(2, nbit_)-1;
  std::cout << species_.get_name_id() << std::endl;
  std::cout << "nbit:" << nbit_ << std::endl;
  std::cout << "one_nbit:" << one_nbit_ << std::endl;
  std::cout << "vac_id:" << vac_id_ << std::endl;
  std::cout << "vac_xor:" << vac_xor_ << std::endl;
}

void Diffuser::populate() {
  if(!D_) { 
    return;
  }
  mols_.resize(box_mols_[0].size());
  thrust::copy(box_mols_[0].begin(), box_mols_[0].end(), mols_.begin());
  thrust::sort(mols_.begin(), mols_.end());
  //col=even, layer=even
  offsets_[0] = -1;
  offsets_[1] = 1;
  offsets_[2] = -NUM_ROW-1;
  offsets_[3] = -NUM_ROW;
  offsets_[4] = NUM_ROW-1;
  offsets_[5] = NUM_ROW;
  offsets_[6] = -NUM_COLROW-NUM_ROW;
  offsets_[7] = -NUM_COLROW-1;
  offsets_[8] = -NUM_COLROW;
  offsets_[9] = NUM_COLROW-NUM_ROW;
  offsets_[10] = NUM_COLROW-1;
  offsets_[11] = NUM_COLROW;

  //col=even, layer=odd +24 = %layer*24
  offsets_[24] = -1;
  offsets_[25] = 1;
  offsets_[26] = -NUM_ROW;
  offsets_[27] = -NUM_ROW+1;
  offsets_[28] = NUM_ROW;
  offsets_[29] = NUM_ROW+1;
  offsets_[30] = -NUM_COLROW;
  offsets_[31] = -NUM_COLROW+1;
  offsets_[32] = -NUM_COLROW+NUM_ROW;
  offsets_[33] = NUM_COLROW;
  offsets_[34] = NUM_COLROW+1;
  offsets_[35] = NUM_COLROW+NUM_ROW;

  //col=odd, layer=even +12 = %col*12
  offsets_[12] = -1;
  offsets_[13] = 1;
  offsets_[14] = -NUM_ROW;
  offsets_[15] = -NUM_ROW+1;
  offsets_[16] = NUM_ROW;
  offsets_[17] = NUM_ROW+1;
  offsets_[18] = -NUM_COLROW-NUM_ROW;
  offsets_[19] = -NUM_COLROW;
  offsets_[20] = -NUM_COLROW+1;
  offsets_[21] = NUM_COLROW-NUM_ROW;
  offsets_[22] = NUM_COLROW;
  offsets_[23] = NUM_COLROW+1;

  //col=odd, layer=odd +36 = %col*12 + %layer*24
  offsets_[36] = -1;
  offsets_[37] = 1;
  offsets_[38] = -NUM_ROW-1;
  offsets_[39] = -NUM_ROW;
  offsets_[40] = NUM_ROW-1;
  offsets_[41] = NUM_ROW;
  offsets_[42] = -NUM_COLROW-1;
  offsets_[43] = -NUM_COLROW; //a
  offsets_[44] = -NUM_COLROW+NUM_ROW;
  offsets_[45] = NUM_COLROW-1;
  offsets_[46] = NUM_COLROW;
  offsets_[47] = NUM_COLROW+NUM_ROW;
}


__host__ __device__
unsigned int hash(unsigned int a)
{
    a = (a+0x7ed55d16) + (a<<12);
    a = (a^0xc761c23c) ^ (a>>19);
    a = (a+0x165667b1) + (a<<5);
    a = (a+0xd3a2646c) ^ (a<<9);
    a = (a+0xfd7046c5) + (a<<3);
    a = (a^0xb55a4f09) ^ (a>>16);
    return a;
}


/*
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
*/

struct generate {
  __host__ __device__ generate(mol_t* _offsets):
    offsets(_offsets) {} 
  __device__ umol_t operator()(const unsigned n, const umol_t vdx) const {
    curandState s;
    curand_init(n, 0, 0, &s);
    float ranf(curand_uniform(&s)*11.999999);
    const unsigned rand((unsigned)truncf(ranf));
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets[rand+(24&(-odd_lay))+(12&(-odd_col))]);
    if(val < 0 || val > NUM_VOXEL) {
      return vdx;
    }
    return val;
  }
  mol_t* offsets;
};

/*
struct generate {
  __host__ __device__ generate(mol_t* _offsets):
    offsets(_offsets) {} 
  __device__ umol_t operator()(const unsigned n, const umol_t vdx) const {
    thrust::default_random_engine rng(hash(n));
    thrust::uniform_int_distribution<int> uniform(0, 11);
    const bool odd_lay((vdx/NUM_COLROW)&1);
    const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
    mol2_t val(mol2_t(vdx)+offsets[uniform(rng)+(24&(-odd_lay))+
        (12&(-odd_col))]);
    if(val < 0 || val > NUM_VOXEL) {
      return vdx;
    }
    return val;
  }
  mol_t* offsets;
};
*/

void Diffuser::walk() {
  tars_.resize(mols_.size());
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+mols_.size()),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::sort(thrust::device, tars_.begin(), tars_.end());
  thrust::device_vector<umol_t> collisions(mols_.size());
  thrust::set_intersection(thrust::device, mols_.begin(), mols_.end(), tars_.begin(),
      tars_.end(), collisions.begin());
  //if(!collisions.size()) { 
    thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //}
  thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}

