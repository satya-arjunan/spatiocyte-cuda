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
//'iptables -I INPUT -p tcp --dport 1022 -j ACCEPT'
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
  mols_(species_.get_mols()),
  voxels_(species_.get_compartment().get_lattice().get_voxels()),
  offsets_(species_.get_compartment().get_offsets()),
  species_id_(species_.get_id()),
  vac_id_(species_.get_vac_id()),
  seed_(0) {
}

void Diffuser::initialize() {
  std::cout << "init diffuser of:" << species_.get_name_id() << std::endl;
}

double Diffuser::getD() {
  return D_;
}

struct generate {
  __host__ __device__ generate(const mol_t* _offsets):
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
  const mol_t* offsets;
};

struct is_occupied {
  __device__ bool operator()(const voxel_t voxel) {
    return ((bool)voxel);
  }
};

struct update {
  __device__ umol_t operator()(const umol_t mol) const {
    return mol;
  }
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> stencil(voxels_.begin(),
        tars_.begin());
  thrust::transform_if(
      mols_.begin(),
      mols_.end(),
      stencil,
      tars_.begin(),
      update(),
      is_occupied());
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> vacants(voxels_.begin(),
        mols_.begin());
  thrust::fill_n(thrust::device, vacants, size, false);
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> occupieds(voxels_.begin(),
        tars_.begin());
  thrust::fill_n(thrust::device, occupieds, size, true);
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  seed_ += size;
}

/*
//Simplified full collision with transform_if: very small performance
//improvement
struct is_occupied {
  __device__ bool operator()(const voxel_t voxel) {
    return ((bool)voxel);
  }
};

struct update {
  __device__ umol_t operator()(const umol_t mol) const {
    return mol;
  }
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> stencil(lattice_.begin(),
        tars_.begin());
  thrust::transform_if(
      mols_.begin(),
      mols_.end(),
      stencil,
      tars_.begin(),
      update(),
      is_occupied());
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> vacants(lattice_.begin(),
        mols_.begin());
  thrust::fill_n(thrust::device, vacants, size, false);
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> occupieds(lattice_.begin(),
        tars_.begin());
  thrust::fill_n(thrust::device, occupieds, size, true);
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += size;
}
*/


/*
//Full collision lattice with transform_if: same performance
struct is_vacant {
  __device__ bool operator()(const voxel_t voxel) {
    return (bool)voxel;
  }
};

struct update {
  __device__ umol_t operator()(const voxel_t tar, const umol_t mol) const {
    return mol;
  }
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> stencil(lattice_.begin(),
        tars_.begin());
  thrust::transform_if(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      stencil,
      tars_.begin(),
      update(),
      is_vacant());
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> vacants(lattice_.begin(),
        mols_.begin());
  thrust::fill_n(thrust::device, vacants, size, true);
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> occupieds(lattice_.begin(),
        tars_.begin());
  thrust::fill_n(thrust::device, occupieds, size, false);
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += size;
}
*/


/*
//Full collision check with lattice: 0.3 s (1,600,000 molecules, 42.5 s)
struct update {
  __host__ __device__ update(const voxel_t* _lattice):
    lattice(_lattice) {} 
  __device__ umol_t operator()(const umol_t tar, const umol_t mol) const {
    if(lattice[tar]) {
      return tar;
    }
    return mol;
  }
  const voxel_t* lattice;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::transform(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      tars_.begin(),
      update(thrust::raw_pointer_cast(&lattice_[0])));
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> vacants(lattice_.begin(),
        mols_.begin());
  thrust::fill_n(thrust::device, vacants, mols_.size(), true);
  thrust::permutation_iterator<thrust::device_vector<voxel_t>::iterator,
    thrust::device_vector<umol_t>::iterator> occupieds(lattice_.begin(),
        tars_.begin());
  thrust::fill_n(thrust::device, occupieds, tars_.size(), false);
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/

/*
//Partial collision check with lattice: 0.2 s
struct update {
  __host__ __device__ update(bool* _lattice):
    lattice(_lattice) {} 
  __device__ umol_t operator()(const umol_t tar, const umol_t mol) const {
    if(lattice[tar]) {
      lattice[tar] = false;
      lattice[mol] = true;
      return tar;
    }
    return mol;
  }
  bool* lattice;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::transform(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      mols_.begin(),
      update(thrust::raw_pointer_cast(&lattice_[0])));
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/

/*
//Full collision check using molecule list only: 4.3s

struct update {
  __host__ __device__ update(const size_t _size, const umol_t* _mols):
    size(_size), mols(_mols) {} 
  __device__ umol_t operator()(const umol_t tar, const umol_t mol) const {
    for(unsigned i(0); i != size; ++i) {
      if(tar == mols[i]) {
        return mol;
      }
    }
    return tar;
  }
  const size_t size;
  const umol_t* mols;
};

void Diffuser::walk() {
  const size_t size(mols_.size());
  tars_.resize(size);
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+size),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::transform(thrust::device, 
      tars_.begin(),
      tars_.end(),
      mols_.begin(),
      tars_.begin(),
      update(size, thrust::raw_pointer_cast(&mols_[0])));
  thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/


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


/*
//Collisions not check, with intersection: 7s
void Diffuser::walk() {
  tars_.resize(mols_.size());
  thrust::transform(thrust::device, 
      thrust::counting_iterator<unsigned>(seed_),
      thrust::counting_iterator<unsigned>(seed_+mols_.size()),
      mols_.begin(),
      tars_.begin(),
      generate(thrust::raw_pointer_cast(&offsets_[0])));
  thrust::sort(thrust::device, tars_.begin(), tars_.end());
  collisions_.resize(mols_.size());
  thrust::set_intersection(thrust::device, mols_.begin(), mols_.end(),
      tars_.begin(), tars_.end(), collisions_.begin());
  //if(!collisions.size()) { 
    thrust::copy(tars_.begin(), tars_.end(), mols_.begin());
  //}
  //thrust::copy(mols_.begin(), mols_.end(), box_mols_[0].begin());
  seed_ += mols_.size();
}
*/

/*
//Sequential original 10.5 s
void Diffuser::walk(umol_t* mols, const unsigned size) {
  for (unsigned i(0); i != size; ++i) {
    umol_t tar(compartment_.get_tar(mols[i], rng_.Ran16_12()));
    for(unsigned j(0); j != size; ++j) {
      if(mols[j] == tar) {
        goto next;
      }
    }
    mols[i] = tar;
next:
    continue;
  }
}
*/
