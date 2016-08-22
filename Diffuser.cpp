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
  rng_(0, 11, 10000000, time(0)) {}

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


//t(10000, 16) = 56 s
//t(10000, 1) = 1.87 s
void Diffuser::walk() {
  for(unsigned box(0), n(box_mols_.size()); box != n; ++box) {
    walk(box_mols_[box].data(), box_mols_[box].size());
  }
}

void Diffuser::walk(umol_t* mols, const unsigned size) {
  for (unsigned i(0); i != size; ++i) {
    umol_t tar(compartment_.get_tar(mols[i], rng_.ran()));
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


