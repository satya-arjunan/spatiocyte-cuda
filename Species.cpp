//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Spatiocyte package
//
//        Copyright (C) 2006-2009 Keio University
//        Copyright (C) 2010-2013 RIKEN
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

#include <iostream>
#include <string>
#include <time.h>
#include <Species.hpp>
#include <Compartment.hpp>
#include <Model.hpp>

Species::Species(const std::string name, const unsigned nmols, const double D,
    Model& model, Compartment& compartment, Species& vacant,
    const bool is_structure_species)
    : compartment_(compartment),
      vacant_(vacant),
      name_(get_init_name(name)),
      init_nmols_(nmols),
      is_structure_species_(is_structure_species),
      id_(model.push_species(*this)),
      vac_id_(vacant_.get_id()),
      vac_xor_(vac_id_^id_),
      diffuser_(D, *this),
      rng_(time(0)) {
  get_box_mols().resize(get_compartment().get_lattice().get_num_box());
}

void Species::initialize() {
  nbit_ = get_compartment().get_model().get_nbit();
  diffuser_.initialize();
}

void Species::populate() {
  /*
  for(unsigned box(0), n(get_box_mols().size()); box != n; ++box) {
    for(unsigned i(0); i != init_nmols_; ++i) {
      umol_t mol(vacant_.get_random_valid_mol(box));
      //Ensure mols are not populated on the box edges.
      //When getting target coords during walk, we assume mols have additional
      //space at the edge to be added or subtracted.
      //This is needed because the coords are unsigned values.
      while(get_compartment().get_lattice().is_mol_at_box_edge(mol)) {
        mol = vacant_.get_random_valid_mol(box);
      }
      populate_mol(box, mol);
    }
  }
  */
  for(unsigned box(0), n(get_box_mols().size()); box != n; ++box) {
    for(unsigned i(0); i != init_nmols_; ++i) {
      populate_mol(box, vacant_.get_random_valid_mol(box));
    }
  }
  diffuser_.populate();
}

void Species::populate_mol(const unsigned box, const umol_t vdx) {
  /*
  get_compartment().get_lattice().get_box_voxels()[box][vdx] = get_id();
  Coord clr(get_compartment().get_lattice().box_mol_to_box_coord(vdx));
  get_box_mols()[box].push_back(clr);
  */
  get_compartment().get_lattice().get_box_voxels()[box][vdx] = get_id();
  get_box_mols()[box].push_back(vdx);
}

umol_t Species::get_random_valid_mol(const unsigned box) {
  /*
  std::vector<Coord>& mols(get_box_mols()[box]);
  umol_t mol(get_compartment().get_lattice().box_coord_to_box_mol(
                                  mols[rng_.IRan(0, mols.size())]));
  while(get_compartment().get_lattice().get_box_voxels()[box][mol] != 
        get_id()) {
    mol = get_compartment().get_lattice().box_coord_to_box_mol(
                                 mols[rng_.IRan(0, mols.size())]);
  }
  return mol;
  */
  std::vector<umol_t>& mols(get_box_mols()[box]);
  umol_t mol(mols[rng_.ran(0, mols.size())]);
  while(get_compartment().get_lattice().get_box_voxels()[box][mol] !=
        get_id()) {
    mol = mols[rng_.ran(0, mols.size())];
  }
  return mol;
}

bool Species::is_structure_species() const {
  return is_structure_species_;
}

bool Species::is_root_structure_species() const {
  return (this == &vacant_);
}

voxel_t Species::get_id() const {
  return id_;
}

voxel_t Species::get_vac_id() const {
  return vac_id_;
}

voxel_t Species::get_vac_xor() const {
  return vac_xor_;
}

Diffuser& Species::get_diffuser() {
  return diffuser_;
}

Compartment& Species::get_compartment() {
  return compartment_;
}

Species& Species::get_vacant() {
  return vacant_;
}

const std::string& Species::get_name() const {
  return name_;
}

const std::string Species::get_name_id() const {
  std::stringstream sid;
  sid << get_id();
  return std::string(get_name()+":"+sid.str());
  //return std::string(get_name()+":id:"+std::to_string(get_id())); c++11
}

const std::string Species::get_init_name(const std::string name) const {
  if(is_root_structure_species()) {
    return std::string(compartment_.get_name()+"/"+name);
  }
  return std::string(vacant_.get_name()+"/"+name);
}

std::vector<unsigned> Species::get_mols() {
  /*
  std::vector<unsigned> mols;
  for(unsigned box(0), m(get_box_mols().size()); box != m; ++box) {
    for(unsigned i(0), n(get_box_mols()[box].size()); i != n; ++i) {
      mols.push_back(get_compartment().get_lattice().box_coord_to_mol(box,
          get_box_mols()[box][i]));
    }
  }
  return mols;
  */
  std::vector<unsigned> mols;
  for(unsigned box(0), m(get_box_mols().size()); box != m; ++box) {
    for(unsigned i(0), n(get_box_mols()[box].size()); i != n; ++i) {
      mols.push_back(get_compartment().get_lattice().box_mol_to_mol(box,
          get_box_mols()[box][i]));
    }
  }
  return mols;
}

std::vector<std::vector<umol_t> >& Species::get_box_mols() {
  return box_mols_;
}
