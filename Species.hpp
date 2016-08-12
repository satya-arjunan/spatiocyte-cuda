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


#ifndef __Species_hpp
#define __Species_hpp

//#include <RandomLib/Random.hpp>
#include <Common.hpp>
#include <Diffuser.hpp>
#include <Random.hpp>

class Species
{ 
public: 
  Species(const std::string, const unsigned, const double, Model&, Compartment&,
      Species& vacant, const bool is_structure_species=false);
  ~Species() {}
  void initialize();
  void populate();
  void populate_mol(const unsigned, const umol_t);
  bool is_structure_species() const;
  bool is_root_structure_species() const;
  voxel_t get_id() const;
  voxel_t get_vac_id() const;
  voxel_t get_vac_xor() const;
  umol_t get_random_valid_mol(const unsigned);
  Diffuser& get_diffuser();
  Species& get_vacant();
  Compartment& get_compartment();
  const std::string& get_name() const;
  const std::string get_name_id() const;
  std::vector<unsigned> get_mols();
  std::vector<std::vector<umol_t> >& get_box_mols();
private:
  const std::string get_init_name(const std::string) const;
private:
  Compartment& compartment_;
  Species& vacant_;
  const std::string name_;
  const unsigned init_nmols_;
  const bool is_structure_species_;
  const voxel_t id_;
  const voxel_t vac_id_;
  const voxel_t vac_xor_;
  std::vector<std::vector<umol_t> > box_mols_;
  Diffuser diffuser_;
  Random rng_;
  voxel_t nbit_;
};

#endif /* __Species_hpp */

