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


#ifndef __Lattice_hpp
#define __Lattice_hpp

#include <Common.hpp>

class Lattice {
 public: 
  Lattice(const unsigned, const Vector<unsigned>&, const unsigned);
  ~Lattice() {}
  void initialize();
  bool is_mol_at_box_edge(const umol_t) const;
  unsigned get_num_box() const;
  unsigned get_num_box_voxel() const;
  unsigned get_num_voxel() const;
  unsigned box_coord_to_mol(const unsigned, const Coord) const;
  unsigned box_mol_to_mol(const unsigned, const umol_t) const;
  unsigned coord_to_mol(const Vector<unsigned>&) const;
  umol_t box_coord_to_box_mol(const Coord) const;
  Coord box_mol_to_box_coord(const umol_t) const;
  Vector<unsigned> box_mol_to_coord(const umol_t) const;
  Vector<unsigned> box_to_coord(const unsigned) const;
  const Vector<unsigned>& get_dimensions() const;
  const Vector<unsigned>& get_box_dimensions() const;
  const Vector<unsigned>& get_box_voxel_dimensions() const;
  voxel_t* get_voxels();
  voxel_t** get_box_voxels();
 private:
  const unsigned num_box_;
  const Vector<unsigned> box_dimensions_;
  const Vector<unsigned> box_voxel_dimensions_;
  const Vector<unsigned> dimensions_;
  const unsigned num_box_voxel_;
  voxel_t* voxels_;
  voxel_t** box_voxels_;
};

#endif /* __Lattice_hpp */

