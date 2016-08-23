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


#ifndef __Diffuser_hpp
#define __Diffuser_hpp

#include <thrust/device_vector.h>
#include <Random.hpp>
#include <Common.hpp>

class Diffuser
{ 
public: 
  Diffuser(const double, Species&);
  ~Diffuser() {}
  void initialize();
  void walk();
  void populate();
private:
  void set_offsets();
  double D_;
  Species& species_;
  Compartment& compartment_;
  std::vector<std::vector<umol_t> >& box_mols_;
  voxel_t** box_voxels_;
  const voxel_t species_id_;
  const voxel_t vac_id_;
  const voxel_t vac_xor_;
  //RandomGPU rng_;
  voxel_t nbit_;
  voxel_t one_nbit_;
  unsigned seed_;
  thrust::device_vector<mol_t> offsets_;
  thrust::device_vector<umol_t> mols_;
  thrust::device_vector<umol_t> tars_;
};

#endif /* __Diffuser_hpp */

