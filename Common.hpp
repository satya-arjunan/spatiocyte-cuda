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

#ifndef __Common_hpp
#define __Common_hpp

#include <iostream>
#include <vector>
#include <bitset>
#include <climits>
#include <stdint.h>
#include <emmintrin.h> //SSE2 intrinsics
#include <immintrin.h> //AVX2 intrinsics
#include <Vector.hpp>

class Compartment;
class Diffuser;
class Model;
class Species;
class Stepper;
class VisualLogger;

#define ADJS 12

typedef union
{
  __m128i m128i[4];
  __m256i m256i[2];
  uint8_t uint8[64];
  uint16_t uint16[32];
  uint32_t uint32[16];
  int8_t int8[64];
  int16_t int16[32];
  int32_t int32[16];
} union512;

typedef union
{
  __m128i m128i[2];
  __m256i m256i;
  uint8_t uint8[32];
  uint16_t uint16[16];
  uint32_t uint32[8];
  int8_t int8[32];
  int16_t int16[16];
  int32_t int32[8];
} union256;

typedef union
{
  __m128i m128i;
  uint8_t uint8[16];
  uint16_t uint16[8];
  uint32_t uint32[4];
  int8_t int8[16];
  int16_t int16[8];
  int32_t int32[4];
} union128;

typedef uint8_t voxel_t;
#define WORD (sizeof(voxel_t)*8)

typedef int16_t mol_t;
typedef uint16_t umol_t;
typedef int32_t mol2_t;
typedef uint32_t umol2_t;

struct Coord
{
  umol_t x:5;
  umol_t y:5;
  umol_t z:5;
};

struct CoordInt
{
  int x;
  int y;
  int z;
};

template<typename T>
void cout_binary(const T a, const std::string str) {
  std::cout << str << " msb -> lsb:" << std::endl;
  const char* beg(reinterpret_cast<const char*>(&a));
  const char* end(beg + sizeof(a));
  while(1) {
    std::cout << std::bitset<CHAR_BIT>(*--end) << ' ';
    if(end == beg) {
      break;
    }
  }
  std::cout << std::endl;
}

template<typename Register, typename CastType>
void cout_uint(const Register reg, const std::string title) {
  std::cout << title << " msb -> lsb:" << std::endl;
  for(int i(sizeof(Register)/sizeof(CastType)-1); i >= 0; --i) {
    std::cout << (uint32_t)((CastType*)&reg)[i] << " ";
  }
  std::cout << std::endl;
}

#endif /* __Common_hpp */
