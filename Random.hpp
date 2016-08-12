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
// based on the SFMT code by Agner Fog
//

#ifndef __Random_hpp
#define __Random_hpp

#include <stdint.h>
#include <Common.hpp>

// Choose one of the possible Mersenne exponents.
// Higher values give longer cycle length and use more memory:
//#define MEXP   607
//#define MEXP  1279
//#define MEXP  2281
//#define MEXP  4253
//#define MEXP 11213
//#define MEXP 19937
#define MEXP 44497

// Define constants for the selected Mersenne exponent:
#if MEXP == 44497
#define SFMT_N    348 // Size of state vector
#define SFMT_M    330 // Position of intermediate feedback
#define SFMT_SL1    5 // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	  3 // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1    9 // Right shift of W[M], 32-bit words
#define SFMT_SR2	  3 // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xeffffffb,0xdfbebfff,0xbfbf7bef,0x9ffd7bff // AND mask
#define SFMT_PARITY 1,0,0xa3ac4000,0xecc1327a   // Period certification vector

#elif MEXP == 19937
#define SFMT_N    156
#define SFMT_M    122
#define SFMT_SL1   18
#define SFMT_SL2	  1
#define SFMT_SR1   11
#define SFMT_SR2	  1
#define SFMT_MASK	  0xdfffffef,0xddfecb7f,0xbffaffff,0xbffffff6
#define SFMT_PARITY 1,0,0,0x13c9e684

#elif MEXP == 11213
#define SFMT_N    88
#define SFMT_M    68
#define SFMT_SL1	14
#define SFMT_SL2	 3
#define SFMT_SR1	 7
#define SFMT_SR2	 3
#define SFMT_MASK	 0xeffff7fb,0xffffffef,0xdfdfbfff,0x7fffdbfd 
#define SFMT_PARITY 1,0,0xe8148000,0xd0c7afa3

#elif MEXP == 4253
#define SFMT_N    34
#define SFMT_M    17
#define SFMT_SL1	20
#define SFMT_SL2	 1
#define SFMT_SR1	 7
#define SFMT_SR2	 1
#define SFMT_MASK	 0x9f7bffff, 0x9fffff5f, 0x3efffffb, 0xfffff7bb
#define SFMT_PARITY 0xa8000001, 0xaf5390a3, 0xb740b3f8, 0x6c11486d

#elif MEXP == 2281
#define SFMT_N    18
#define SFMT_M    12
#define SFMT_SL1	19
#define SFMT_SL2	 1
#define SFMT_SR1	 5
#define SFMT_SR2	 1
#define SFMT_MASK	 0xbff7ffbf, 0xfdfffffe, 0xf7ffef7f, 0xf2f7cbbf
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x41dfa600

#elif MEXP == 1279
#define SFMT_N    10
#define SFMT_M     7
#define SFMT_SL1	14
#define SFMT_SL2	 3
#define SFMT_SR1	 5
#define SFMT_SR2	 1
#define SFMT_MASK	  0xf7fefffd, 0x7fefcfff, 0xaff3ef3f, 0xb5ffff7f
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x20000000

#elif MEXP == 607
#define SFMT_N     5
#define SFMT_M     2
#define SFMT_SL1	15
#define SFMT_SL2	 3
#define SFMT_SR1	13
#define SFMT_SR2	 3
#define SFMT_MASK	  0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x5986f054
#endif

class Random {
 public:
  Random(int seed); 
  void RandomInit(int seed);
  void RandomInitByArray(int const seeds[], int NumSeeds);
  int IRan(int min, int max);
  __m256i Ran16();
  int IRanX(int min, int max);
  uint32_t RanUint32_12();
  uint16_t Ran16_12();
  uint8_t RanUint8_12();
  double Ran();
  uint32_t BRan();
  uint16_t BRan16();
  uint8_t BRan8();
  __m128i BinRan128();
  __m256i BinRan256();
 private:
  void Init2();
  void Generate();
  uint32_t ix;         
  uint32_t LastInterval;
  uint32_t RLimit;     
  __m128i mask;     
  __m128i state[SFMT_N];
};

#endif /* __Random_hpp */ 
