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

#include <Common.hpp>
#include <Random.hpp>

Random::Random(const int seed)
    : LastInterval(0) {
  RandomInit(seed);
}

void Random::RandomInit(int seed) {
  // Re-seed
  uint32_t i;                         // Loop counter
  uint32_t y = seed;                  // Temporary
  uint32_t statesize = SFMT_N*4;      // Size of state vector

  // Fill state vector with random numbers from seed
  ((uint32_t*)state)[0] = y;
  const uint32_t factor = 1812433253U;// Multiplication factor

  for (i = 1; i < statesize; i++) {
     y = factor * (y ^ (y >> 30)) + i;
     ((uint32_t*)state)[i] = y;
  }

  // Further initialization and period certification
  Init2();
}

// Functions used by Random::RandomInitByArray
static uint32_t func1(uint32_t x) {
    return (x ^ (x >> 27)) * 1664525U;
}

static uint32_t func2(uint32_t x) {
    return (x ^ (x >> 27)) * 1566083941U;
}

void Random::RandomInitByArray(int const seeds[], int NumSeeds) {
   // Seed by more than 32 bits
   uint32_t i, j, count, r, lag;

   if (NumSeeds < 0) NumSeeds = 0;

   const uint32_t size = SFMT_N*4; // number of 32-bit integers in state

   // Typecast state to uint32_t *
   uint32_t * sta = (uint32_t*)state;

   if (size >= 623) {
      lag = 11;} 
   else if (size >= 68) {
      lag = 7;}
   else if (size >= 39) {
      lag = 5;}
   else {
      lag = 3;
   }
   const uint32_t mid = (size - lag) / 2;

   if ((uint32_t)NumSeeds + 1 > size) {
      count = (uint32_t)NumSeeds;
   }
   else {
      count = size - 1;
   }
#if 0
   // Original code. Argument to func1 is constant!
   for (i = 0; i < size; i++) sta[i] = 0x8B8B8B8B;
   r = func1(sta[0] ^ sta[mid] ^ sta[size - 1]);
   sta[mid] += r;
   r += NumSeeds;
   sta[mid + lag] += r;
   sta[0] = r;
#else
   // 1. loop: Fill state vector with random numbers from NumSeeds
   const uint32_t factor = 1812433253U;// Multiplication factor
   r = (uint32_t)NumSeeds;
   for (i = 0; i < SFMT_N*4; i++) {
      r = factor * (r ^ (r >> 30)) + i;
      sta[i] = r;
   }

#endif

   // 2. loop: Fill state vector with random numbers from seeds[]
   for (i = 1, j = 0; j < count; j++) {
      r = func1(sta[i] ^ sta[(i + mid) % size] ^ sta[(i + size - 1) % size]);
      sta[(i + mid) % size] += r;
      if (j < (uint32_t)NumSeeds) r += (uint32_t)seeds[j];
      r += i;
      sta[(i + mid + lag) % size] += r;
      sta[i] = r;
      i = (i + 1) % size;
   }

   // 3. loop: Randomize some more
   for (j = 0; j < size; j++) {
      r = func2(sta[i] + sta[(i + mid) % size] + sta[(i + size - 1) % size]);
      sta[(i + mid) % size] ^= r;
      r -= i;
      sta[(i + mid + lag) % size] ^= r;
      sta[i] = r;
      i = (i + 1) % size;
   }
   
   // Further initialization and period certification
   Init2();
}


void Random::Init2() {
   // Various initializations and period certification
   uint32_t i, j, temp;

   // Initialize mask
   static const uint32_t maskinit[4] = {SFMT_MASK};
   mask = _mm_loadu_si128((__m128i*)maskinit);

   // Period certification
   // Define period certification vector
   static const uint32_t parityvec[4] = {SFMT_PARITY};

   // Check if parityvec & state[0] has odd parity
   temp = 0;
   for (i = 0; i < 4; i++) {
      temp ^= parityvec[i] & ((uint32_t*)state)[i];
   }
   for (i = 16; i > 0; i >>= 1) temp ^= temp >> i;
   if (!(temp & 1)) {
      // parity is even. Certification failed
      // Find a nonzero bit in period certification vector
      for (i = 0; i < 4; i++) {
         if (parityvec[i]) {
            for (j = 1; j; j <<= 1) {
               if (parityvec[i] & j) {
                  // Flip the corresponding bit in state[0] to change parity
                  ((uint32_t*)state)[i] ^= j;
                  // Done. Exit i and j loops
                  i = 5;  break;
               }
            }
         }
      }
   }
   // Generate first random numbers and set ix = 0
   Generate();
}


// Subfunction for the sfmt algorithm
static inline __m128i sfmt_recursion(__m128i const &a, __m128i const &b, 
__m128i const &c, __m128i const &d, __m128i const &mask) {
    __m128i a1, b1, c1, d1, z1, z2;
    b1 = _mm_srli_epi32(b, SFMT_SR1);
    a1 = _mm_slli_si128(a, SFMT_SL2);
    c1 = _mm_srli_si128(c, SFMT_SR2);
    d1 = _mm_slli_epi32(d, SFMT_SL1);
    b1 = _mm_and_si128(b1, mask);
    z1 = _mm_xor_si128(a, a1);
    z2 = _mm_xor_si128(b1, d1);
    z1 = _mm_xor_si128(z1, c1);
    z2 = _mm_xor_si128(z1, z2);
    return z2;
}

void Random::Generate() {
   // Fill state array with new random numbers
   int i;
   __m128i r, r1, r2;

   r1 = state[SFMT_N - 2];
   r2 = state[SFMT_N - 1];
   for (i = 0; i < SFMT_N - SFMT_M; i++) {
      r = sfmt_recursion(state[i], state[i + SFMT_M], r1, r2, mask);
      state[i] = r;
      r1 = r2;
      r2 = r;
   }
   for (; i < SFMT_N; i++) {
      r = sfmt_recursion(state[i], state[i + SFMT_M - SFMT_N], r1, r2, mask);
      state[i] = r;
      r1 = r2;
      r2 = r;
   }
   ix = 0;
}

int Random::IRan(int min, int max) {
   // Assume 64 bit integers supported. Use multiply and shift method
   uint32_t interval;                  // Length of interval
   uint64_t longran;                   // Random bits * interval
   uint32_t iran;                      // Longran / 2^32

   interval = (uint32_t)(max - min + 1);
   longran  = (uint64_t)BRan() * interval;
   iran = (uint32_t)(longran >> 32);
   // Convert back to signed and return result
   return (int32_t)iran + min;
}


uint32_t Random::RanUint32_12() {
   return (uint32_t)(((uint64_t)BRan()*12) >> 32);
}

/*
//Doesn't work, not accurate
uint8_t Random::RanUint8_12() {
   return (uint8_t)(((uint16_t)BRan8()*12) >> 8);
}

uint8_t Random::BRan8() {
   // Output 32 random bits
   uint8_t y;

   if (ix >= SFMT_N*16) {
      Generate();
   }
   y = ((uint8_t*)state)[ix++];
   return y;
}
*/

uint16_t Random::BRan16() {
   // Output 32 random bits
   uint16_t y;

   if (ix >= SFMT_N*8) {
      Generate();
   }
   y = ((uint16_t*)state)[ix++];
   return y;
}

uint32_t Random::BRan() {
   // Output 32 random bits
   uint32_t y;

   if (ix >= SFMT_N*4) {
      Generate();
   }
   y = ((uint32_t*)state)[ix++];
   return y;
}

__m128i Random::BinRan128() {
  if (ix >= SFMT_N) {
    Generate();
  }
  return state[ix++];
}

__m256i Random::BinRan256() {
  if (ix >= SFMT_N/2) {
    Generate();
  }
  return ((__m256i*)state)[ix++];
}

uint16_t Random::Ran16_12() {
   return (uint16_t)(((uint32_t)BRan16()*12) >> 16);
}

__m256i Random::Ran16() {
  return _mm256_mulhi_epu16(BinRan256(), _mm256_set1_epi16(12));
}

/*
//VPMULHUW __m256i _mm256_mulhi_epu16 ( __m256i a, __m256i b)
//VPMADDUBSW __m256i _mm256_maddubs_epi16 (__m256i a, __m256i b)
//VPSRLW __m256i _mm_srli_epi16 (__m256i m, int count) (V)PSRLW
__m256i Random::Ran16() {
  //BinRan128(): 8bit*16numbers = 128 bits
  //Cast uint8_t to uint16_t
  __m256i ran(_mm256_cvtepu8_epi16(BinRan128()));
  //Multiply by 12
  ran = _mm256_maddubs_epi16(ran, m256i_12_);
  _mm256_mulhi_epu16(BinRan256(), m256i_12_);
  //Shift right logical by 8 counts
  return _mm256_srli_epi16(ran, 8);
}
*/

/*
//Doesn't work (not accurate):
uint8_t Random::RanUint8_12() {
  uint8_t ran8(BRan8());
  //std::cout << "ran8:";
  //cout_binary(ran8);
  uint16_t ran16((uint16_t)ran8);
  //std::cout << "ran16:";
  //cout_binary(ran16);
  ran16 = ran16*12;
  //std::cout << "ran16*12:";
  //cout_binary(ran16);
  ran16 = ran16 >> 8;
  //std::cout << "ran16 >> 8:";
  //cout_binary(ran16);
  ran8 = (uint8_t)ran16;
  //std::cout << "ran8 final:";
  //cout_binary(ran8);
  return ran8;
}
*/


/*
union256i_uint16 Random::Ran16() {
  //8bit*16numbers = 128 bits
  //Cast uint8_t to uint16_t
  union256i_uint16 y;
  union128i_uint8 x;
  x.x = BinRan128();
  y.x = _mm256_cvtepu8_epi16(x.x);
  //Multiply with 12
  y.x = _mm256_maddubs_epi16(y.x, const12_.x);
  //Shift right logical by 8 counts
  y.x = _mm256_srli_epi16(y.x, 8);
  std::cout << "before: ";
  cout_binary(x.x);
  for(unsigned i(0); i != 16; ++i)
    {
      std::cout << "i:" << i << std::endl;
      uint8_t ran8(x.a[i]);
      std::cout << "ran8:";
      cout_binary(ran8);
      uint16_t ran16((uint16_t)ran8);
      std::cout << "ran16:";
      cout_binary(ran16);
      ran16 = ran16*12;
      std::cout << "ran16*12:";
      cout_binary(ran16);
      ran16 = ran16 >> 8;
      std::cout << "ran16 >> 8:";
      cout_binary(ran16);
      ran8 = (uint8_t)ran16;
      std::cout << "ran8 final: " << (uint16_t)ran8 << " ";
      cout_binary(ran8);
    }
  std::cout << "after" << std::endl;
  return y;
}
*/

/*
union256i_uint16 Random::Ran16() {
  union256i_uint16 y;
  y.x = _mm256_mulhi_epu16(BinRan256(), const12_.x);
  return y;
}
*/

/*
union256i_uint16 Random::Ran16() {
  union256i_uint16 y;
  union256i_uint16 x;
  y.x = BinRan256();
  x.x = y.x;
  y.x = _mm256_mulhi_epu16 (y.x, const12_.x);
  std::cout << "before: ";
  cout_binary(x.x);
  for(unsigned i(0); i != 16; ++i)
    {
      std::cout << "i:" << i << std::endl;
      uint16_t ran16(x.a[i]);
      std::cout << "ran16:";
      cout_binary(ran16);
      uint32_t ran32((uint32_t)ran16);
      std::cout << "ran32:";
      cout_binary(ran32);
      ran32 = ran32*12;
      std::cout << "ran32*12:";
      cout_binary(ran32);
      ran32 = ran32 >> 16;
      std::cout << "ran32 >> 16:";
      cout_binary(ran32);
      ran16 = (uint16_t)ran32;
      std::cout << "ran16 final:" << ran16 << " ";
      cout_binary(ran16);
      //uint16_t a((uint16_t)(((uint32_t)(x.x[i])*12) >> 16));
    }
  std::cout << "after" << std::endl;
  return y;
}
*/

int  Random::IRanX (int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) {
         return min;                   // max == min. Only one possible value
      }
      else {
         return 0x80000000;            // max < min. Error output
      }
   }
   // Assume 64 bit integers supported. Use multiply and shift method
   uint32_t interval;                  // Length of interval
   uint64_t longran;                   // Random bits * interval
   uint32_t iran;                      // Longran / 2^32
   uint32_t remainder;                 // Longran % 2^32

   interval = (uint32_t)(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder = 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then.
      RLimit = (uint32_t)(((uint64_t)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64_t)BRan() * interval;
      iran = (uint32_t)(longran >> 32);
      remainder = (uint32_t)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int32_t)iran + min;
}

double Random::Ran() {
   // Output random floating point number
   if (ix >= SFMT_N*4-1) {
      // Make sure we have at least two 32-bit numbers
      Generate();
   }
   uint64_t r = *(uint64_t*)((uint32_t*)state+ix);
   ix += 2;
   // 53 bits resolution:
   // return (int64_t)(r >> 11) * (1./(67108864.0*134217728.0)); 
   // (r >> 11)*2^(-53)
   // 52 bits resolution for compatibility with assembly version:
   // (r >> 12)*2^(-52)
   return (int64_t)(r >> 12) * (1./(67108864.0*67108864.0));
}
