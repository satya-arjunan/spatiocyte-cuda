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
// based on the GCC source expmed.c and constant unsigned division multiplier
// by Denis Vlasenko 

#ifndef __Optimization_hpp
#define __Optimization_hpp

#include <Common.hpp>
#define HOST_WIDE_INT mol_t 
#define UNSIGNED_WIDE_INT umol_t
#define HOST_DOUBLE_WIDE_INT mol2_t 
#define UNSIGNED_DOUBLE_WIDE_INT umol2_t
/*
#define HOST_WIDE_INT int32_t 
#define UNSIGNED_WIDE_INT uint32_t
#define HOST_DOUBLE_WIDE_INT int64_t 
#define UNSIGNED_DOUBLE_WIDE_INT uint64_t
*/
#define MAX_HOST_WIDE_INT ((UNSIGNED_WIDE_INT) (~(HOST_WIDE_INT)0))
#define HOST_BITS_PER_WIDE_INT (sizeof(HOST_WIDE_INT)*8)
#define HOST_BITS_PER_DOUBLE_INT (2 * HOST_BITS_PER_WIDE_INT)
#define OVERFLOW_SUM_SIGN(a, b, sum) ((~((a) ^ (b)) & ((a) ^ (sum))) < 0)
#define LOWPART(x) \
  ((x) & (((UNSIGNED_WIDE_INT) 1 << (HOST_BITS_PER_WIDE_INT / 2)) - 1))
#define HIGHPART(x) \
  ((UNSIGNED_WIDE_INT) (x) >> HOST_BITS_PER_WIDE_INT / 2)
#define BASE ((UNSIGNED_WIDE_INT) 1 << HOST_BITS_PER_WIDE_INT / 2)
#define EXACT_POWER_OF_2_OR_ZERO_P(x) \
  (((x) & ((x) - (UNSIGNED_WIDE_INT) 1)) == 0)

int floor_log2 (UNSIGNED_WIDE_INT x);
void encode (HOST_WIDE_INT *words, UNSIGNED_WIDE_INT low,
    HOST_WIDE_INT hi);
void decode (HOST_WIDE_INT *words, UNSIGNED_WIDE_INT *low,
    HOST_WIDE_INT *hi);
int neg_double (UNSIGNED_WIDE_INT l1, HOST_WIDE_INT h1,
    UNSIGNED_WIDE_INT *lv, HOST_WIDE_INT *hv);
int add_double (UNSIGNED_WIDE_INT l1, HOST_WIDE_INT h1, UNSIGNED_WIDE_INT l2,
    HOST_WIDE_INT h2, UNSIGNED_WIDE_INT *lv, HOST_WIDE_INT *hv);
int mul_double (UNSIGNED_WIDE_INT l1, HOST_WIDE_INT h1, UNSIGNED_WIDE_INT l2,
    HOST_WIDE_INT h2, UNSIGNED_WIDE_INT *lv, HOST_WIDE_INT *hv);
int div_and_round_double (int uns, UNSIGNED_WIDE_INT lnum_orig,
    HOST_WIDE_INT hnum_orig, UNSIGNED_WIDE_INT lden_orig,
    HOST_WIDE_INT hden_orig, UNSIGNED_WIDE_INT *lquo, HOST_WIDE_INT *hquo,
    UNSIGNED_WIDE_INT *lrem, HOST_WIDE_INT *hrem);
int set_const_division_param(UNSIGNED_WIDE_INT divisor, 
    UNSIGNED_WIDE_INT *multiplier, UNSIGNED_WIDE_INT *nshift);

#endif /* __Optimization_hpp */
