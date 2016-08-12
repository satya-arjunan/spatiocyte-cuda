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

#include <string.h>
#include <Optimization.hpp>

int floor_log2 (UNSIGNED_WIDE_INT x) {
  int t(0);
  if (x == 0) {
    return -1;
  }
  if (HOST_BITS_PER_WIDE_INT > 64) {
    if (x >= (UNSIGNED_WIDE_INT) 1 << (t + 64)) {
      t += 64;
    }
  }
  if (HOST_BITS_PER_WIDE_INT > 32) {
    if (x >= ((UNSIGNED_WIDE_INT) 1) << (t + 32)) {
      t += 32;
    }
  }
  if (x >= ((UNSIGNED_WIDE_INT) 1) << (t + 16)) {
    t += 16;
  }
  if (x >= ((UNSIGNED_WIDE_INT) 1) << (t + 8)) {
    t += 8;
  }
  if (x >= ((UNSIGNED_WIDE_INT) 1) << (t + 4)) {
    t += 4;
  }
  if (x >= ((UNSIGNED_WIDE_INT) 1) << (t + 2)) {
    t += 2;
  }
  if (x >= ((UNSIGNED_WIDE_INT) 1) << (t + 1)) {
    t += 1;
  }
  return t;
}

void encode (HOST_WIDE_INT *words, UNSIGNED_WIDE_INT low,
    HOST_WIDE_INT hi) {
  words[0] = LOWPART (low);
  words[1] = HIGHPART (low);
  words[2] = LOWPART (hi);
  words[3] = HIGHPART (hi);
}

void decode (HOST_WIDE_INT *words, UNSIGNED_WIDE_INT *low,
    HOST_WIDE_INT *hi) {
  *low = words[0] + words[1] * BASE;
  *hi = words[2] + words[3] * BASE;
}

int neg_double (UNSIGNED_WIDE_INT l1, HOST_WIDE_INT h1,
    UNSIGNED_WIDE_INT *lv, HOST_WIDE_INT *hv) {
  if (l1 == 0) {
    *lv = 0;
    *hv = - (UNSIGNED_WIDE_INT) h1;
    return (*hv & h1) < 0;
  } else {
    *lv = -l1;
    *hv = ~h1;
    return 0;
  }
}

int add_double (UNSIGNED_WIDE_INT l1, HOST_WIDE_INT h1, UNSIGNED_WIDE_INT l2,
    HOST_WIDE_INT h2, UNSIGNED_WIDE_INT *lv, HOST_WIDE_INT *hv) {
  UNSIGNED_WIDE_INT l;
  HOST_WIDE_INT h;
  l = l1 + l2;
  h = h1 + h2 + (l < l1);
  *lv = l;
  *hv = h;
  return OVERFLOW_SUM_SIGN(h1, h2, h);
}

int mul_double (UNSIGNED_WIDE_INT l1, HOST_WIDE_INT h1, UNSIGNED_WIDE_INT l2,
    HOST_WIDE_INT h2, UNSIGNED_WIDE_INT *lv, HOST_WIDE_INT *hv) {
  HOST_WIDE_INT arg1[4];
  HOST_WIDE_INT arg2[4];
  HOST_WIDE_INT prod[4 * 2];
  UNSIGNED_WIDE_INT carry;
  int i, j, k;
  UNSIGNED_WIDE_INT toplow, neglow;
  HOST_WIDE_INT tophigh, neghigh;
  encode (arg1, l1, h1);
  encode (arg2, l2, h2);
  memset (prod, 0, sizeof prod);
  for (i = 0; i < 4; i++) {
    carry = 0;
    for (j = 0; j < 4; j++) {
      k = i + j;
      /* This product is <= 0xFFFE0001, the sum <= 0xFFFF0000.  */
      carry += arg1[i] * arg2[j];
      /* Since prod[p] < 0xFFFF, this sum <= 0xFFFFFFFF.  */
      carry += prod[k];
      prod[k] = LOWPART (carry);
      carry = HIGHPART (carry);
    }
    prod[i + 4] = carry;
  }
  decode (prod, lv, hv);	/* This ignores prod[4] through prod[4*2-1] */
  /* Check for overflow by calculating the top half of the answer in full;
     it should agree with the low half's sign bit.  */
  decode (prod + 4, &toplow, &tophigh);
  if (h1 < 0) {
    neg_double (l2, h2, &neglow, &neghigh);
    add_double (neglow, neghigh, toplow, tophigh, &toplow, &tophigh);
  }
  if (h2 < 0) {
    neg_double (l1, h1, &neglow, &neghigh);
    add_double (neglow, neghigh, toplow, tophigh, &toplow, &tophigh);
  }
  return (*hv < 0 ? ~(toplow & tophigh) : toplow | tophigh) != 0;
}

int div_and_round_double (/* enum tree_code code,*/ int uns,
		      UNSIGNED_WIDE_INT lnum_orig, /* num == numerator == dividend */
		      HOST_WIDE_INT hnum_orig,
		      UNSIGNED_WIDE_INT lden_orig, /* den == denominator == divisor */
		      HOST_WIDE_INT hden_orig,
		      UNSIGNED_WIDE_INT *lquo,
		      HOST_WIDE_INT *hquo, UNSIGNED_WIDE_INT *lrem,
		      HOST_WIDE_INT *hrem)
{
  int quo_neg = 0;
  HOST_WIDE_INT num[4 + 1];	/* extra element for scaling.  */
  HOST_WIDE_INT den[4], quo[4];
  int i, j;
  UNSIGNED_WIDE_INT work;
  UNSIGNED_WIDE_INT carry = 0;
  UNSIGNED_WIDE_INT lnum = lnum_orig;
  HOST_WIDE_INT hnum = hnum_orig;
  UNSIGNED_WIDE_INT lden = lden_orig;
  HOST_WIDE_INT hden = hden_orig;
  int overflow = 0;
  if (hden == 0 && lden == 0) {
    overflow = 1, lden = 1;
  }
  /* Calculate quotient sign and convert operands to unsigned.  */
  if (!uns) {
    if (hnum < 0) {
      quo_neg = ~ quo_neg;
      /* (minimum integer) / (-1) is the only overflow case.  */
      if (neg_double (lnum, hnum, &lnum, &hnum) &&
          ((HOST_WIDE_INT) lden & hden) == -1) {
        overflow = 1;
      }
    }
    if (hden < 0) {
      quo_neg = ~ quo_neg;
      neg_double (lden, hden, &lden, &hden);
    }
  }
  if (hnum == 0 && hden == 0) {
    /* single precision */
    *hquo = *hrem = 0;
    /* This unsigned division rounds toward zero.  */
    *lquo = lnum / lden;
    goto finish_up;
  }
  if (hnum == 0) {
    /* trivial case: dividend < divisor */
    /* hden != 0 already checked.  */
    *hquo = *lquo = 0;
    *hrem = hnum;
    *lrem = lnum;
    goto finish_up;
  }
  memset (quo, 0, sizeof quo);
  memset (num, 0, sizeof num);	/* to zero 9th element */
  memset (den, 0, sizeof den);
  encode (num, lnum, hnum);
  encode (den, lden, hden);
  /* Special code for when the divisor < BASE.  */
  if (hden == 0 && lden < (UNSIGNED_WIDE_INT) BASE) {
    /* hnum != 0 already checked.  */
    for (i = 4 - 1; i >= 0; i--) {
      work = num[i] + carry * BASE;
      quo[i] = work / lden;
      carry = work % lden;
    }
  } else {
    /* Full double precision division,
     * with thanks to Don Knuth's "Seminumerical Algorithms".  */
    int num_hi_sig, den_hi_sig;
    UNSIGNED_WIDE_INT quo_est, scale;
    /* Find the highest nonzero divisor digit.  */
    for (i = 4 - 1;; i--) {
      if (den[i] != 0) {
        den_hi_sig = i;
        break;
      }
    } 
    /* Ensure that the first digit of the divisor is at least BASE/2.
     * This is required by the quotient digit estimation algorithm.  */ 
    scale = BASE / (den[den_hi_sig] + 1);
    if (scale > 1) {
      /* scale divisor and dividend */
      carry = 0;
      for (i = 0; i <= 4 - 1; i++) {
        work = (num[i] * scale) + carry;
        num[i] = LOWPART (work);
	      carry = HIGHPART (work);
      } 
      num[4] = carry;
      carry = 0;
      for (i = 0; i <= 4 - 1; i++) {
	      work = (den[i] * scale) + carry;
	      den[i] = LOWPART (work);
	      carry = HIGHPART (work);
	      if (den[i] != 0) {
          den_hi_sig = i;
        }
      }
    } 
    num_hi_sig = 4;
    /* Main loop */
    for (i = num_hi_sig - den_hi_sig - 1; i >= 0; i--) {
      /* Guess the next quotient digit, quo_est, by dividing the first
       * two remaining dividend digits by the high order quotient digit.
       * quo_est is never low and is at most 2 high.  */
      UNSIGNED_WIDE_INT tmp; 
      num_hi_sig = i + den_hi_sig + 1;
      work = num[num_hi_sig] * BASE + num[num_hi_sig - 1];
      if (num[num_hi_sig] != den[den_hi_sig]) {
        quo_est = work / den[den_hi_sig];
      } else {
        quo_est = BASE - 1;
      } 
      /* Refine quo_est so it's usually correct, and at most one high.  */
      tmp = work - quo_est * den[den_hi_sig];
      if (tmp < BASE && (den[den_hi_sig - 1] *
                         quo_est > (tmp * BASE + num[num_hi_sig - 2]))) {
        quo_est--;
      } 
      /* Try QUO_EST as the quotient digit, by multiplying the
       * divisor by QUO_EST and subtracting from the remaining dividend.
       * Keep in mind that QUO_EST is the I - 1st digit.  */ 
      carry = 0;
      for (j = 0; j <= den_hi_sig; j++) {
        work = quo_est * den[j] + carry;
	      carry = HIGHPART (work);
	      work = num[i + j] - LOWPART (work);
	      num[i + j] = LOWPART (work);
	      carry += HIGHPART (work) != 0;
	    } 
      /* If quo_est was high by one, then num[i] went negative and
       * we need to correct things.  */
      if (num[num_hi_sig] < (HOST_WIDE_INT) carry) {
        quo_est--;
	      carry = 0;		/* add divisor back in */
	      for (j = 0; j <= den_hi_sig; j++) {
          work = num[i + j] + den[j] + carry;
          carry = HIGHPART (work);
          num[i + j] = LOWPART (work);
        } 
        num [num_hi_sig] += carry;
      } 
      /* Store the quotient digit.  */
      quo[i] = quo_est;
    }
  } 
  decode (quo, lquo, hquo);
finish_up:
  /* If result is negative, make it so.  */
  if (quo_neg) {
    neg_double (*lquo, *hquo, lquo, hquo);
  }
  /* Compute trial remainder:  rem = num - (quo * den)  */
  mul_double (*lquo, *hquo, lden_orig, hden_orig, lrem, hrem);
  neg_double (*lrem, *hrem, lrem, hrem);
  add_double (lnum_orig, hnum_orig, *lrem, *hrem, lrem, hrem);
  return overflow;
}

int set_const_division_param(UNSIGNED_WIDE_INT B,
    UNSIGNED_WIDE_INT *multiplier_ptr, UNSIGNED_WIDE_INT *shift_ptr) {
  UNSIGNED_WIDE_INT max_A = MAX_HOST_WIDE_INT;
  UNSIGNED_WIDE_INT K, last_K, d_LB, bad_A, udummy;
  HOST_WIDE_INT high_A, dummy;
  int lgdn; 
  /* L/B should fit into uint, so max L is floor_log2(B) << 
   * HOST_BITS_PER_WIDE_INT
   * L is impicitly = 1 << (lgdn+HOST_BITS_PER_WIDE_INT) */
  lgdn = floor_log2(B); 
  /* K = L/B + 1 */
  div_and_round_double(1, 0,1<<lgdn, B,(HOST_WIDE_INT)0, &K,&dummy, &udummy,
                       &dummy);
  K++; 
  /* d_LB = ((L/B) * B + B - L) */
  /* since (HOST_WIDE_INT)L == 0, subtracting it can be optimized out */
  d_LB = K * B;
  /* bad_A = L / d_LB */
  div_and_round_double(1, 0,1<<lgdn, d_LB,(HOST_WIDE_INT)0, &bad_A,&high_A,
                       &udummy,&dummy); 
  bad_A = (bad_A / B) * B; /* ... + B-1 later */
  if(high_A || bad_A > MAX_HOST_WIDE_INT - (B-1)) {
    high_A = 1; /* we aren't interested much in true value if high_A!=0 */
  } else {
    bad_A += B-1;
  }
  if (!high_A && bad_A <= max_A) {
    /* K doesn't work for all required A's.
     * However, there always exist
     * (potentially HOST_BITS_PER_WIDE_INT+1 bit) K which works
     * for all A's and it equals K*2-1.
     * We return lower HOST_BITS_PER_WIDE_INT bits in *multiplier_ptr. */
    *multiplier_ptr = (K-1)*2 + 1; /* = K*2-1 but overflow-safe */
    *shift_ptr = lgdn+1;
    /* if K > 0x800..000 then K*2-1 overflows into 32th bit,
     * and we need to return 33-bit value.
     * this code returns 1/0 = K*2-1 is 33bit / K*2-1 is 32bit or less */
    return K > ((UNSIGNED_WIDE_INT)1 << (HOST_BITS_PER_WIDE_INT-1));
  }
  /* There is not much point in multiplying by even number
  ** and then shifting right. Reduce K & L to avoid it: */
  while (!(K & 1) && lgdn) {
    lgdn--, K >>= 1;
  }
  /* K is good, but maybe we can reduce it? */
  while(1) {
    last_K = K;
    if (--lgdn < 0) {
      break;
    }
    /* K = L/B + 1 */
    div_and_round_double(1, 0, 1<<lgdn, B, (HOST_WIDE_INT) 0, &K, &dummy,
                         &udummy, &dummy);
    K++;
    /* d_LB = ((L/B) * B + B - L) */
    /* since (HOST_WIDE_INT)L == 0, subtracting it can be optimized out */
    d_LB = K * B;
    /* bad_A = L / d_LB */
    div_and_round_double(1, 0, 1<<lgdn, d_LB, (HOST_WIDE_INT) 0, &bad_A,
                         &high_A, &udummy, &dummy);
    bad_A = (bad_A / B) * B;
    if (high_A || bad_A > MAX_HOST_WIDE_INT - (B-1) ) {
      high_A = 1; /* we aren't interested much in true value if high_A!=0 */
    } else {
      bad_A += B-1;
    }
    if (!high_A && bad_A <= max_A) {
      break;
    }
  }
  /* Report previous parameters because the current ones are bad */
  *multiplier_ptr = last_K;
  *shift_ptr = lgdn+1;
  return 0; /* 32-bit or less */
}

void test() { 
  UNSIGNED_WIDE_INT type, multiplier, nshift;
  UNSIGNED_WIDE_INT quotient(0); // result
  std::cout << "MAX_HOST_WIDE_INT:" << MAX_HOST_WIDE_INT << std::endl;
  for (UNSIGNED_WIDE_INT divisor(1); divisor <= MAX_HOST_WIDE_INT; ++divisor) {
    std::cout << "divisor:" << divisor << std::endl;
    //n = numerator
    for (UNSIGNED_WIDE_INT numerator(0); numerator <= MAX_HOST_WIDE_INT;
         ++numerator) {
      if (EXACT_POWER_OF_2_OR_ZERO_P(divisor)) {
        type = 2;
        multiplier = 1;
		    nshift = floor_log2(divisor);
        quotient = UNSIGNED_DOUBLE_WIDE_INT(multiplier)*numerator;
        quotient = quotient >> nshift;
      } else {
        type = set_const_division_param(divisor, &multiplier, &nshift); 
        if (type) {
          UNSIGNED_WIDE_INT tmp(0);
          --nshift;
          quotient = (UNSIGNED_DOUBLE_WIDE_INT(multiplier)*numerator) >> 
            HOST_BITS_PER_WIDE_INT;
          tmp = ((numerator-quotient) >> 1) + quotient;
          quotient = tmp >> nshift;
        } else {
          quotient = (UNSIGNED_DOUBLE_WIDE_INT(multiplier)*numerator) >>
            HOST_BITS_PER_WIDE_INT;
          quotient = quotient >> nshift;
        }
      }
      if (quotient != numerator/divisor || 
          numerator%divisor != numerator-divisor*quotient) {
        std::cout << "divisor:" << divisor << " numerator:" <<
          numerator << std::endl;
        std::cout << "host int bits:" << HOST_BITS_PER_WIDE_INT << std::endl;
        std::cout << "type:" << type << " multiplier:" << multiplier <<
          " nshift:" << nshift << std::endl;
        std::cout << "actual quotient:" << numerator/divisor << " modulus:" <<
          numerator%divisor << std::endl;
        std::cout << "computed quotient:" << quotient << " modulus:" <<
          numerator-divisor*quotient << std::endl;
        exit(0);
      }
    }
  }
}

void print_multiplier(UNSIGNED_WIDE_INT divisor) { 
  UNSIGNED_WIDE_INT type, multiplier, nshift;
  if (EXACT_POWER_OF_2_OR_ZERO_P(divisor)) {
    type = 2;
    multiplier = 1;
    nshift = floor_log2(divisor);
    /*
    quotient = UNSIGNED_DOUBLE_WIDE_INT(multiplier)*numerator;
    quotient = quotient >> nshift;
    */
  } else {
    type = set_const_division_param(divisor, &multiplier, &nshift); 
    if (type) {
      --nshift;
      /*
      UNSIGNED_WIDE_INT tmp(0);
      quotient = (UNSIGNED_DOUBLE_WIDE_INT(multiplier)*numerator) >> 
        HOST_BITS_PER_WIDE_INT;
      tmp = ((numerator-quotient) >> 1) + quotient;
      quotient = tmp >> nshift;
      */
    } else {
      /*
      quotient = (UNSIGNED_DOUBLE_WIDE_INT(multiplier)*numerator) >>
        HOST_BITS_PER_WIDE_INT;
      quotient = quotient >> nshift;
      */
    }
  }
  std::cout << "divisor:" << divisor << std::endl;
  std::cout << "host int bits:" << HOST_BITS_PER_WIDE_INT << std::endl;
  std::cout << "type:" << type << " multiplier:" << multiplier <<
    " nshift:" << nshift << std::endl;
}

/*
int main() { 
  //test();
  print_multiplier(1440);
  return 0;
}
*/


