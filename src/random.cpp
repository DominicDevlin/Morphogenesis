/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <iostream>
#include "random.h"
#include "parameter.h"
#include <cmath>
#include <stdint.h>

extern Parameter par;



uint64_t splitmix(uint64_t x) 
{
	uint64_t z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}


uint64_t Seed(int org_num)
{
  // Set the seed according to the local time * prg num
  struct timeb t;
  long seed;
  ftime(&t);
  seed=static_cast<uint64_t>(abs((int)((t.time*t.millitm)%655337)) * org_num);
  return splitmix(seed);
}


// thread safe xooshiro256**

uint64_t rotl(uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}


uint64_t next(uint64_t s[4]) 
{
	const uint64_t result = rotl(s[1] * 5, 7) * 9;

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}


void jump(uint64_t s[4]) 
{
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 1;
	uint64_t s1 = 1;
	uint64_t s2 = 1;
	uint64_t s3 = 1;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			next(s);	
		}
		
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}


double to_double(uint64_t x) 
{
   const union { uint64_t i; double d; } u = {.i = UINT64_C(0x3FF) << 52 | x >> 12 };
   return u.d - 1.0;
}


double RANDOM(uint64_t s[4])
{
  return to_double(next(s));
}



long RandomNumber(long max, uint64_t s[4])
{
  return ((long)(RANDOM(s)*max+1));
}


// OLD RNG BELOW

/*! \param An integer random seed
  \return the random seed
**/
// int Seed(int seed)
// {
//   if (seed < 0) {
// 	  std::cerr << "Randomizing random generator, seed is ";
//     int rseed=Randomize();
//     std::cerr << rseed << "\n";
//     return rseed;
//   } else {
//     int i;
//     idum = -seed;
//     for (i=0; i <100; i++)
//       RANDOM();
//     return seed;
//   }
// }



// static int idum = -1;

// /*! \return A random double between 0 and 1
// **/
// double RANDOM(void)
// /* Knuth's substrative method, see Numerical Recipes */
// {
//   static int inext,inextp;
//   static long ma[56];
//   static int iff=0;
//   long mj,mk;
//   int i,ii,k;

//   if (idum < 0 || iff == 0) {
//     iff=1;
//     mj=MSEED-(idum < 0 ? -idum : idum);
//     mj %= MBIG;
//     ma[55]=mj;
//     mk=1;
//     i=1;
//     do {
//       ii=(21*i) % 55;
//       ma[ii]=mk;
//       mk=mj-mk;
//       if (mk < MZ) mk += MBIG;
//       mj=ma[ii];
//     } while ( ++i <= 54 );
//     k=1;
//     do {
//       i=1;
//       do {
//         ma[i] -= ma[1+(i+30) % 55];
//         if (ma[i] < MZ) ma[i] += MBIG;
//       } while ( ++i <= 55 );
//     } while ( ++k <= 4 );
//     inext=0;
//     inextp=31;
//     idum=1;
//   }
//   if (++inext == 56) inext=1;
//   if (++inextp == 56) inextp=1;
//   mj=ma[inext]-ma[inextp];
//   if (mj < MZ) mj += MBIG;
//   ma[inext]=mj;
//   return mj*FAC;
// }


/*! Returns a random integer value between 1 and 'max'
  \param The maximum value (long)
  \return A random integer (long)
**/


// long RandomNumber(long max)
// {
//    return((long)(RANDOM()*max+1));
// }

/*! Interactively ask for the seed
\param void
\return void
**/
// void AskSeed(void)
// {
//   int seed;
//   printf("Please enter a random seed: ");
//   scanf("%d",&seed);
//   printf("\n");
//   Seed(seed);
// }


/*! Make a random seed based on the local time
\param void
\return void
**/

// int Randomize(void) {
  
//   // Set the seed according to the local time
//   struct timeb t;
//   int seed;

//   ftime(&t);
  
//   seed=abs((int)((t.time*t.millitm)%655337));
//   Seed(seed);
//   fprintf(stderr,"Random seed is %d\n",seed);
//   return seed;
// }


