/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: bitencode.c,v 1.8 2005/01/18 09:26:41 bzfpfend Exp $"

/**@file   bitencode.c
 * @brief  packing single and dual bit values
 * @author Thorsten Koch
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "bitencode.h"


/** encode a single bit vector into packed format */
void SCIPencodeSingleBit(
   const int*       inp,                /**< unpacked input vector */
   SINGLEPACKET*    out,                /**< buffer to store the packed vector */
   int              count               /**< number of elements */
   )
{
   static const SINGLEPACKET mask[SINGLEPACKETSIZE][2] = {   /* if the packetsize changes, the mask has to be updated */
      {0x00000000, 0x00000001},
      {0x00000000, 0x00000002},
      {0x00000000, 0x00000004},
      {0x00000000, 0x00000008},
      {0x00000000, 0x00000010},
      {0x00000000, 0x00000020},
      {0x00000000, 0x00000040},
      {0x00000000, 0x00000080},
      {0x00000000, 0x00000100},
      {0x00000000, 0x00000200},
      {0x00000000, 0x00000400},
      {0x00000000, 0x00000800},
      {0x00000000, 0x00001000},
      {0x00000000, 0x00002000},
      {0x00000000, 0x00004000},
      {0x00000000, 0x00008000},
      {0x00000000, 0x00010000},
      {0x00000000, 0x00020000},
      {0x00000000, 0x00040000},
      {0x00000000, 0x00080000},
      {0x00000000, 0x00100000},
      {0x00000000, 0x00200000},
      {0x00000000, 0x00400000},
      {0x00000000, 0x00800000},
      {0x00000000, 0x01000000},
      {0x00000000, 0x02000000},
      {0x00000000, 0x04000000},
      {0x00000000, 0x08000000},
      {0x00000000, 0x10000000},
      {0x00000000, 0x20000000},
      {0x00000000, 0x40000000},
      {0x00000000, 0x80000000}
   };
   int i;
   int rest = count % SINGLEPACKETSIZE;
   int nfull = count - rest;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(SINGLEPACKETSIZE == 32);

   for( i = 0; i < nfull; i += SINGLEPACKETSIZE )
   {
#ifndef NDEBUG
      unsigned int j;
      for( j = 0; j < SINGLEPACKETSIZE; ++j )
         assert(0 <= inp[j] && inp[j] <= 1);
#endif
      *out++ =
	 mask[0][inp[0]] | mask[1][inp[1]] | mask[2][inp[2]] | mask[3][inp[3]]
	 | mask[4][inp[4]] | mask[5][inp[5]] | mask[6][inp[6]] | mask[7][inp[7]]
         | mask[8][inp[8]] | mask[9][inp[9]] | mask[10][inp[10]] | mask[11][inp[11]]
         | mask[12][inp[12]] | mask[13][inp[13]] | mask[14][inp[14]] | mask[15][inp[15]]
         | mask[16][inp[16]] | mask[17][inp[17]] | mask[18][inp[18]] | mask[19][inp[19]]
         | mask[20][inp[20]] | mask[21][inp[21]] | mask[22][inp[22]] | mask[23][inp[23]]
         | mask[24][inp[24]] | mask[25][inp[25]] | mask[26][inp[26]] | mask[27][inp[27]]
         | mask[28][inp[28]] | mask[29][inp[29]] | mask[30][inp[30]] | mask[31][inp[31]];
      inp += SINGLEPACKETSIZE;
   }

   if( rest > 0 )
   {
      SINGLEPACKET  m = (SINGLEPACKET) 0u;

      for( i = 0; i < rest; i++ )
	 m |= mask[i][inp[i]];
      *out = m;
   }
}

/** decode a packed single bit vector into unpacked format */
void SCIPdecodeSingleBit(
   const SINGLEPACKET* inp,             /**< packed input vector */
   int*             out,                /**< buffer to store unpacked vector */
   int              count               /**< number of elements */
   )
{
   SINGLEPACKET m;
   int rest = count % SINGLEPACKETSIZE;
   int nfull = count - rest;
   int i;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(SINGLEPACKETSIZE == 32);

   for( i = 0; i < nfull; i += SINGLEPACKETSIZE )
   {
      m = *inp++;

      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      m >>= 1;
      *out++ = m & 1;
      assert(m >> 1 == 0);
   }

   if( rest > 0 )
   {
      m = *inp;
      for( i = 0; i < rest; i++ )
      {
	 *out++ = m & 1;
	 m >>= 1;
      }
   }
}

/** encode a dual bit vector into packed format */
void SCIPencodeDualBit(
   const int*       inp,                /**< unpacked input vector */
   DUALPACKET*      out,                /**< buffer to store the packed vector */
   int              count               /**< number of elements */
   )
{
   static const DUALPACKET mask[DUALPACKETSIZE][4] = {   /* if the packetsize changes, the mask has to be updated */
      {0x00000000, 0x00000001, 0x00000002, 0x00000003},
      {0x00000000, 0x00000004, 0x00000008, 0x0000000C},
      {0x00000000, 0x00000010, 0x00000020, 0x00000030},
      {0x00000000, 0x00000040, 0x00000080, 0x000000C0},
      {0x00000000, 0x00000100, 0x00000200, 0x00000300},
      {0x00000000, 0x00000400, 0x00000800, 0x00000C00},
      {0x00000000, 0x00001000, 0x00002000, 0x00003000},
      {0x00000000, 0x00004000, 0x00008000, 0x0000C000},
      {0x00000000, 0x00010000, 0x00020000, 0x00030000},
      {0x00000000, 0x00040000, 0x00080000, 0x000C0000},
      {0x00000000, 0x00100000, 0x00200000, 0x00300000},
      {0x00000000, 0x00400000, 0x00800000, 0x00C00000},
      {0x00000000, 0x01000000, 0x02000000, 0x03000000},
      {0x00000000, 0x04000000, 0x08000000, 0x0C000000},
      {0x00000000, 0x10000000, 0x20000000, 0x30000000},
      {0x00000000, 0x40000000, 0x80000000, 0xC0000000}
   };
   int i;
   int rest = count % DUALPACKETSIZE;
   int nfull = count - rest;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(DUALPACKETSIZE == 16);

   for( i = 0; i < nfull; i += DUALPACKETSIZE, inp += DUALPACKETSIZE )
   {
#ifndef NDEBUG
      unsigned int j;
      for( j = 0; j < DUALPACKETSIZE; ++j )
         assert(0 <= inp[j] && inp[j] <= 3);
#endif
      *out++ =
	 mask[0][inp[0]] | mask[1][inp[1]] | mask[2][inp[2]] | mask[3][inp[3]]
	 | mask[4][inp[4]] | mask[5][inp[5]] | mask[6][inp[6]]
         | mask[7][inp[7]] | mask[8][inp[8]] | mask[9][inp[9]] 
         | mask[10][inp[10]] | mask[11][inp[11]] | mask[12][inp[12]] 
         | mask[13][inp[13]] | mask[14][inp[14]] | mask[15][inp[15]];
   }

   if( rest > 0 )
   {
      DUALPACKET  m = (DUALPACKET) 0u;

      for( i = 0; i < rest; i++ )
	 m |= mask[i][inp[i]];
      *out = m;
   }
}

/** decode a packed dual bit vector into unpacked format */
void SCIPdecodeDualBit(
   const DUALPACKET* inp,               /**< packed input vector */
   int*             out,                /**< buffer to store unpacked vector */
   int              count               /**< number of elements */
   )
{
   DUALPACKET m;
   int rest = count % DUALPACKETSIZE;
   int nfull = count - rest;
   int i;

   assert(inp != NULL || count == 0);
   assert(out != NULL || count == 0);
   assert(count >= 0);
   assert(DUALPACKETSIZE == 16);

   for( i = 0; i < nfull; i += DUALPACKETSIZE )
   {
      m = *inp++;

      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      m >>= 2;
      *out++ = m & 3;
      assert(m >> 2 == 0);
   }

   if( rest > 0 )
   {
      m = *inp;
      for( i = 0; i < rest; i++ )
      {
	 *out++ = m & 3;
	 m >>= 2;
      }
   }
}

