/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bitencode.c
 * @brief  packing single and dual bit values
 * @author Thorsten Koch
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BITENCODE_H__
#define __BITENCODE_H__

typedef int SINGLEPACKET;                         /**< storing single bits in packed format */
#define SINGLEPACKETSIZE (sizeof(SINGLEPACKET)*8) /**< each entry needs one bit of information */
typedef int DUALPACKET;                           /**< storing bit pairs in packed format */
#define DUALPACKETSIZE   (sizeof(DUALPACKET)*4)   /**< each entry needs two bits of information */


extern
void SCIPencodeSingleBit(               /**< encode a single bit vector into packed format */
   const int*       inp,                /**< unpacked input vector */
   SINGLEPACKET*    out,                /**< buffer to store the packed vector */
   int              count               /**< number of elements */
   );

extern
void SCIPdecodeSingleBit(               /**< decode a packed single bit vector into unpacked format */
   const SINGLEPACKET* inp,             /**< packed input vector */
   int*             out,                /**< buffer to store unpacked vector */
   int              count               /**< number of elements */
   );

extern
void SCIPencodeDualBit(                 /**< encode a dual bit vector into packed format */
   const int*       inp,                /**< unpacked input vector */
   DUALPACKET*      out,                /**< buffer to store the packed vector */
   int              count               /**< number of elements */
   );

extern
void SCIPdecodeDualBit(                 /**< decode a packed dual bit vector into unpacked format */
   const DUALPACKET* inp,               /**< packed input vector */
   int*             out,                /**< buffer to store unpacked vector */
   int              count               /**< number of elements */
   );


#endif
