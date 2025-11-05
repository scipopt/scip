/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   rectlu_num.c
 * @brief  rectlu memory functions
 * @author Leon Eifler
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
/*
 * include build configuration flags
 */
#include "scip/config.h"

#include "rectlu/rectlu_factor.h"
#include "rectlu/rectlu.h"

#pragma GCC diagnostic ignored "-Wpedantic"

#ifdef SCIP_WITH_GMP

/* allocates array with size elements of QSnum_type */
QSnum_type* QSnum_AllocArray(int size)
{
   int i;
   QSnum_type* res = (QSnum_type *) malloc(size*sizeof(QSnum_type));
   if (res)
   {
      for (i = 0; i < size; i++)
         mpq_init(res[i]);
   }
   return res;
}

/* frees array ea with size elements of QSnum_type */
void QSnum_FreeArray(QSnum_type* ea, int size)
{
   int i;
   if (ea)
   {
      for (i = 0; i < size; i++)
         mpq_clear(ea[i]);
   }
   free(ea);
}

#endif
