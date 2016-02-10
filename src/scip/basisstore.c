/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   basisstore.c
 * @brief  methods for storing separated cuts
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/scip.h"
#include "scip/basisstore.h"
#include "scip/struct_basisstore.h"


/**
 * Methods for SCIP_BASESSTORE
 */

/** initialize the basis storage */
static
SCIP_RETCODE initializeBasisstore(
   SCIP_BASISSTORE*      basisstore,
   BMS_BLKMEM*           blkmem
   )
{
   assert(basisstore != NULL);
   assert(basisstore->basessize == 0);
   assert(blkmem != NULL);

   SCIPdebugMessage("initialize basisstore\n");

   SCIP_ALLOC( BMSallocMemoryArray(&basisstore->bases, 1) );
   basisstore->basessize = 1;

   return SCIP_OKAY;
}

/** ensure that the store is large enoug */
static
SCIP_RETCODE checkMem(
   SCIP_BASISSTORE*      basisstore,
   BMS_BLKMEM*           blkmem
   )
{
   assert(basisstore != NULL);
   assert(blkmem != NULL);

   if( basisstore->nbases == basisstore->basessize )
   {
      SCIP_ALLOC( BMSreallocMemoryArray(&basisstore->bases, 2*basisstore->basessize) );
      basisstore->basessize *= 2;
   }

   return SCIP_OKAY;
}

/** creates basis storage */
SCIP_RETCODE SCIPbasisstoreCreate(
   SCIP_BASISSTORE**     basisstore           /**< pointer to store basis storage */
   )
{
   assert(basisstore != NULL);

   SCIP_ALLOC( BMSallocMemory(basisstore) );

   (*basisstore)->bases = NULL;
   (*basisstore)->basessize = 0;
   (*basisstore)->nbases = 0;

   return SCIP_OKAY;
}

/** frees basis storage */
SCIP_RETCODE SCIPbasisstoreFree(
   SCIP_BASISSTORE**     basisstore,           /**< pointer to store basis storage */
   BMS_BLKMEM*           blkmem
   )
{
   assert(basisstore != NULL);
   assert(*basisstore != NULL);

   if( *basisstore == NULL)
      return SCIP_OKAY;

   while( (*basisstore)->nbases > 0 )
   {
      int b;

      b = (*basisstore)->nbases-1;

      SCIPdebugMessage("free basis %d\n", b);

      BMSfreeMemoryArray(&(*basisstore)->bases[b]->cstat);
      BMSfreeMemoryArray(&(*basisstore)->bases[b]->vstat);
      BMSfreeMemoryArray(&(*basisstore)->bases[b]->conss);
      BMSfreeMemoryArray(&(*basisstore)->bases[b]->vars);
      BMSfreeMemoryArrayNull(&(*basisstore)->bases[b]);
      (*basisstore)->nbases--;
   }

   if( (*basisstore)->bases != NULL )
   {
      BMSfreeMemoryArray(&(*basisstore)->bases);
   }

   BMSfreeMemoryNull(&(*basisstore));

   return SCIP_OKAY;
}

/** add a basis to the storage */
SCIP_RETCODE SCIPbasisstoreAddBasis(
   SCIP_BASISSTORE*      basisstore,
   BMS_BLKMEM*           blkmem,
   SCIP_VAR**            vars,
   SCIP_CONS**           conss,
   int*                  vstat,
   int*                  cstat,
   int                   nvars,
   int                   nconss
   )
{
   assert(basisstore != NULL);
   assert(blkmem != NULL);
   assert(vars != NULL);
   assert(conss != NULL);
   assert(vstat != NULL);
   assert(cstat != NULL);

   /* initialize the storage */
   if( basisstore->basessize == 0 )
   {
      SCIP_CALL( initializeBasisstore(basisstore, blkmem) );
   }

   /* ensure that the storage is large enough */
   if( basisstore->nbases == basisstore->basessize )
   {
      SCIP_CALL( checkMem(basisstore, blkmem) );
   }

   SCIPdebugMessage("add basis: pos %d, nvars %d, nconss %d\n", basisstore->nbases, nvars, nconss);

   /* store the basis */
   SCIP_ALLOC( BMSallocMemory(&basisstore->bases[basisstore->nbases]) );
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->vars, vars, nvars);
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->conss, conss, nconss);
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->vstat, vstat, nvars);
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->cstat, cstat, nconss);
   basisstore->bases[basisstore->nbases]->nvars = nvars;
   basisstore->bases[basisstore->nbases]->nconss = nconss;
   basisstore->nbases++;

   return SCIP_OKAY;
}

/** return the number of stored basis */
int SCIPbasisstoreGetNBasis(
   SCIP_BASISSTORE*      basisstore
   )
{
   assert(basisstore != NULL);

   return basisstore->nbases;
}

/** returns the stored basis */
SCIP_BASIS* SCIPbasestoreGetBasis(
   SCIP_BASISSTORE*      basestore
   )
{
   assert(basestore != NULL);
   assert(basestore->nbases == 1); /* tmp: need to be changed if we store more than one basis */

   return basestore->bases[0];
}

/** copy the basis */
SCIP_RETCODE SCIPbasisstoreCopy(
   SCIP_BASISSTORE*      basisstore,
   SCIP_BASISSTORE*      cpybasisstore,
   BMS_BLKMEM*           blkmem
   )
{
   int b;

   assert(basisstore != NULL);
   assert(cpybasisstore != NULL);

   SCIPdebugMessage("copy basis %p to %p\n", (void*) basisstore, (void*) cpybasisstore);

   for(b = 0; b < basisstore->nbases; b++)
   {
      assert(basisstore->bases[b]->vstat != NULL);
      assert(basisstore->bases[b]->cstat != NULL);

      SCIP_CALL( SCIPbasisstoreAddBasis(cpybasisstore, blkmem,
            basisstore->bases[b]->vars, basisstore->bases[b]->conss, basisstore->bases[b]->vstat,
            basisstore->bases[b]->cstat, basisstore->bases[b]->nvars, basisstore->bases[b]->nconss) );
   }

   return SCIP_OKAY;
}


void SCIPbasestoreRemoveBasis(
   SCIP_BASISSTORE*      basesstore
   )
{
   assert(basesstore != NULL);

   SCIPdebugMessage("remove last basis\n");

   BMSfreeMemoryArray(&basesstore->bases[0]->cstat);
   BMSfreeMemoryArray(&basesstore->bases[0]->vstat);
   BMSfreeMemoryArray(&basesstore->bases[0]->conss);
   BMSfreeMemoryArray(&basesstore->bases[0]->vars);
   BMSfreeMemory(&basesstore->bases[0]);

   basesstore->nbases--;
}

/**
 * Methods for SCIP_BASIS
 */

/* returns the data of the given basis */
void SCIPbaseGetData(
   SCIP_BASIS*           basis,
   SCIP_VAR***           vars,
   SCIP_CONS***          conss,
   int**                 vstat,
   int**                 cstat,
   int*                  nvars,
   int*                  nconss
   )
{
   assert(basis != NULL);

   (*vars) = basis->vars;
   (*conss) = basis->conss;
   (*vstat) = basis->vstat;
   (*cstat) = basis->cstat;
   (*nvars) = basis->nvars;
   (*nconss) = basis->nconss;
}
