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
 * @brief  methods for storing LP starting bases
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/scip.h"
#include "scip/basisstore.h"
#include "scip/struct_basisstore.h"


/**
 * Methods for SCIP_BASISSTORE
 */

/** initialize the basis storage */
static
SCIP_RETCODE initializeBasisstore(
   SCIP_BASISSTORE*      basisstore          /**< basis storage */
   )
{
   assert(basisstore != NULL);
   assert(basisstore->basessize == 0);

   SCIPdebugMessage("initialize basisstore %p\n", (void*) basisstore);

   SCIP_ALLOC( BMSallocMemoryArray(&basisstore->bases, 1) );
   basisstore->basessize = 1;

   return SCIP_OKAY;
}

/** ensure that the store is large enoug */
static
SCIP_RETCODE checkMem(
   SCIP_BASISSTORE*      basisstore          /**< basis storage */
   )
{
   assert(basisstore != NULL);

   if( basisstore->nbases == basisstore->basessize )
   {
      SCIP_ALLOC( BMSreallocMemoryArray(&basisstore->bases, 2*basisstore->basessize) );
      basisstore->basessize *= 2;
   }

   return SCIP_OKAY;
}

/** creates basis storage */
SCIP_RETCODE SCIPbasisstoreCreate(
   SCIP_BASISSTORE**     basisstore          /**< basis storage */
   )
{
   assert(basisstore != NULL);

   SCIP_ALLOC( BMSallocMemory(basisstore) );

   (*basisstore)->bases = NULL;
   (*basisstore)->basessize = 0;
   (*basisstore)->nbases = 0;

   SCIPdebugMessage("create basisstore %p\n", (void*)(*basisstore));

   return SCIP_OKAY;
}

/** frees basis storage */
SCIP_RETCODE SCIPbasisstoreFree(
   SCIP_BASISSTORE**     basisstore           /**< basis storage */
   )
{
   assert(basisstore != NULL);
   assert(*basisstore != NULL);

   if( (*basisstore)->nbases > 0 )
   {
      SCIPbasistoreRemoveBasis(*basisstore);
   }
   assert((*basisstore)->nbases == 0);

   if( (*basisstore)->bases != NULL )
   {
      BMSfreeMemoryArray(&(*basisstore)->bases);
   }

   SCIPdebugMessage("free basisstore %p.\n", (void*)(*basisstore));
   BMSfreeMemoryNull(&(*basisstore));
   *basisstore = NULL;

   return SCIP_OKAY;
}

/** add a basis to the storage */
SCIP_RETCODE SCIPbasisstoreAddBasis(
   SCIP_BASISSTORE*      basisstore,          /**< basis storage */
   SCIP_VAR**            vars,                /**< array of problem variables */
   SCIP_CONS**           conss,               /**< array of problem constraints */
   int*                  varstat,             /**< array of variable basis status */
   int*                  consstat,            /**< array of constraint basis status */
   int                   nvars,               /**< number of problem variables */
   int                   nconss               /**< number of problem constraints */
   )
{
   assert(basisstore != NULL);
   assert(vars != NULL);
   assert(conss != NULL);
   assert(varstat != NULL);
   assert(consstat != NULL);

   /* initialize the storage */
   if( basisstore->basessize == 0 )
   {
      SCIP_CALL( initializeBasisstore(basisstore) );
   }

   /* ensure that the storage is large enough */
   if( basisstore->nbases == basisstore->basessize )
   {
      SCIP_CALL( checkMem(basisstore) );
   }

   SCIPdebugMessage("add basis to basisstore %p: pos %d, nvars %d, nconss %d\n", (void*) basisstore, basisstore->nbases,
         nvars, nconss);

   /* store the basis */
   SCIP_ALLOC( BMSallocMemory(&basisstore->bases[basisstore->nbases]) );
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->vars, vars, nvars);
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->conss, conss, nconss);
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->varstat, varstat, nvars);
   BMSduplicateMemoryArray(&basisstore->bases[basisstore->nbases]->consstat, consstat, nconss);
   basisstore->bases[basisstore->nbases]->nvars = nvars;
   basisstore->bases[basisstore->nbases]->nconss = nconss;
   basisstore->nbases++;

   return SCIP_OKAY;
}

/** return the number of stored basis */
int SCIPbasisstoreGetNBasis(
   SCIP_BASISSTORE*      basisstore           /**< basis storage */
   )
{
   assert(basisstore != NULL);

   return basisstore->nbases;
}

/** returns the stored basis */
SCIP_BASIS* SCIPbasistoreGetBasis(
   SCIP_BASISSTORE*      basistore           /**< basis storage */
   )
{
   assert(basistore != NULL);
   assert(basistore->bases != NULL);
   assert(basistore->nbases == 1); /* tmp: need to be changed if we store more than one basis */

   return basistore->bases[0];
}

/** copy basis storage */
SCIP_RETCODE SCIPbasisstoreCopy(
   SCIP_BASISSTORE*      basisstore,         /**< source basis storage */
   SCIP_BASISSTORE*      targetbasisstore    /**< target basis storage */
   )
{
   int b;

   assert(basisstore != NULL);
   assert(targetbasisstore != NULL);

   SCIPdebugMessage("copy basisstore: %p -> %p\n", (void*) basisstore, (void*) targetbasisstore);

   for(b = 0; b < basisstore->nbases; b++)
   {
      assert(basisstore->bases[b]->varstat != NULL);
      assert(basisstore->bases[b]->consstat != NULL);

      SCIP_CALL( SCIPbasisstoreAddBasis(targetbasisstore,
            basisstore->bases[b]->vars, basisstore->bases[b]->conss, basisstore->bases[b]->varstat,
            basisstore->bases[b]->consstat, basisstore->bases[b]->nvars, basisstore->bases[b]->nconss) );
   }

   return SCIP_OKAY;
}

/** remove last added basis */
void SCIPbasistoreRemoveBasis(
   SCIP_BASISSTORE*      basisstore         /**< source basis storage */
   )
{
   assert(basisstore != NULL);
   assert(basisstore->nbases > 0);

   basisstore->nbases--;
   SCIPdebugMessage("delete last basis at position %d in basisstore %p\n", basisstore->nbases, (void*) basisstore);

   BMSfreeMemoryArray(&basisstore->bases[basisstore->nbases]->consstat);
   BMSfreeMemoryArray(&basisstore->bases[basisstore->nbases]->varstat);
   BMSfreeMemoryArray(&basisstore->bases[basisstore->nbases]->conss);
   BMSfreeMemoryArray(&basisstore->bases[basisstore->nbases]->vars);
   BMSfreeMemory(&basisstore->bases[basisstore->nbases]);
}

/**
 * Methods for SCIP_BASIS
 */

/* returns the data of the given basis */
void SCIPbasisGetData(
   SCIP_BASIS*           basis,         /**< source basis storage */
   SCIP_VAR***           vars,          /**< array to store the problem variables */
   SCIP_CONS***          conss,         /**< array to store the problem constraints */
   int**                 varstat,       /**< array to store basis status of variables */
   int**                 consstat,      /**< array to store basis status of constraints */
   int*                  nvars,         /**< pointer to store number of variabls */
   int*                  nconss         /**< pointer to store number of constraints */
   )
{
   assert(basis != NULL);
   assert(vars != NULL);
   assert(conss != NULL);
   assert(varstat != NULL);
   assert(consstat != NULL);
   assert(nvars != NULL);
   assert(nconss != NULL);

   (*vars) = basis->vars;
   (*conss) = basis->conss;
   (*varstat) = basis->varstat;
   (*consstat) = basis->consstat;
   (*nvars) = basis->nvars;
   (*nconss) = basis->nconss;
}
