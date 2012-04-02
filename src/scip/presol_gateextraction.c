/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_gateextraction.c
 * @brief  gateextraction presolver
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/presol_gateextraction.h"
#include "scip/cons_setppc.h"
#include "scip/cons_logicor.h"
#include "scip/cons_and.h"


#define PRESOL_NAME            "gateextraction"
#define PRESOL_DESC            "presolver template"
#define PRESOL_PRIORITY         1000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY               TRUE /**< should presolver be delayed, if other presolvers found reductions? */

#define HASHSIZE_LOGICORCONS     131101 /**< minimal size of hash table in logicor constraint tables */
#define HASHSIZE_SETPPCCONS      131101 /**< minimal size of hash table in setppc constraint tables */

#define DEFAULT_ONLYSETPART       FALSE  /**< should only set-partitioning constraints be extrated and no and-constraints */


/*
 * Data structures
 */


/** data object to compare constraint easier */
struct HashData
{
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   int nvars;
};
typedef struct HashData HASHDATA;


/** presolver data */
struct SCIP_PresolData
{
   HASHDATA** setppchashdatas;
   HASHDATA* setppchashdatastore;
   SCIP_HASHTABLE* hashdatatable;
   SCIP_HASHTABLE* setppchashtable;
   SCIP_HASHTABLE* logicorhashtable;
   SCIP_CONS** usefullogicor;
   int nusefullogicor;
   int susefullogicor;
   int nsetppchashdatas;
   int ssetppchashdatas;
   int ngates;
   int firstchangedlogicor;
   SCIP_Bool usefulsetppcexist;
   SCIP_Bool usefullogicorexist;
   SCIP_Bool newsetppchashdatas;
   SCIP_Bool initialized;
   SCIP_Bool onlysetpart;
};


/*
 * Local methods
 */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyCons)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same pointer */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqCons)
{
#ifndef NDEBUG
   SCIP* scip;

   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   return (key1 == key2);
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValCons)
{  /*lint --e{715}*/
   /* the key is used as the keyvalue too */
   return (unsigned int)(size_t) key;
}


/* put your local methods here, and declare them static */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashdataGetKeyCons)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same pointer */
static
SCIP_DECL_HASHKEYEQ(hashdataKeyEqCons)
{
#ifndef NDEBUG
   SCIP* scip;
#endif
   HASHDATA* hashdata1;
   HASHDATA* hashdata2;
   int v;

   hashdata1 = (HASHDATA*)key1;
   hashdata2 = (HASHDATA*)key2;
#ifndef NDEBUG
   scip = (SCIP*)userptr;
   assert(scip != NULL);
#endif

   /* check data structure */
   assert(hashdata1->nvars == 2);
   assert(hashdata2->nvars == 2);
   /* at least one data object needs to be have a real set packing constraint */
   assert(hashdata1->cons != NULL || hashdata2->cons != NULL);

   for( v = 1; v >= 0; --v )
   {
      /* tests if variables are equal */
      if( hashdata1->vars[v] != hashdata2->vars[v] )
	 return FALSE;

      assert(SCIPvarCompare(hashdata1->vars[v], hashdata2->vars[v]) == 0);
   }

   /* a hashdata object is only equal if it has the same constraint pointer, or one has no constraint pointer, latter
    * means that this object is a form a logicor constraint derived hashdata object
    */
   if( hashdata1->cons == NULL || hashdata2->cons == NULL || hashdata1->cons == hashdata2->cons )
      return TRUE;
   else
      return FALSE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashdataKeyValCons)
{  /*lint --e{715}*/
   HASHDATA* hashdata;
   unsigned int hashval;
#if 0
   int v;
#endif

   hashdata = (HASHDATA*)key;
   assert(hashdata != NULL);
   assert(hashdata->vars != NULL);
   assert(hashdata->nvars == 2);

#if 0
   hashval = 0;
   for( v = 1; v >= 0; --v )
      hashval |= (1<<(SCIPvarGetIndex(hashdata->vars[v])%32));
#else
   hashval = (SCIPvarGetIndex(hashdata->vars[1])<<16) + SCIPvarGetIndex(hashdata->vars[0]); /*lint !e701*/
#endif

   return hashval;
}



/** initialize gateextraction presolver data */
static
SCIP_RETCODE presoldataInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< data object of presolver */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);

   presoldata->usefullogicor = NULL;
   presoldata->nusefullogicor = 0;
   presoldata->susefullogicor = 0;
   presoldata->firstchangedlogicor = -1;
   presoldata->nsetppchashdatas = 0;
   presoldata->ssetppchashdatas = 0;
   presoldata->ngates = 0;
   presoldata->usefulsetppcexist = FALSE;
   presoldata->usefullogicorexist = FALSE;
   presoldata->newsetppchashdatas = FALSE;
   presoldata->initialized = FALSE;

   SCIP_CALL( SCIPhashtableCreate(&(presoldata->hashdatatable), SCIPblkmem(scip), HASHSIZE_SETPPCCONS, hashdataGetKeyCons, hashdataKeyEqCons, hashdataKeyValCons, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&(presoldata->setppchashtable), SCIPblkmem(scip), HASHSIZE_SETPPCCONS, hashGetKeyCons, hashKeyEqCons, hashKeyValCons, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&(presoldata->logicorhashtable), SCIPblkmem(scip), HASHSIZE_LOGICORCONS, hashGetKeyCons, hashKeyEqCons, hashKeyValCons, (void*) scip) );

   return SCIP_OKAY;
}


/** create useful set-packing information by adding new set-packing constraints with two variables */
static
SCIP_RETCODE createPresoldata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< data object of presolver */
   SCIP_CONS**           setppcs,            /**< active setppc constraints */
   int                   nsetppcs,           /**< number of active setppc constraints */
   SCIP_CONS**           logicors,           /**< active logicor constraints */
   int                   nlogicors           /**< number of active logicor constraints */
   )
{
   SCIP_CONS** usefulconss;
   int nusefulconss = 0;
   int size;
   int c;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(setppcs != NULL);
   assert(nsetppcs > 0);
   assert(logicors != NULL);
   assert(nlogicors > 0);
   assert(presoldata->setppchashtable != NULL);
   assert(presoldata->logicorhashtable != NULL);

   presoldata->initialized = TRUE;

   size = MAX(nsetppcs, nlogicors);

   /* temporary memory for collecting set-packing constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &usefulconss, size) );

   if( !presoldata->usefulsetppcexist )
   {
      /* find set-packing constraints with exactly two varibales */
      for( c = 0; c < nsetppcs; ++c )
      {
	 assert(SCIPconsIsActive(setppcs[c]));

	 if( SCIPgetTypeSetppc(scip, setppcs[c]) == SCIP_SETPPCTYPE_PACKING && SCIPgetNVarsSetppc(scip, setppcs[c]) == 2 && !SCIPconsIsModifiable(setppcs[c]) )
	 {
	    /* insert new element in hashtable */
	    SCIP_CALL( SCIPhashtableInsert(presoldata->setppchashtable, (void*) setppcs[c]) );

	    usefulconss[nusefulconss] = setppcs[c];
	    ++nusefulconss;
	 }
      }

      /* add usefulconss constraints to hashdata elements */
      if( nusefulconss > 0 )
      {
	 SCIP_Bool negated[2];
	 int h;

	 presoldata->usefulsetppcexist = TRUE;
	 presoldata->ssetppchashdatas = nusefulconss;

	 SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->setppchashdatas), nusefulconss) );
	 SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(presoldata->setppchashdatastore), nusefulconss) );

	 h = 0;
	 for( c = 0; c < nusefulconss; ++c )
	 {
	    assert(SCIPconsIsActive(usefulconss[c]));
	    assert(SCIPgetNVarsSetppc(scip, usefulconss[c]) == 2);
	    presoldata->setppchashdatas[h] = &(presoldata->setppchashdatastore[h]);

	    SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(presoldata->setppchashdatas[h]->vars), SCIPgetVarsSetppc(scip, usefulconss[c]), 2) );

	    SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[h]->vars[0], &(presoldata->setppchashdatas[h]->vars[0]), &(negated[0])) );
	    SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[h]->vars[1], &(presoldata->setppchashdatas[h]->vars[1]), &(negated[1])) );

	    if( SCIPvarGetStatus(presoldata->setppchashdatas[h]->vars[0]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[h]->vars[0]) == SCIP_VARSTATUS_MULTAGGR
	       || SCIPvarGetStatus(presoldata->setppchashdatas[h]->vars[1]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[h]->vars[1]) == SCIP_VARSTATUS_MULTAGGR )
	    {
	       SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[h]->vars), 2);
	       continue;
	    }

	    presoldata->setppchashdatas[h]->nvars = 2;

	    /* capture variables */
	    SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[h]->vars[0]) );
	    SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[h]->vars[1]) );

	    /* order the variables after their index */
	    if( SCIPvarGetIndex(presoldata->setppchashdatas[h]->vars[0]) > SCIPvarGetIndex(presoldata->setppchashdatas[h]->vars[1]) )
	    {
	       SCIP_VAR* tmp = presoldata->setppchashdatas[h]->vars[0];
	       presoldata->setppchashdatas[h]->vars[0] = presoldata->setppchashdatas[h]->vars[1];
	       presoldata->setppchashdatas[h]->vars[1] = tmp;
	    }

	    presoldata->setppchashdatas[h]->cons = usefulconss[c];

	    SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[h]) );
	    SCIP_CALL( SCIPcaptureCons(scip, usefulconss[c]) );

	    ++h;
	 }
	 presoldata->nsetppchashdatas = h;
      }
   }

   nusefulconss = 0;

   if( !presoldata->usefullogicorexist )
   {
      /* find logicor constraints with exactly three varibales */
      for( c = 0; c < nlogicors; ++c )
      {
	 assert(SCIPconsIsActive(logicors[c]));

	 if( SCIPgetNVarsLogicor(scip, logicors[c]) == 3  && !SCIPconsIsModifiable(logicors[c]) )
	 {
	    /* insert new element in hashtable */
	    SCIP_CALL( SCIPhashtableInsert(presoldata->logicorhashtable, (void*) logicors[c]) );
	    SCIP_CALL( SCIPcaptureCons(scip, logicors[c]) );

	    usefulconss[nusefulconss] = logicors[c];
	    ++nusefulconss;
	 }
      }

      /* no usefulconss constraints */
      if( nusefulconss > 0 )
      {
	 presoldata->firstchangedlogicor = 0;
	 presoldata->usefullogicorexist = TRUE;
	 presoldata->susefullogicor = nusefulconss;
	 presoldata->nusefullogicor = nusefulconss;
	 SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &presoldata->usefullogicor, usefulconss, presoldata->susefullogicor) );
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &usefulconss);

   return SCIP_OKAY;
}


/** remove old setppchashdatas objects, so that the allocated memory will stay low */
static
SCIP_RETCODE cleanupHashDatas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata          /**< data object of presolver */
   )
{
   assert(scip != NULL);
   assert(presoldata != NULL);

   if( presoldata->usefulsetppcexist )
   {
      int c;

      assert(presoldata->setppchashdatas != NULL || presoldata->nsetppchashdatas == 0);

      for( c = presoldata->nsetppchashdatas - 1; c >= 0; --c )
      {
	 SCIP_Bool removeentry = FALSE;

	 assert(presoldata->setppchashdatas[c] != NULL);
	 assert(presoldata->setppchashdatas[c]->cons != NULL);

	 if( SCIPconsIsDeleted(presoldata->setppchashdatas[c]->cons) || SCIPconsIsModifiable(presoldata->setppchashdatas[c]->cons)
	    || SCIPgetTypeSetppc(scip, presoldata->setppchashdatas[c]->cons) != SCIP_SETPPCTYPE_PACKING || SCIPgetNVarsSetppc(scip, presoldata->setppchashdatas[c]->cons) != 2 )
	 {
	    removeentry = TRUE;
	 }
	 else
	 {
	    SCIP_VAR* vars[2];
	    SCIP_Bool negated[2];

	    SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[c]->vars[0], &(vars[0]), &(negated[0])) );
	    SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[c]->vars[1], &(vars[1]), &(negated[1])) );

	    if( SCIPvarGetStatus(vars[0]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(vars[0]) == SCIP_VARSTATUS_MULTAGGR
	       || SCIPvarGetStatus(vars[1]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(vars[1]) == SCIP_VARSTATUS_MULTAGGR
	       || presoldata->setppchashdatas[c]->vars[0] != vars[0] || presoldata->setppchashdatas[c]->vars[1] != vars[1] )
	    {
	       removeentry = TRUE;
	    }
	 }

	 if( removeentry )
	 {
	    /* remove constraint from setppc-hashtable */
	    assert(SCIPhashtableExists(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c]->cons));
	    SCIP_CALL( SCIPhashtableRemove(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c]->cons) );

	    /* remove hashdata entry from hashtable */
	    SCIP_CALL( SCIPhashtableRemove(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[c]) );

	    /* release old constraints */
	    SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->setppchashdatas[c]->cons)) );

	    /* release variables */
	    SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c]->vars[0])) );
	    SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c]->vars[1])) );

	    /* free memory for variables */
	    SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[c]->vars), 2);

	    if( c < presoldata->nsetppchashdatas - 1 )
	    {
	       /* remove old hashdata entry from hashtable */
	       SCIP_CALL( SCIPhashtableRemove(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1]) );
	    }

	    /* move last content to free position */
	    presoldata->setppchashdatas[c]->cons = presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1]->cons;
	    presoldata->setppchashdatas[c]->vars = presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1]->vars;
	    presoldata->setppchashdatas[c]->nvars = presoldata->setppchashdatas[presoldata->nsetppchashdatas - 1]->nvars;

	    if( c < presoldata->nsetppchashdatas - 1 )
	    {
	       /* add new hashdata entry from hashtable */
	       SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[c]) );
	    }
	    --(presoldata->nsetppchashdatas);
	 }
      }

#ifndef NDEBUG
      for( c = presoldata->nsetppchashdatas - 1; c >= 0; --c )
      {
	 assert(presoldata->setppchashdatas[c] != NULL);
	 assert(presoldata->setppchashdatas[c]->nvars == 2);
	 assert(presoldata->setppchashdatas[c]->vars != NULL);
	 assert(presoldata->setppchashdatas[c]->vars[0] != NULL);
	 assert(presoldata->setppchashdatas[c]->vars[1] != NULL);
	 assert(presoldata->setppchashdatas[c]->cons != NULL);
	 assert(SCIPconsIsActive(presoldata->setppchashdatas[c]->cons));
	 assert(SCIPhashtableExists(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[c]));
	 assert(SCIPhashtableExists(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c]->cons));
      }
#endif
   }

   return SCIP_OKAY;
}

/** refresh useful set-packing information, delete redundant constraints and add new constraints */
static
SCIP_RETCODE correctPresoldata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< data object of presolver */
   SCIP_CONS**           setppcs,            /**< active setppc constraints */
   int                   nsetppcs,           /**< number of active setppc constraints */
   SCIP_CONS**           logicors,           /**< active setppc constraints */
   int                   nlogicors           /**< number of active setppc constraints */
   )
{
   int oldnsetppchashdatas;
   int c;

   assert(scip != NULL);
   assert(presoldata != NULL);
   assert(setppcs != NULL);
   assert(nsetppcs > 0);
   assert(logicors != NULL);
   assert(nlogicors > 0);
   assert(presoldata->initialized);
   assert(presoldata->setppchashtable != NULL);
   assert(presoldata->logicorhashtable != NULL);

   /* check if there already exist some set-packing and some logicor constraints with the right amount of variables */
   if( !presoldata->usefulsetppcexist || !presoldata->usefullogicorexist )
   {
      SCIP_CALL( createPresoldata(scip, presoldata, setppcs, nsetppcs, logicors, nlogicors) );

      /* no correct logicor or set-packing constraints available, so abort */
      if( !presoldata->usefulsetppcexist || !presoldata->usefullogicorexist )
	 return SCIP_OKAY;
   }

   /* correct old data */
   SCIP_CALL( cleanupHashDatas(scip, presoldata) );

   oldnsetppchashdatas = presoldata->nsetppchashdatas;

   /* first update setppc part */
   /* add new setppc constraints */
   for( c = nsetppcs - 1; c >= 0; --c )
   {
      assert(SCIPconsIsActive(setppcs[c]));

      if( SCIPgetTypeSetppc(scip, setppcs[c]) == SCIP_SETPPCTYPE_PACKING && SCIPgetNVarsSetppc(scip, setppcs[c]) == 2 && !SCIPconsIsModifiable(setppcs[c]) )
      {
	 /* check if constraint is new, and correct array size if necessary */
	 if( !SCIPhashtableExists(presoldata->setppchashtable, (void*) setppcs[c]) )
	 {
	    SCIP_Bool negated[2];

	    /* resize array if necessary */
	    if( presoldata->nsetppchashdatas == presoldata->ssetppchashdatas )
	    {
	       int newsize;
	       int d;

	       newsize = SCIPcalcMemGrowSize(scip, presoldata->nsetppchashdatas + 1);

	       /* array already at maximal size */
	       if( newsize <= presoldata->ssetppchashdatas )
		  return SCIP_NOMEMORY;

	       /* correct hashtable, remove old elements */
	       SCIPhashtableClear(presoldata->hashdatatable);

	       SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(presoldata->setppchashdatas), presoldata->ssetppchashdatas, newsize) );
	       SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(presoldata->setppchashdatastore), presoldata->ssetppchashdatas, newsize) );
	       presoldata->ssetppchashdatas = newsize;

	       /* correct pointers in array, to point to the new storage, due to allocation, and add all elements to the
		* hashtable again
		*/
	       for( d = presoldata->nsetppchashdatas - 1; d >= 0; --d )
	       {
		  presoldata->setppchashdatas[d] = &(presoldata->setppchashdatastore[d]);
		  SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[d]) );
	       }
	    }

	    /* insert new element in hashtable */
	    SCIP_CALL( SCIPhashtableInsert(presoldata->setppchashtable, (void*) setppcs[c]) );

	    assert(SCIPgetNVarsSetppc(scip, setppcs[c]) == 2);
	    presoldata->setppchashdatas[presoldata->nsetppchashdatas] = &(presoldata->setppchashdatastore[presoldata->nsetppchashdatas]);

	    SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars), SCIPgetVarsSetppc(scip, setppcs[c]), 2) );
	    SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0], &(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0]), &(negated[0])) );
	    SCIP_CALL( SCIPgetBinvarRepresentative(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1], &(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1]), &(negated[1])) );

	    if( SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0]) == SCIP_VARSTATUS_MULTAGGR
	       || SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1]) == SCIP_VARSTATUS_FIXED || SCIPvarGetStatus(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1]) == SCIP_VARSTATUS_MULTAGGR )
	    {
	       SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars), 2);
	       continue;
	    }

	    presoldata->setppchashdatas[presoldata->nsetppchashdatas]->nvars = 2;

	    /* capture variables */
	    SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0]) );
	    SCIP_CALL( SCIPcaptureVar(scip, presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1]) );

	    /* order the variables after their index */
	    if( SCIPvarGetIndex(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0]) > SCIPvarGetIndex(presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1]) )
	    {
	       SCIP_VAR* tmp = presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0];
	       presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[0] = presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1];
	       presoldata->setppchashdatas[presoldata->nsetppchashdatas]->vars[1] = tmp;
	    }

	    presoldata->setppchashdatas[presoldata->nsetppchashdatas]->cons = setppcs[c];

	    SCIP_CALL( SCIPhashtableInsert(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[presoldata->nsetppchashdatas]) );
	    SCIP_CALL( SCIPcaptureCons(scip, setppcs[c]) );

	    ++(presoldata->nsetppchashdatas);
	 }
      }
   }

   /* if we found new set-packing constraints, we want to check against all logicors */
   if( oldnsetppchashdatas < presoldata->nsetppchashdatas )
      presoldata->newsetppchashdatas = TRUE;

   /* now logicor part */
   /* removed last deleted logicor constraints from local presolver data */
   while( presoldata->nusefullogicor > 0 && !SCIPconsIsActive(presoldata->usefullogicor[presoldata->nusefullogicor - 1]) )
   {
      SCIP_CALL( SCIPhashtableRemove(presoldata->logicorhashtable, (void*) presoldata->usefullogicor[presoldata->nusefullogicor - 1]) );
      SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->usefullogicor[presoldata->nusefullogicor - 1])) );

      --(presoldata->nusefullogicor);
   }

   /* remove old inactive logicor constraints */
   for( c = presoldata->nusefullogicor - 1; c >= 0; --c )
   {
      if( !SCIPconsIsActive(presoldata->usefullogicor[c]) || SCIPgetNVarsLogicor(scip, presoldata->usefullogicor[c]) != 3 )
      {
	 SCIP_CALL( SCIPhashtableRemove(presoldata->logicorhashtable, (void*) presoldata->usefullogicor[c]) );
	 SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->usefullogicor[c])) );

	 presoldata->usefullogicor[c] = presoldata->usefullogicor[presoldata->nusefullogicor - 1];
	 --(presoldata->nusefullogicor);
      }
   }

   presoldata->firstchangedlogicor = presoldata->nusefullogicor;
   assert(presoldata->firstchangedlogicor >= 0);

   /* add new logicor constraints */
   for( c = nlogicors - 1; c >= 0; --c )
   {
      assert(SCIPconsIsActive(logicors[c]));

      if( SCIPgetNVarsLogicor(scip, logicors[c]) == 3 && !SCIPconsIsModifiable(logicors[c]) )
      {
	 /* check if constraint is new, and correct array size if necessary */
	 if( !SCIPhashtableExists(presoldata->logicorhashtable, (void*) logicors[c]) )
	 {
	    /* resize array if necessary */
	    if( presoldata->nusefullogicor == presoldata->susefullogicor )
	    {
	       int newsize;

	       newsize = SCIPcalcMemGrowSize(scip, presoldata->nusefullogicor + 1);

	       /* array already at maximal size */
	       if( newsize <= presoldata->susefullogicor )
		  return SCIP_NOMEMORY;

	       SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(presoldata->usefullogicor), presoldata->susefullogicor, newsize) );
	       presoldata->susefullogicor = newsize;
	    }

	    /* insert new element in hashtable */
	    SCIP_CALL( SCIPhashtableInsert(presoldata->logicorhashtable, (void*) logicors[c]) );
	    SCIP_CALL( SCIPcaptureCons(scip, logicors[c]) );

	    presoldata->usefullogicor[presoldata->nusefullogicor] = logicors[c];
	    ++(presoldata->nusefullogicor);
	 }
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyGateextraction)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolGateextraction(scip) );

   return SCIP_OKAY;
}


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeGateextraction)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPhashtableFree(&(presoldata->logicorhashtable));
   SCIPhashtableFree(&(presoldata->setppchashtable));
   SCIPhashtableFree(&(presoldata->hashdatatable));

   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** initialization method of presolver (called after problem was transformed) */
#define presolInitGateextraction NULL


/** deinitialization method of presolver (called before transformed problem is freed) */
static
SCIP_DECL_PRESOLEXIT(presolExitGateextraction)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   int c;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* release old constraints */
   for( c = presoldata->nusefullogicor - 1; c >= 0; --c )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->usefullogicor[c])) );
   }

   if( presoldata->usefullogicorexist )
   {
      SCIPfreeBlockMemoryArray(scip, &presoldata->usefullogicor, presoldata->susefullogicor);
   }

   if( presoldata->usefulsetppcexist )
   {
      assert(presoldata->setppchashdatas != NULL || presoldata->nsetppchashdatas == 0);
      for( c = presoldata->nsetppchashdatas - 1; c >= 0; --c )
      {
	 assert(presoldata->setppchashdatas[c] != NULL);
	 assert(presoldata->setppchashdatas[c]->cons != NULL);
	 assert(presoldata->setppchashdatas[c]->vars != NULL);

	 /* remove constraint from setppc-hashtable */
	 assert(SCIPhashtableExists(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c]->cons));
	 SCIP_CALL( SCIPhashtableRemove(presoldata->setppchashtable, (void*) presoldata->setppchashdatas[c]->cons) );

	 /* release old constraints */
	 SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->setppchashdatas[c]->cons)) );

	 /* remove hashdata entry from hashtable */
	 SCIP_CALL( SCIPhashtableRemove(presoldata->hashdatatable, (void*) presoldata->setppchashdatas[c]) );

	 /* release variables */
	 SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c]->vars[0])) );
	 SCIP_CALL( SCIPreleaseVar(scip, &(presoldata->setppchashdatas[c]->vars[1])) );

	 /* free memory for variables */
	 SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas[c]->vars), 2);
      }

      SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatastore), presoldata->ssetppchashdatas);
      SCIPfreeBlockMemoryArray(scip, &(presoldata->setppchashdatas), presoldata->ssetppchashdatas);
   }

   presoldata->nusefullogicor = 0;
   presoldata->susefullogicor = 0;
   presoldata->nsetppchashdatas = 0;
   presoldata->ssetppchashdatas = 0;
   presoldata->firstchangedlogicor = -1;
   presoldata->ngates = 0;
   presoldata->usefullogicorexist = FALSE;
   presoldata->usefulsetppcexist = FALSE;
   presoldata->newsetppchashdatas = FALSE;
   presoldata->initialized = FALSE;

   return SCIP_OKAY;
}


/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreGateextraction)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreGateextraction)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}



#define HASHTABLESIZE_FACTOR 5

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecGateextraction)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_CONSHDLR* conshdlrsetppc;
   SCIP_CONSHDLR* conshdlrlogicor;
   SCIP_CONS** setppcconss;
   SCIP_CONS** logicorconss;
   int nsetppcconss;
   int nlogicorconss;
   int c;
   SCIP_Bool paramvalue;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get necessary constraint handlers */
   conshdlrsetppc = SCIPfindConshdlr(scip, "setppc");
   conshdlrlogicor = SCIPfindConshdlr(scip, "logicor");

   if( conshdlrsetppc == NULL || conshdlrlogicor == NULL )
      return SCIP_OKAY;

   /* get number of active constraints */
   nsetppcconss = SCIPconshdlrGetNActiveConss(conshdlrsetppc);
   assert(nsetppcconss >= 0);
   nlogicorconss = SCIPconshdlrGetNActiveConss(conshdlrlogicor);
   assert(nlogicorconss >= 0);

   if( nsetppcconss == 0 || nlogicorconss == 0 )
      return SCIP_OKAY;

   paramvalue = FALSE;
   if( SCIPgetBoolParam(scip, "constraints/and/linearize", &paramvalue) == SCIP_OKAY )
      if( paramvalue )
      {
	 SCIPwarningMessage(scip, "Gate-presolving is the counterpart of linearizing all and constraints, so enabling both presolving steps at ones does not make sense.\n");
      }

   *result = SCIP_DIDNOTFIND;

   /* get active constraints */
   nsetppcconss = SCIPconshdlrGetNActiveConss(conshdlrsetppc);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &setppcconss, SCIPconshdlrGetConss(conshdlrsetppc), nsetppcconss) );

   assert(setppcconss != NULL);
   logicorconss = SCIPconshdlrGetConss(conshdlrlogicor);
   assert(logicorconss != NULL);

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   presoldata->newsetppchashdatas = FALSE;

   /* @todo - we only extract and-gates with two operands, it is possible to expand the functionality to also extract
    *         bigger and-gastes
    *       - improve speed of finding gates by adding a hash value to each constraint and do a hash comparison
    */

   if( !presoldata->initialized )
   {
      assert(presoldata->usefullogicor == NULL);

      /* create useful set-packing information by adding new set-packing constraints with two variables */
      SCIP_CALL( createPresoldata(scip, presoldata, setppcconss, nsetppcconss, logicorconss, nlogicorconss) );
   }
   else
   {
      /* refresh useful set-packing information, delete redundant constraints and add new constraints */
      SCIP_CALL( correctPresoldata(scip, presoldata, setppcconss, nsetppcconss, logicorconss, nlogicorconss) );
   }
   assert(presoldata->initialized);

   /* we do not have any useful set-packing or logicor constraint, or since last run did not get any new constraints, so abort */
   if( presoldata->nsetppchashdatas == 0 || presoldata->nusefullogicor == 0 || presoldata->firstchangedlogicor == presoldata->nusefullogicor )
      goto TERMINATE;

   assert(presoldata->usefullogicor != NULL);
   assert(presoldata->nusefullogicor > 0);
   assert(presoldata->firstchangedlogicor >= 0);
   assert(presoldata->nsetppchashdatas > 0);

   if( presoldata->nsetppchashdatas > 0 )
   {
      SCIP_HASHMAP* varmap;
      SCIP_CONS* logicor;
      int size;
      int endloop;

      /* due to negations we can have at most two times the number of binary variables */
      size = 2 * (SCIPgetNBinVars(scip) + SCIPgetNImplVars(scip));

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(scip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * size)) );

      /* if we found new setppcs we want to check all logicors again */
      if( presoldata->newsetppchashdatas )
	 endloop = 0;
      else
	 endloop = presoldata->firstchangedlogicor;

      /* check all (new) logicors against all set-packing constraints constraints */
      for( c = presoldata->nusefullogicor - 1; c >= endloop; --c )
      {
	 logicor = presoldata->usefullogicor[c];
	 assert(logicor != NULL);
	 assert(SCIPgetNVarsLogicor(scip, logicor) == 3);
	 assert(!SCIPconsIsModifiable(logicor));

	 if( SCIPconsIsActive(logicor) )
	 {
	    HASHDATA* hashdata;
	    HASHDATA* hashmaphashdata;
	    SCIP_CONS* gateconss[3];
	    SCIP_VAR** logicorvars;
	    SCIP_VAR* activevars[3];
	    SCIP_VAR* tmpvars[2];
	    SCIP_Bool negated[3];
	    int nfound = 0;
	    int d;

	    /* logicor constraint has the form: x + y + z >= 1
	     *
	     * find set-packing constraints:  (~x + ~y >= 1 and ~x + ~z >= 1)  <=>  (x + y <= 1 and x + z <= 1)
	     *                           or:  (~y + ~x >= 1 and ~y + ~z >= 1)  <=>  (y + x <= 1 and y + z <= 1)
	     *                           or:  (~z + ~x >= 1 and ~z + ~y >= 1)  <=>  (z + x <= 1 and z + y <= 1)
	     *
	     * these three constraints are aquivalent to: x = ~y * ~z (x = AND(~y,~z))
	     */

	    logicorvars = SCIPgetVarsLogicor(scip, logicor);
	    assert(logicorvars != NULL);

	    /* get active logicor variables */
	    for( d = 0; d < 3; ++d )
	    {
	       activevars[d] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, logicorvars[d]);
	       if( activevars[d] == NULL )
	       {
		  SCIP_CALL( SCIPgetBinvarRepresentative(scip, logicorvars[d], &(activevars[d]), &(negated[d])) );
		  SCIP_CALL( SCIPhashmapInsert(varmap, logicorvars[d], activevars[d]) );
	       }
	    }
	    SCIPsortPtr((void**)activevars, SCIPvarComp, 3);

	    assert(SCIPvarGetIndex(activevars[0]) <= SCIPvarGetIndex(activevars[1]));
	    assert(SCIPvarGetIndex(activevars[1]) <= SCIPvarGetIndex(activevars[2]));

	    /* check that we have really three different variables, if not remove the constraint from the hashmap and the data storage */
	    if( SCIPvarGetIndex(activevars[0]) == SCIPvarGetIndex(activevars[1]) || SCIPvarGetIndex(activevars[1]) == SCIPvarGetIndex(activevars[2]) )
	    {
	       SCIP_CALL( SCIPhashtableRemove(presoldata->logicorhashtable, (void*) presoldata->usefullogicor[c]) );
	       SCIP_CALL( SCIPreleaseCons(scip, &(presoldata->usefullogicor[c])) );

	       presoldata->usefullogicor[c] = presoldata->usefullogicor[presoldata->nusefullogicor - 1];
	       --(presoldata->nusefullogicor);

	       continue;
	    }

	    BMSclearMemoryArray(gateconss, 3);

	    SCIP_CALL( SCIPallocBuffer(scip, &hashdata) );
	    hashdata->nvars = 2;
	    hashdata->cons = NULL;

	    for( d = 0; d < 3; ++d )
	    {

	       if( d == 0 )
	       {
		  tmpvars[0] = activevars[0];
		  tmpvars[1] = activevars[1];
	       }
	       else if ( d == 1 )
	       {
		  tmpvars[0] = activevars[0];
		  tmpvars[1] = activevars[2];
	       }
	       else
	       {
		  tmpvars[0] = activevars[1];
		  tmpvars[1] = activevars[2];
	       }

	       hashdata->vars = tmpvars;

	       hashmaphashdata = (HASHDATA*) SCIPhashtableRetrieve(presoldata->hashdatatable, (void*) hashdata);

	       if( hashmaphashdata != NULL && SCIPconsIsActive(hashmaphashdata->cons) )
	       {
		  gateconss[d] = hashmaphashdata->cons;
		  ++nfound;
	       }
	    }

	    SCIPfreeBuffer(scip, &hashdata);

	    /* did we find three set-packing constraints for the upgrade to a set-partitioning constraint */
	    if( nfound == 3 )
	    {
	       SCIP_CONS* newcons;
	       char name[SCIP_MAXSTRLEN];
	       SCIP_Bool initial;
	       SCIP_Bool separate;
	       SCIP_Bool enforce;
	       SCIP_Bool check;
	       SCIP_Bool propagate;
	       SCIP_Bool local;
	       SCIP_Bool modifiable;
	       SCIP_Bool dynamic;
	       SCIP_Bool removable;
	       SCIP_Bool stickingatnode;
	       int i;

	       initial = SCIPconsIsInitial(logicor);
	       separate = SCIPconsIsSeparated(logicor);
	       enforce = SCIPconsIsEnforced(logicor);
	       check = SCIPconsIsChecked(logicor);
	       propagate = SCIPconsIsPropagated(logicor);
	       local = SCIPconsIsLocal(logicor);
	       modifiable = SCIPconsIsModifiable(logicor);
	       dynamic = SCIPconsIsDynamic(logicor);
	       removable = SCIPconsIsRemovable(logicor);
	       stickingatnode = SCIPconsIsStickingAtNode(logicor);

	       SCIPdebugMessage("Following four constraints form an set-partitioning constraint.\n");

	       for( i = 2; i >= 0; --i )
	       {
		  assert(gateconss[i] != NULL);

		  initial |= SCIPconsIsInitial(gateconss[i]);
		  separate |= SCIPconsIsSeparated(gateconss[i]);
		  enforce |= SCIPconsIsEnforced(gateconss[i]);
		  check |= SCIPconsIsChecked(gateconss[i]);
		  propagate |= SCIPconsIsPropagated(gateconss[i]);
		  local &= SCIPconsIsLocal(gateconss[i]);
		  modifiable &= SCIPconsIsModifiable(gateconss[i]);
		  dynamic &= SCIPconsIsDynamic(gateconss[i]);
		  removable &= SCIPconsIsRemovable(gateconss[i]);
		  stickingatnode &= SCIPconsIsStickingAtNode(gateconss[i]);

		  SCIPdebug( SCIP_CALL( SCIPprintCons(scip, gateconss[i], NULL) ) );

		  SCIP_CALL( SCIPdelCons(scip, gateconss[i]) );
		  ++(*ndelconss);
	       }

	       SCIPdebug( SCIP_CALL( SCIPprintCons(scip, logicor, NULL) ) );

	       /* create and add "and" constraint for the extracted gate */
	       (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "setpart_%d", presoldata->ngates);
	       SCIP_CALL( SCIPcreateConsSetpart(scip, &newcons, name, 3, activevars,
		     initial, separate, enforce, check, propagate,
		     local, modifiable, dynamic, removable, stickingatnode) );

	       SCIP_CALL( SCIPaddCons(scip, newcons) );
	       SCIPdebugMessage("-------------->\n");
	       SCIPdebug( SCIP_CALL( SCIPprintCons(scip, newcons, NULL) ) );

	       ++(*naddconss);
	       ++(presoldata->ngates);

	       SCIP_CALL( SCIPdelCons(scip, logicor) );
	       ++(*ndelconss);

	       /* @todo: maybe remove the deleted logicor from the hashmap too */

	       SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
	    }
	    /* did we find two set-packing constraints for the upgrade to the gate constraint */
	    else if( !presoldata->onlysetpart && nfound == 2 )
	    {
	       SCIP_VAR* vars[2];
	       SCIP_VAR* resvar = NULL;
	       SCIP_CONS* newcons;
	       char name[SCIP_MAXSTRLEN];
	       SCIP_Bool initial;
	       SCIP_Bool separate;
	       SCIP_Bool enforce;
	       SCIP_Bool check;
	       SCIP_Bool propagate;
	       SCIP_Bool local;
	       SCIP_Bool modifiable;
	       SCIP_Bool dynamic;
	       SCIP_Bool removable;
	       SCIP_Bool stickingatnode;
	       int i;

	       if( gateconss[0] != NULL )
	       {
		  if( gateconss[1] != NULL )
		  {
		     assert(gateconss[2] == NULL);

		     resvar = activevars[0];
		     SCIP_CALL( SCIPgetNegatedVar(scip, activevars[1], &vars[0]) );
		     SCIP_CALL( SCIPgetNegatedVar(scip, activevars[2], &vars[1]) );
		  }
		  else if( gateconss[2] != NULL )
		  {
		     assert(gateconss[1] == NULL);

		     resvar = activevars[1];
		     SCIP_CALL( SCIPgetNegatedVar(scip, activevars[0], &vars[0]) );
		     SCIP_CALL( SCIPgetNegatedVar(scip, activevars[2], &vars[1]) );
		  }
	       }
	       else if( gateconss[1] != NULL && gateconss[2] != NULL )
	       {
		  assert(gateconss[0] == NULL);

		  resvar = activevars[2];
		  SCIP_CALL( SCIPgetNegatedVar(scip, activevars[0], &vars[0]) );
		  SCIP_CALL( SCIPgetNegatedVar(scip, activevars[1], &vars[1]) );
	       }
	       assert(resvar != NULL);

	       initial = SCIPconsIsInitial(logicor);
	       separate = SCIPconsIsSeparated(logicor);
	       enforce = SCIPconsIsEnforced(logicor);
	       check = SCIPconsIsChecked(logicor);
	       propagate = SCIPconsIsPropagated(logicor);
	       local = SCIPconsIsLocal(logicor);
	       modifiable = SCIPconsIsModifiable(logicor);
	       dynamic = SCIPconsIsDynamic(logicor);
	       removable = SCIPconsIsRemovable(logicor);
	       stickingatnode = SCIPconsIsStickingAtNode(logicor);

	       SCIPdebugMessage("Following three constraints form an and constraint.\n");

	       for( i = 2; i >= 0; --i )
	       {
		  if( gateconss[i] != NULL )
		  {
		     initial |= SCIPconsIsInitial(gateconss[i]);
		     separate |= SCIPconsIsSeparated(gateconss[i]);
		     enforce |= SCIPconsIsEnforced(gateconss[i]);
		     check |= SCIPconsIsChecked(gateconss[i]);
		     propagate |= SCIPconsIsPropagated(gateconss[i]);
		     local &= SCIPconsIsLocal(gateconss[i]);
		     modifiable &= SCIPconsIsModifiable(gateconss[i]);
		     dynamic &= SCIPconsIsDynamic(gateconss[i]);
		     removable &= SCIPconsIsRemovable(gateconss[i]);
		     stickingatnode &= SCIPconsIsStickingAtNode(gateconss[i]);

		     SCIPdebug( SCIP_CALL( SCIPprintCons(scip, gateconss[i], NULL) ) );
		     SCIP_CALL( SCIPdelCons(scip, gateconss[i]) );
		     ++(*ndelconss);
		  }
	       }
	       SCIPdebug( SCIP_CALL( SCIPprintCons(scip, logicor, NULL) ) );

	       /* create and add "and" constraint for the extracted gate */
	       (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "andgate_%d", presoldata->ngates);
	       SCIP_CALL( SCIPcreateConsAnd(scip, &newcons, name, resvar, 2, vars,
		     initial, separate, enforce, check, propagate,
		     local, modifiable, dynamic, removable, stickingatnode) );

	       SCIP_CALL( SCIPaddCons(scip, newcons) );
	       SCIPdebugMessage("-------------->\n");
	       SCIPdebug( SCIP_CALL( SCIPprintCons(scip, newcons, NULL) ) );

	       ++(*naddconss);
	       ++(presoldata->ngates);


	       SCIP_CALL( SCIPdelCons(scip, logicor) );
	       ++(*ndelconss);

	       /* @todo: maybe remove the deleted logicor from the hashmap too */

	       SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
	    }
	 }
      }
      SCIPhashmapFree(&varmap);
   }

 TERMINATE:
   SCIPfreeBufferArray(scip, &setppcconss);

   /* remove old setppchashdatas objects */
   SCIP_CALL( cleanupHashDatas(scip, presoldata) );

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the gateextraction presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolGateextraction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* alloc presolve data object */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   /* initialize gateextraction presolver data */
   SCIP_CALL( presoldataInit(scip, presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolCopyGateextraction,
         presolFreeGateextraction, presolInitGateextraction, presolExitGateextraction,
         presolInitpreGateextraction, presolExitpreGateextraction, presolExecGateextraction,
         presoldata) );

   /* add gateextraction presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/"PRESOL_NAME"/onlysetpart",
         "should we only try to extract set-partitioning constraints and no and-constraints",
         &presoldata->onlysetpart, TRUE, DEFAULT_ONLYSETPART, NULL, NULL) );

   return SCIP_OKAY;
}
