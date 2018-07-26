/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_orbitalfixing.c
 * @brief  propagator for orbital fixing
 * @author Marc Pfetsch
 *
 * This propagator implements orbital fixing as introduced by
 *
 * F. Margot: Exploiting orbits in symmetric ILP. Math. Program., 98(1-3):3–21, 2003.
 *
 * The method obtains symmetries from the symmetry presolver and then computes orbits of variables with respect to the
 * subgroup of the symmetry group that stabilizes the variables branched to 1. Then one can fix all variables in an
 * orbit to 0 or 1 if one of the other variables in the orbit is fixed to 0 or 1, respectively. Different from Margot,
 * the subgroup is obtained by filtering out generators that do not individually stabilize the variables branched to 1.
 *
 * @todo Possibly turn off propagator in subtrees.
 * @todo Check application of conflict resolution.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "prop_orbitalfixing.h"

#include <scip/pub_tree.h>
#include <scip/pub_table.h>

#include "presol_symmetry.h"
#include "presol_symbreak.h"

/* propagator properties */
#define PROP_NAME              "orbitalfixing"
#define PROP_DESC              "propagator for orbital fixing"
#define PROP_TIMING    SCIP_PROPTIMING_BEFORELP   /**< propagation timing mask */
#define PROP_PRIORITY          -1000000           /**< propagator priority */
#define PROP_FREQ                     1           /**< propagator frequency */
#define PROP_DELAY                FALSE           /**< should propagation method be delayed, if other propagators found reductions? */

#define PROP_PRESOL_PRIORITY    -1000000          /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers) */
#define PROP_PRESOLTIMING   SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1           /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */

/* parameters */
#define DEFAULT_SYMCOMPTIMING         2           /**< timing of symmetry computation for orbital fixing
                                                   *   (0 = before presolving, 1 = during presolving, 2 = at first call) */
#define DEFAULT_PERFORMPRESOLVING     FALSE       /**< Run orbital fixing during presolving? */
#define DEFAULT_ENABLEDAFTERRESTARTS  FALSE       /**< Run orbital fixing after a restart has occured? */


/* output table properties */
#define TABLE_NAME_ORBITALFIXING        "orbitalfixing"
#define TABLE_DESC_ORBITALFIXING        "orbital fixing statistics"
#define TABLE_POSITION_ORBITALFIXING    7001                    /**< the position of the statistics table */
#define TABLE_EARLIEST_ORBITALFIXING    SCIP_STAGE_SOLVING      /**< output of the statistics table is only printed from this stage onwards */


/*
 * Data structures
 */

/** propagator data for orbital branching */
struct SCIP_PropData
{
   int                   npermvars;          /**< pointer to store number of variables for permutations */
   SCIP_VAR**            permvars;           /**< pointer to store variables on which permutations act */
   SCIP_HASHMAP*         permvarmap;         /**< map of variables to indices in permvars array */
   int                   nperms;             /**< pointer to store number of permutations */
   int**                 perms;              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Bool             enabled;            /**< run orbital branching? */
   SCIP_Bool             performpresolving;  /**< Run orbital fixing during presolving? */
   SCIP_Bool             enabledafterrestarts; /**< Run orbital fixing after a restart has occured? */
   int                   symcomptiming;      /**< timing of symmetry computation for orbital fixing
                                              *   (0 = before presolving, 1 = during presolving, 2 = at first call) */
   int                   lastrestart;        /**< last restart for which symmetries have been computed */
   int                   nfixedzero;         /**< number of variables fixed to 0 */
   int                   nfixedone;          /**< number of variables fixed to 1 */
   SCIP_Longint          nodenumber;         /**< number of node where propagation has been last applied */
};



/*
 * Table callback methods
 */

/** table data */
struct SCIP_TableData
{
   SCIP_PROPDATA*        propdata;           /** pass data of propagator for table output function */
};


/** output method of orbital fixing propagator statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputOrbitalfixing)
{
   SCIP_TABLEDATA* tabledata;

   assert( scip != NULL );
   assert( table != NULL );

   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );
   assert( tabledata->propdata != NULL );

   if ( tabledata->propdata->nperms > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "Orbital fixing     :\n");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  vars fixed to 0  :%11d\n", tabledata->propdata->nfixedzero);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, file, "  vars fixed to 1  :%11d\n", tabledata->propdata->nfixedone);
   }

   return SCIP_OKAY;
}


/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeOrbitalfixing)
{
   SCIP_TABLEDATA* tabledata;
   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   SCIPfreeBlockMemory(scip, &tabledata);

   return SCIP_OKAY;
}

/*
 * Local methods
 */


/** possibly get symmetries */
static
SCIP_RETCODE getSymmetries(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int v;

   assert( scip != NULL );
   assert( propdata != NULL );

   if ( propdata->nperms < 0 || SCIPgetNRuns(scip) > propdata->lastrestart )
   {
      /* recompute symmetries after a restart */
      if ( SCIPgetNRuns(scip) > propdata->lastrestart )
      {
         /* reset symmetry information */
         if ( propdata->permvarmap != NULL )
         {
            SCIPhashmapFree(&propdata->permvarmap);
         }
         propdata->nperms = -1;
         propdata->perms = NULL;
         propdata->permvars = NULL;
         propdata->npermvars = -1;
         propdata->permvarmap = NULL;

         /* recompute symmetries and update restart counter */
         SCIP_CALL( SCIPgetGeneratorsSymmetry(scip, SYM_SPEC_BINARY, SYM_SPEC_INTEGER, TRUE,
               &(propdata->npermvars), &(propdata->permvars), &(propdata->nperms), &(propdata->perms), NULL, NULL) );

         propdata->lastrestart = SCIPgetNRuns(scip);
      }
      else
      {
         SCIP_CALL( SCIPgetGeneratorsSymmetry(scip, SYM_SPEC_BINARY, SYM_SPEC_INTEGER, FALSE,
               &(propdata->npermvars), &(propdata->permvars), &(propdata->nperms), &(propdata->perms), NULL, NULL) );
      }

      if ( propdata->nperms == 0 )
         return SCIP_OKAY;

      /* create hashmap for storing the indices of variables */
      assert( propdata->permvarmap == NULL );
      SCIP_CALL( SCIPhashmapCreate(&propdata->permvarmap, SCIPblkmem(scip), propdata->npermvars) );

      /* insert variables */
      for (v = 0; v < propdata->npermvars; ++v)
      {
         SCIP_CALL( SCIPhashmapInsert(propdata->permvarmap, propdata->permvars[v], (void*) (size_t) v) );
      }
   }

   return SCIP_OKAY;
}

/** perform orbital fixing
 *
 *  Note that we do not have to distinguish between variables that have been fixed or branched to 1, since the
 *  stabilizer is with respect to the variables that have been branched to 1. Thus, if an orbit contains a variable that
 *  has been branched to 1, the whole orbit only contains variables that have been branched to 1 - and nothing can be
 *  fixed.
 */
static
SCIP_RETCODE performOrbitalFixing(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_VAR**            permvars,           /**< variables */
   int                   npermvars,          /**< number of variables */
   int*                  orbits,             /**< array of non-trivial orbits */
   int*                  orbitbegins,        /**< array containing begin positions of new orbits in orbits array */
   int                   norbits,            /**< number of orbits */
   SCIP_Bool*            infeasible,         /**< pointer to store whether problem is infeasible */
   int*                  nfixedzero,         /**< pointer to store number of variables fixed to 0 */
   int*                  nfixedone           /**< pointer to store number of variables fixed to 1 */
   )
{
   SCIP_Bool tightened;
   int i;

   assert( scip != NULL );
   assert( permvars != NULL );
   assert( orbits != NULL );
   assert( orbitbegins != NULL );
   assert( infeasible != NULL );
   assert( nfixedzero != NULL );
   assert( nfixedone != NULL );
   assert( norbits > 0 );
   assert( orbitbegins[0] == 0 );

   *infeasible = FALSE;
   *nfixedzero = 0;
   *nfixedone = 0;

   SCIPdebugMsg(scip, "Perform orbital fixing on %d orbits.\n", norbits);

   /* check all orbits */
   for (i = 0; i < norbits; ++i)
   {
      SCIP_Bool havefixedone = FALSE;
      SCIP_Bool havefixedzero = FALSE;
      SCIP_VAR* var;
      int j;

      /* we only have nontrivial orbits */
      assert( orbitbegins[i+1] - orbitbegins[i] >= 2 );

      /* check all variables in the orbit */
      for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
      {
         assert( 0 <= orbits[j] && orbits[j] < npermvars );
         var = permvars[orbits[j]];
         assert( var != NULL );

         /* check whether variable is not binary (and not implicit integer!) */
         if ( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
         {
            /* skip orbit if there are non-binary variables */
            havefixedone = FALSE;
            havefixedzero = FALSE;
            break;
         }

         /* if variable is fixed to 1 -> can fix all variables in orbit to 1 */
         if ( SCIPvarGetLbLocal(var) > 0.5 )
            havefixedone = TRUE;

         /* check for zero-fixed variables */
         if ( SCIPvarGetUbLocal(var) < 0.5 )
            havefixedzero = TRUE;
      }

      /* check consistency */
      if ( havefixedone && havefixedzero )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }

      /* fix all variables to 0 if there is one variable fixed to 0 */
      if ( havefixedzero )
      {
         assert( ! havefixedone );

         for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         {
            assert( 0 <= orbits[j] && orbits[j] < npermvars );
            var = permvars[orbits[j]];
            assert( var != NULL );

            /* only variables that are not yet fixed to 0 */
            if ( SCIPvarGetUbLocal(var) > 0.5 )
            {
               SCIPdebugMsg(scip, "can fix <%s> (index %d) to 0.\n", SCIPvarGetName(var), orbits[j]);
               assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
               /* due to aggregation, var might already be fixed to 1, so do not put assert here */

               /* do not use SCIPinferBinvarProp(), since conflict analysis is not valid */
               SCIP_CALL( SCIPtightenVarUb(scip, var, 0.0, FALSE, infeasible, &tightened) );
               if ( *infeasible )
                  return SCIP_OKAY;
               if ( tightened )
                  ++(*nfixedzero);
            }
         }
      }

      /* fix all variables to 1 if there is one variable fixed to 1 */
      if ( havefixedone )
      {
         assert( ! havefixedzero );

         for (j = orbitbegins[i]; j < orbitbegins[i+1]; ++j)
         {
            assert( 0 <= orbits[j] && orbits[j] < npermvars );
            var = permvars[orbits[j]];
            assert( var != NULL );

            /* only variables that are not yet fixed to 1 */
            if ( SCIPvarGetLbLocal(var) < 0.5)
            {
               SCIPdebugMsg(scip, "can fix <%s> (index %d) to 1.\n", SCIPvarGetName(var), orbits[j]);
               assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
               /* due to aggregation, var might already be fixed to 0, so do not put assert here */

               /* do not use SCIPinferBinvarProp(), since conflict analysis is not valid */
               SCIP_CALL( SCIPtightenVarLb(scip, var, 1.0, FALSE, infeasible, &tightened) );
               if ( *infeasible )
                  return SCIP_OKAY;
               if ( tightened )
                  ++(*nfixedone);
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** Get branching variables on the path to root */
static
SCIP_RETCODE computeBranchingVariables(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varmap,             /**< map of variables to indices in vars array */
   SCIP_Shortbool*       b1,                 /**< bitset marking the variables branched to 1 */
   SCIP_Bool*            success             /**< pointer to store whether branching variables were computed successfully */
   )
{
   SCIP_NODE* node;

   assert( scip != NULL );
   assert( varmap != NULL );
   assert( b1 != NULL );
   assert( success != NULL );

   *success = TRUE;

   /* get current node */
   node = SCIPgetCurrentNode(scip);

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPprintNodeRootPath(scip, node, NULL) );
#endif

   /* follow path to the root (in the root no domains were changed due to branching) */
   while ( SCIPnodeGetDepth(node) != 0 )
   {
      SCIP_BOUNDCHG* boundchg;
      SCIP_DOMCHG* domchg;
      SCIP_VAR* branchvar;
      int nboundchgs;
      int i;

      /* get domain changes of current node */
      domchg = SCIPnodeGetDomchg(node);
      assert( domchg != NULL );

      /* loop through all bound changes */
      nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
      for (i = 0; i < nboundchgs; ++i)
      {
         /* get bound change info */
         boundchg = SCIPdomchgGetBoundchg(domchg, i);
         assert( boundchg != NULL );

         /* branching decisions have to be in the beginning of the bound change array */
         if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
            break;

         /* get corresponding branching variable */
         branchvar = SCIPboundchgGetVar(boundchg);

         /* we only consider binary variables */
         if ( SCIPvarGetType(branchvar) == SCIP_VARTYPE_BINARY )
         {
            /* make sure that branching variable is known, since new binary variables may have
             * been created meanwhile, e.g., by presol_inttobinary */
            if ( ! SCIPhashmapExists(varmap, (void*) branchvar) )
            {
               *success = FALSE;
               return SCIP_OKAY;
            }

            if ( SCIPvarGetLbLocal(branchvar) > 0.5 )
            {
               int branchvaridx;

               branchvaridx = (int) (size_t) SCIPhashmapGetImage(varmap, (void*) branchvar);
               assert( branchvaridx < nvars );
               b1[branchvaridx] = TRUE;
            }
         }
      }

      node = SCIPnodeGetParent(node);
   }

   return SCIP_OKAY;
}


/** propagate orbital fixing */
static
SCIP_RETCODE propagateOrbitalFixing(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the node is detected to be infeasible */
   int*                  nprop               /**< pointer to store the number of propagations */
   )
{
   SCIP_Shortbool* activeperms;
   SCIP_Shortbool* b1;
   SCIP_Bool success = TRUE;
   SCIP_VAR** permvars;
   int* orbitbegins;
   int* orbits;
#ifndef NDEBUG
   SCIP_Real* permvarsobj;
#endif
   int norbits;
   int npermvars;
   int** perms;
   int nperms;
   int p;
   int v;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( infeasible != NULL );
   assert( nprop != NULL );

   *infeasible = FALSE;
   *nprop = 0;

   /* possibly get symmetries */
   SCIP_CALL( getSymmetries(scip, propdata) );

   permvars = propdata->permvars;
   npermvars = propdata->npermvars;
   perms = propdata->perms;
   nperms = propdata->nperms;

   /* return if there is no symmetry available */
   if ( nperms == 0 )
      return SCIP_OKAY;

   assert( propdata->permvars != NULL );
   assert( propdata->npermvars > 0 );
   assert( propdata->permvarmap != NULL );
   assert( propdata->perms != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &activeperms, nperms) );

   /* init bitset for marking variables branched to 1 */
   SCIP_CALL( SCIPallocBufferArray(scip, &b1, npermvars) );
   for (v = 0; v < npermvars; ++v)
      b1[v] = FALSE;

   /* get branching variables */
   SCIP_CALL( computeBranchingVariables(scip, npermvars, propdata->permvarmap, b1, &success) );

   if ( ! success )
   {
      SCIPfreeBufferArray(scip, &b1);
      SCIPfreeBufferArray(scip, &activeperms);
      return SCIP_OKAY;
   }

#ifndef NDEBUG
   SCIP_CALL( SCIPgetPermvarsObjSymmetry(scip, &permvarsobj) );
#endif
   assert( permvarsobj != NULL );

   /* filter out permutations that move variables that are fixed to different values */
   for (p = 0; p < nperms; ++p)
   {
      assert( perms[p] != NULL );

      for (v = 0; v < npermvars; ++v)
      {
         int img;

         img = perms[p][v];

         if ( img != v )
         {
#ifndef NDEBUG
            SCIP_VAR* varv = permvars[v];
            SCIP_VAR* varimg = permvars[img];

            /* check whether moved variables have the same type (might have been aggregated in the meanwhile) */
            assert( SCIPvarGetType(varv) == SCIPvarGetType(varimg) ||
               (SCIPvarIsBinary(varv) && SCIPvarIsBinary(varimg)) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_IMPLINT && SCIPvarGetType(varimg) == SCIP_VARTYPE_CONTINUOUS &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) ||
               (SCIPvarGetType(varv) == SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(varimg) == SCIP_VARTYPE_IMPLINT &&
                  SCIPisEQ(scip, SCIPvarGetLbGlobal(varv), SCIPvarGetLbGlobal(varimg)) &&
                  SCIPisEQ(scip, SCIPvarGetUbGlobal(varv), SCIPvarGetUbGlobal(varimg))) );
#endif
            assert( SCIPisEQ(scip, permvarsobj[v], permvarsobj[img]) );

            /* we are moving a variable branched to 1 to another variable */
            if ( b1[v] && ! b1[img] )
               break;

            /* Global variable fixings during the solving process might arise because parts of the tree are
             * pruned. Since these fixings might be caused by orbital fixing, they can be in conflict with
             * the symmetry handling decisions of orbital fixing in the part of the tree that is not pruned.
             * Thus, we have to take global fixings into account when filtering out symmetries.
             */
            if ( (SCIPvarGetLbGlobal(permvars[v]) > 0.5 && SCIPvarGetLbGlobal(permvars[img]) < 0.5) ||
               (SCIPvarGetLbGlobal(permvars[v]) < 0.5 && SCIPvarGetLbGlobal(permvars[img]) > 0.5) )
               break;
         }
      }

      if ( v >= npermvars )
         activeperms[p] = TRUE;
      else
         activeperms[p] = FALSE;
   }

   SCIPfreeBufferArray(scip, &b1);

   /* compute orbits */
   SCIP_CALL( SCIPallocBufferArray(scip, &orbits, npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &orbitbegins, npermvars) );
   SCIP_CALL( SCIPcomputeGroupOrbitsSymbreak(scip, permvars, npermvars, perms, nperms, activeperms, orbits, orbitbegins, &norbits) );

   SCIPfreeBufferArray(scip, &activeperms);

   if ( norbits > 0 )
   {
      int nfixedzero = 0;
      int nfixedone = 0;

      SCIP_CALL( performOrbitalFixing(scip, permvars, npermvars, orbits, orbitbegins, norbits, infeasible, &nfixedzero, &nfixedone) );

      propdata->nfixedzero += nfixedzero;
      propdata->nfixedone += nfixedone;
      *nprop = nfixedzero + nfixedone;

      SCIPdebugMsg(scip, "Orbital fixings: %d 0s, %d 1s.\n", nfixedzero, nfixedone);
   }

   SCIPfreeBufferArray(scip, &orbitbegins);
   SCIPfreeBufferArray(scip, &orbits);

   return SCIP_OKAY;
}




/*
 * Callback methods of propagator
 */


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeOrbitalfixing)
{  /*lint --e{715,818}*/
   SCIP_PROPDATA* propdata;

   assert( prop != NULL );

   SCIPdebugMsg(scip, "Freeing propagator <%s> ...\n", SCIPpropGetName(prop));

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIPfreeBlockMemory(scip, &propdata);

   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitOrbitalfixing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int usesymmetry;

   assert( prop != NULL );

   SCIPdebugMsg(scip, "Init propagator <%s> ...\n", SCIPpropGetName(prop));

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* check whether we should run */
   SCIP_CALL( SCIPgetIntParam(scip, "misc/usesymmetry", &usesymmetry) );
   if ( usesymmetry == (int) SYM_HANDLETYPE_ORBITALFIXING )
      propdata->enabled = TRUE;
   else
   {
      propdata->enabled = FALSE;
   }

   return SCIP_OKAY;
}


/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitOrbitalfixing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   if ( propdata->permvarmap != NULL )
   {
      SCIPhashmapFree(&propdata->permvarmap);
   }

   /* reset propagator variables */
   propdata->nodenumber = -1;
   propdata->nfixedzero = 0;
   propdata->nfixedone = 0;

   propdata->nperms = -1;
   propdata->perms = NULL;
   propdata->permvars = NULL;
   propdata->npermvars = -1;
   propdata->permvarmap = NULL;
   propdata->lastrestart = 0;

   return SCIP_OKAY;
}


/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreOrbitalfixing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* possibly skip orbital fixing */
   if ( ! propdata->enabled ||  propdata->nperms == 0 )
      return SCIP_OKAY;

   /* stop, if problem has already been solved */
   if ( SCIPgetStatus(scip) != SCIP_STATUS_UNKNOWN )
      return SCIP_OKAY;

   /* run only if timing is correct */
   assert( 0 <= propdata->symcomptiming && propdata->symcomptiming <= 2 );
   if ( propdata->symcomptiming > 0 )
      return SCIP_OKAY;

   assert( SCIPisTransformed(scip) );

   /* possibly get symmetries */
   SCIP_CALL( getSymmetries(scip, propdata) );

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolOrbitalFixing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible = FALSE;
   int nprop = 0;

   assert( scip != NULL );
   assert( nfixedvars != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* check whether we run after a restart */
   if ( propdata->enabled && ! propdata->enabledafterrestarts && SCIPgetNRuns(scip) > 1 )
      propdata->enabled = FALSE;

   /* do not run if not enabled */
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   /* run only if timing is correct */
   assert( 0 <= propdata->symcomptiming && propdata->symcomptiming <= 2 );
   if ( propdata->symcomptiming > 1 )
      return SCIP_OKAY;

   /* run if presolving should be performed */
   if ( propdata->performpresolving )
   {
      /* propagate */
      *result = SCIP_DIDNOTFIND;

      SCIPdebugMsg(scip, "Presolving <%s>.\n", SCIPpropGetName(prop));
      SCIP_CALL( propagateOrbitalFixing(scip, propdata, &infeasible, &nprop) );
      if ( infeasible )
         *result = SCIP_CUTOFF;
      else if ( nprop > 0 )
      {
         *result = SCIP_SUCCESS;
         *nfixedvars += nprop;
      }
   }
   else if ( propdata->symcomptiming == 1 )
   {
      /* otherwise compute symmetry if timing requests it */
      SCIP_CALL( getSymmetries(scip, propdata) );
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecOrbitalfixing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible = FALSE;
   SCIP_Longint nodenumber;
   int nprop = 0;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* do nothing if we are in a probing node */
   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   /* check whether we run after a restart */
   if ( propdata->enabled && ! propdata->enabledafterrestarts && SCIPgetNRuns(scip) > 1 )
      propdata->enabled = FALSE;

   /* do not run if not enabled */
   if ( ! propdata->enabled )
      return SCIP_OKAY;

   /* return if there is no symmetry available */
   if ( propdata->nperms == 0 )
      return SCIP_OKAY;

   /* return if we already ran in this node */
   nodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   if ( nodenumber == propdata->nodenumber )
      return SCIP_OKAY;
   propdata->nodenumber = nodenumber;

   /* propagate */
   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "Propagating <%s>.\n", SCIPpropGetName(prop));
   SCIP_CALL( propagateOrbitalFixing(scip, propdata, &infeasible, &nprop) );
   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( nprop > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator
 *
 *  @todo Implement reverse propagation.
 *
 *  Note that this is relatively difficult to obtain: One needs to include all bounds of variables to would lead to a
 *  different orbit in which the variables that was propagated lies. This includes all variables that are moved by the
 *  permutations which are involved in creating the orbit.
 */
static
SCIP_DECL_PROPRESPROP(propRespropOrbitalfixing)
{  /*lint --e{715,818}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** creates the orbitalfixing propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropOrbitalfixing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create orbitalfixing propagator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );

   propdata->nodenumber = -1;
   propdata->nfixedzero = 0;
   propdata->nfixedone = 0;

   propdata->nperms = -1;
   propdata->perms = NULL;
   propdata->permvars = NULL;
   propdata->npermvars = -1;
   propdata->permvarmap = NULL;
   propdata->lastrestart = 0;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, propExecOrbitalfixing, propdata) );

   /* set callbacks */
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeOrbitalfixing) );
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitOrbitalfixing) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitOrbitalfixing) );
   SCIP_CALL( SCIPsetPropInitpre(scip, prop, propInitpreOrbitalfixing) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropOrbitalfixing) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolOrbitalFixing, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );

   /* include table */
   SCIP_CALL( SCIPallocBlockMemory(scip, &tabledata) );
   tabledata->propdata = propdata;
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME_ORBITALFIXING, TABLE_DESC_ORBITALFIXING, TRUE,
         NULL, tableFreeOrbitalfixing, NULL, NULL, NULL, NULL, tableOutputOrbitalfixing,
         tabledata, TABLE_POSITION_ORBITALFIXING, TABLE_EARLIEST_ORBITALFIXING) );

   /* add parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/" PROP_NAME "/symcomptiming",
         "timing of symmetry computation for orbital fixing (0 = before presolving, 1 = during presolving, 2 = at first call)",
         &propdata->symcomptiming, TRUE, DEFAULT_SYMCOMPTIMING, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/performpresolving",
         "Run orbital fixing during presolving?",
         &propdata->performpresolving, TRUE, DEFAULT_PERFORMPRESOLVING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/" PROP_NAME "/enabledafterrestarts",
         "Run orbital fixing after a restart has occured?",
         &propdata->enabledafterrestarts, TRUE, DEFAULT_ENABLEDAFTERRESTARTS, NULL, NULL) );

   return SCIP_OKAY;
}
