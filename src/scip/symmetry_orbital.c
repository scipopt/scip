/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   symmetry_orbital.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling symmetries by orbital fixing
 * @author Jasper van Doornmalen
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/symmetry_orbital.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/struct_var.h"
#include "scip/type_var.h"
#include "scip/scip.h"
#include "scip/scip_branch.h"
#include "scip/scip_conflict.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"
#include "scip/debug.h"
#include "scip/struct_scip.h"
#include "scip/struct_mem.h"
#include "scip/struct_tree.h"
#include "scip/symmetry.h"
#include "scip/event_shadowtree.h"
#include <ctype.h>
#include <string.h>
#include <symmetry/type_symmetry.h>

#include <memory.h>


/* event handler properties */
#define EVENTHDLR_SYMMETRY_NAME    "symmetry_orbital"
#define EVENTHDLR_SYMMETRY_DESC    "filter global variable fixing event handler for orbital fixing"


/*
 * Data structures
 */


/** data for dynamic orbital fixing component propagator */
struct OrbitalFixingComponentData
{
   SCIP_NODE*            lastnode;           /**< last node processed by dynamic orbital fixing component */
   SCIP_Real*            globalvarlbs;       /**< global variable lower bounds, determined when adding symmetry cons */
   SCIP_Real*            globalvarubs;       /**< global variable upper bounds, determined when adding symmetry cons */
   int**                 perms;              /**< the permutations for orbital fixing */
   int                   nperms;             /**< the number of permutations in perms */
   SCIP_VAR**            permvars;           /**< array consisting of the variables of this component */
   int                   npermvars;          /**< number of vars in this component */
   SCIP_HASHMAP*         permvarmap;         /**< map of variables to indices in permvars array */
};
typedef struct OrbitalFixingComponentData OFCDATA;


/** data for dynamic orbital fixing propagator */
struct SCIP_OrbitalFixingData
{
   SCIP_EVENTHDLR*       shadowtreeeventhdlr;/**< eventhandler for the shadow tree data structure */
   SCIP_EVENTHDLR*       globalfixeventhdlr; /**< event handler for handling global variable fixings */

   OFCDATA**             componentdatas;     /**< array of pointers to individual components for orbital fixing */
   int                   ncomponents;        /**< number of orbital fixing datas in array */
   int                   maxncomponents;     /**< allocated orbital fixing datas array size */
};


/** data structures for provisional disjoint set data structures */

/** provisional disjoint set with minimum/maximum information transaction types */
enum PrDjSetWMinTransaction_Type {
   PRDJSETWMINTRANSACTION_TYPE_PARENT     = 0,
   PRDJSETWMINTRANSACTION_TYPE_SIZE       = 1,
   PRDJSETWMINTRANSACTION_TYPE_LBVALUEMIN = 2,
   PRDJSETWMINTRANSACTION_TYPE_LBVALUEMAX = 3,
   PRDJSETWMINTRANSACTION_TYPE_UBVALUEMIN = 4
};
typedef enum PrDjSetWMinTransaction_Type PRDJSETWMINTRANSACTION_TYPE;

/** provisional disjoint set with minimum/maximum information transaction
 * 
 * The transaction allows for resetting the disjoint set data structure to a former state, before applying unions.
 */
struct PrDjSetWMinTransaction
{
   PRDJSETWMINTRANSACTION_TYPE type;         /**< type of data stored */
   int                   index;              /**< array index of transaction */
   union 
   {
      SCIP_Real          valuereal;          /**< value, if it's a 'real' type */
      int                valueint;           /**< value, if it's an 'int' type */
   }                     data;               /**< data object, either vlauereal or valueint. */
};
typedef struct PrDjSetWMinTransaction PRDJSETWMINTRANSACTION;

/** provisional disjoint set with minimum/maximum information
 * 
 * This is a standard disjoint set data structure, but it also stores certain minima and maxima over arrays in each
 * component. Besides, it is 'provisional' in the sense that performed steps can be undone.
 */
struct PrDjSetWMin
{
   int*                  parents;            /**< array to store the parent node index for every vertex */
   int*                  sizes;              /**< array to store the size of the subtree rooted at each vertex */
   SCIP_Real*            lbminvalues;        /**< array to store the minimal values of the components, for var lbs */
   SCIP_Real*            lbmaxvalues;        /**< array to store the maximal values of the components, for var lbs */
   SCIP_Real*            ubminvalues;        /**< array to store the minimal values of the components, for var ubs */
   PRDJSETWMINTRANSACTION* transactions;     /**< array to store the transactions */
   int                   size;               /**< the number of vertices in the graph */
   int                   ntransactions;      /**< number of unconfirmed transactions */
   int                   maxntransactions;   /**< maximal number of transactions */
};
typedef struct PrDjSetWMin PRDJSETWMIN;


/*
 * Local methods
 */

/** helper functions for LT, GE, LE, GE, EQ, that do take infinity-values into account */
#if SCIP_DISABLED_CODE
static
SCIP_Bool EQ(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return TRUE;
   if ( inf1 != inf2 )
      return FALSE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return TRUE;
   if ( minf1 != minf2 )
      return FALSE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisEQ(scip, val1, val2);
}
#endif

static
SCIP_Bool LE(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return TRUE;
   if ( !inf1 && inf2 )
      return TRUE;
   if ( inf1 && !inf2 )
      return FALSE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return TRUE;
   if ( !minf1 && minf2 )
      return FALSE;
   if ( minf1 && !minf2 )
      return TRUE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisLE(scip, val1, val2);
}

static
SCIP_Bool GE(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return TRUE;
   if ( !inf1 && inf2 )
      return FALSE;
   if ( inf1 && !inf2 )
      return TRUE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return TRUE;
   if ( !minf1 && minf2 )
      return TRUE;
   if ( minf1 && !minf2 )
      return FALSE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisGE(scip, val1, val2);
}

static
SCIP_Bool LT(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return FALSE;
   if ( !inf1 && inf2 )
      return TRUE;
   if ( inf1 && !inf2 )
      return FALSE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return FALSE;
   if ( !minf1 && minf2 )
      return FALSE;
   if ( minf1 && !minf2 )
      return TRUE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisLT(scip, val1, val2);
}

static
SCIP_Bool GT(SCIP* scip, SCIP_Real val1, SCIP_Real val2)
{
   SCIP_Bool inf1;
   SCIP_Bool inf2;
   SCIP_Bool minf1;
   SCIP_Bool minf2;

   inf1 = SCIPisInfinity(scip, val1);
   inf2 = SCIPisInfinity(scip, val2);
   if ( inf1 && inf2 )
      return FALSE;
   if ( !inf1 && inf2 )
      return FALSE;
   if ( inf1 && !inf2 )
      return TRUE;
   assert( !inf1 );
   assert( !inf2 );

   minf1 = SCIPisInfinity(scip, -val1);
   minf2 = SCIPisInfinity(scip, -val2);
   if ( minf1 && minf2 )
      return FALSE;
   if ( !minf1 && minf2 )
      return TRUE;
   if ( minf1 && !minf2 )
      return FALSE;
   assert( !minf1 );
   assert( !minf2 );

   return SCIPisGT(scip, val1, val2);
}


/** clears the provisional disjoint set (union find) with minimum structure \p djset */
static
void provisionalDisjointSetWithMinimumClear(
   PRDJSETWMIN*          djset,              /**< disjoint set (union find) data structure */
   SCIP_Real*            lbminvalues,        /**< values to maintain component minimum of, for lb minima */
   SCIP_Real*            lbmaxvalues,        /**< values to maintain component minimum of, for lb maxima */
   SCIP_Real*            ubminvalues         /**< values to maintain component minimum of, for ub minima */
   )
{
   int i;

   /* reset all components to be unconnected */
   for( i = 0; i < djset->size; i++ )
   {
      djset->parents[i] = i;
      djset->sizes[i] = 1;
      djset->lbminvalues[i] = lbminvalues[i];
      djset->lbmaxvalues[i] = lbmaxvalues[i];
      djset->ubminvalues[i] = ubminvalues[i];
   }
}

/** finds and returns the component identifier of this \p element */
static
SCIP_RETCODE provisionalDisjointSetWithMinimumFind(
   SCIP*                 scip,               /**< SCIP data structure */
   PRDJSETWMIN*          djset,              /**< disjoint set (union find) data structure */
   int                   element,            /**< element to be found */
   int*                  index,              /**< NULL, or to populate with the index to find */
   SCIP_Real*            lbminvalue,         /**< NULL, or to populate with the minimal value in component */
   SCIP_Real*            lbmaxvalue,         /**< NULL, or to populate with the maximal value in component */
   SCIP_Real*            ubminvalue          /**< NULL, or to populate with the minimal value in component */
   )
{
   int newelement;
   int root = element;
   int* parents;
   PRDJSETWMINTRANSACTION* tr;

   assert( scip != NULL );
   assert( djset != NULL );
   assert( element >= 0 );
   assert( element < djset->size );
   assert( index != NULL || lbminvalue != NULL || lbmaxvalue != NULL || ubminvalue != NULL );

   parents = djset->parents;

   /* find root of this element */
   while( root != parents[root] )
      root = parents[root];

   /* compress the path to make future queries faster */
   while( element != root )
   {
      newelement = parents[element];

      if ( djset->ntransactions + 1 > djset->maxntransactions )
      {
         int newsize = SCIPcalcMemGrowSize(scip, djset->ntransactions + 1);
         assert( newsize > djset->maxntransactions );

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &djset->transactions, djset->maxntransactions, newsize) );
         djset->maxntransactions = newsize;
      }

      /* store parent old data transaction */
      assert( djset->ntransactions < djset->maxntransactions );
      tr = &djset->transactions[djset->ntransactions++];
      tr->type = PRDJSETWMINTRANSACTION_TYPE_PARENT;
      tr->index = element;
      tr->data.valueint = parents[element];

      parents[element] = root;
      element = newelement;
   }

   if ( index != NULL )
      *index = root;
   if ( lbminvalue != NULL )
      *lbminvalue = djset->lbminvalues[root];
   if ( lbmaxvalue != NULL )
      *lbmaxvalue = djset->lbmaxvalues[root];
   if ( ubminvalue != NULL )
      *ubminvalue = djset->ubminvalues[root];

   return SCIP_OKAY;
}

/** merges the components containing the elements \p p and \p q */
static
SCIP_RETCODE provisionalDisjointSetWithMinimumUnion(
   SCIP*                 scip,               /**< SCIP data structure */
   PRDJSETWMIN*          djset,              /**< disjoint set (union find) data structure */
   int                   p,                  /**< first element */
   int                   q                   /**< second element */
   )
{
   int idp;
   int idq;
   int swap;
   SCIP_Real lbminvalp;
   SCIP_Real lbminvalq;
   SCIP_Real lbmaxvalp;
   SCIP_Real lbmaxvalq;
   SCIP_Real ubminvalp;
   SCIP_Real ubminvalq;
   PRDJSETWMINTRANSACTION* tr;

   assert(djset != NULL);
   assert(0 <= p);
   assert(0 <= q);
   assert(djset->size > p);
   assert(djset->size > q);

   SCIP_CALL( provisionalDisjointSetWithMinimumFind(scip, djset, p, &idp, &lbminvalp, &lbmaxvalp, &ubminvalp) );
   SCIP_CALL( provisionalDisjointSetWithMinimumFind(scip, djset, q, &idq, &lbminvalq, &lbmaxvalq, &ubminvalq) );

   /* if p and q lie in the same component, there is nothing to be done */
   if( idp == idq )
      return SCIP_OKAY;

   if ( djset->sizes[idp] >= djset->sizes[idq] )
   {
      /* swap idp and idq */
      swap = idp;
      idp = idq;
      idq = swap;
   }

   if ( djset->ntransactions + 5 > djset->maxntransactions )
   {
      int newsize = SCIPcalcMemGrowSize(scip, djset->ntransactions + 5);
      assert( newsize > djset->maxntransactions );

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &djset->transactions, djset->maxntransactions, newsize) );
      djset->maxntransactions = newsize;
   }

   /* update parents, sizes and values arrays, and store transaction. */
   assert( djset->ntransactions < djset->maxntransactions );
   tr = &djset->transactions[djset->ntransactions++];
   tr->type = PRDJSETWMINTRANSACTION_TYPE_PARENT;
   tr->index = idp;
   tr->data.valueint = djset->parents[idp];
   djset->parents[idp] = idq;

   assert( djset->ntransactions < djset->maxntransactions );
   tr = &djset->transactions[djset->ntransactions++];
   tr->type = PRDJSETWMINTRANSACTION_TYPE_SIZE;
   tr->index = idq;
   tr->data.valueint = djset->sizes[idq];
   djset->sizes[idq] += djset->sizes[idp];

   assert( djset->ntransactions < djset->maxntransactions );
   tr = &djset->transactions[djset->ntransactions++];
   tr->type = PRDJSETWMINTRANSACTION_TYPE_LBVALUEMIN;
   tr->index = idq;
   tr->data.valuereal = djset->lbminvalues[idq];
   djset->lbminvalues[idq] = MIN(lbminvalp, lbminvalq);

   assert( djset->ntransactions < djset->maxntransactions );
   tr = &djset->transactions[djset->ntransactions++];
   tr->type = PRDJSETWMINTRANSACTION_TYPE_LBVALUEMAX;
   tr->index = idq;
   tr->data.valuereal = djset->lbmaxvalues[idq];
   djset->lbmaxvalues[idq] = MAX(lbmaxvalp, lbmaxvalq);

   assert( djset->ntransactions < djset->maxntransactions );
   tr = &djset->transactions[djset->ntransactions++];
   tr->type = PRDJSETWMINTRANSACTION_TYPE_UBVALUEMIN;
   tr->index = idq;
   tr->data.valuereal = djset->ubminvalues[idq];
   djset->ubminvalues[idq] = MIN(ubminvalp, ubminvalq);

   return SCIP_OKAY;
}

/* undo all unconfirmed changes */
static
SCIP_RETCODE provisionalDisjointSetWithMinimumUndo(
   PRDJSETWMIN*          djset               /**< disjoint set (union find) data structure */
)
{
   PRDJSETWMINTRANSACTION* tr;

   while (djset->ntransactions > 0)
   {
      tr = &djset->transactions[--djset->ntransactions];
      switch (tr->type)
      {
         case PRDJSETWMINTRANSACTION_TYPE_PARENT:
            djset->parents[tr->index] = tr->data.valueint;
            break;
         case PRDJSETWMINTRANSACTION_TYPE_SIZE:
            djset->sizes[tr->index] = tr->data.valueint;
            break;
         case PRDJSETWMINTRANSACTION_TYPE_LBVALUEMIN:
            djset->lbminvalues[tr->index] = tr->data.valuereal;
            break;
         case PRDJSETWMINTRANSACTION_TYPE_LBVALUEMAX:
            djset->lbmaxvalues[tr->index] = tr->data.valuereal;
            break;
         case PRDJSETWMINTRANSACTION_TYPE_UBVALUEMIN:
            djset->ubminvalues[tr->index] = tr->data.valuereal;
            break;
         default:
            return SCIP_ERROR;
      }
   }
   return SCIP_OKAY;
}

/** confirm all unconfirmed changes */
static
void provisionalDisjointSetWithMinimumConfirm(
   PRDJSETWMIN*          djset               /**< disjoint set (union find) data structure */
)
{
   djset->ntransactions = 0;
}

/** initialized provisional disjoint set with minimum */
static
SCIP_RETCODE provisionalDisjointSetWithMinimumCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   PRDJSETWMIN**         djset,              /**< disjoint set (union find) data structure */
   SCIP_Real*            lbminvalues,        /**< array to compute minimum over, for var lower bounds */
   SCIP_Real*            lbmaxvalues,        /**< array to compute maximum over, for var lower bounds */
   SCIP_Real*            ubminvalues,        /**< array to compute minimum over, for var upper bounds */
   int                   ncomponents,        /**< number of components */
   int                   maxntransactions    /**< maximal number of transactions before confirmation */
)
{
   assert( scip != NULL );
   assert( djset != NULL );
   assert( ncomponents > 0 );

   /* allocate the necessary memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, djset) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*djset)->parents), ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*djset)->sizes), ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*djset)->lbminvalues), ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*djset)->lbmaxvalues), ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*djset)->ubminvalues), ncomponents) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*djset)->transactions), maxntransactions) );
   (*djset)->size = ncomponents;
   (*djset)->ntransactions = 0;
   (*djset)->maxntransactions = maxntransactions;

   /* clear the data structure */
   provisionalDisjointSetWithMinimumClear(*djset, lbminvalues, lbmaxvalues, ubminvalues);

   return SCIP_OKAY;
}

/** frees the disjoint set (union find) data structure */
static
void provisionalDisjointSetWithMinimumFree(
   SCIP*                 scip,               /**< SCIP data structure */
   PRDJSETWMIN**         djset               /**< pointer to disjoint set (union find) data structure */
   )
{
   PRDJSETWMIN* dsptr;

   assert(djset != NULL);
   assert(*djset != NULL);

   dsptr = *djset;

   SCIPfreeBlockMemoryArray(scip, &dsptr->transactions, dsptr->maxntransactions);
   SCIPfreeBlockMemoryArray(scip, &dsptr->ubminvalues, dsptr->size);
   SCIPfreeBlockMemoryArray(scip, &dsptr->lbmaxvalues, dsptr->size);
   SCIPfreeBlockMemoryArray(scip, &dsptr->lbminvalues, dsptr->size);
   SCIPfreeBlockMemoryArray(scip, &dsptr->sizes, dsptr->size);
   SCIPfreeBlockMemoryArray(scip, &dsptr->parents, dsptr->size);
   SCIPfreeBlockMemory(scip, djset);
}


/** populates chosenperms with a generating set of the symmetry group stabilizing the branching decisions */
static
SCIP_RETCODE orbitalFixingDynamicGetSymmetryStabilizerSubgroup(
   SCIP*              scip,               /**< pointer to SCIP data structure */
   OFCDATA*           ofdata,             /**< pointer to data for orbital fixing data */
   int**              chosenperms,        /**< pointer to permutations that are chosen */
   int*               nchosenperms,       /**< pointer to store the number of chosen permutations */
   SCIP_Real*         varlbs,             /**< array of ofdata->permvars variable LBs. If NULL, use SCIPvarGetLbLocal */
   SCIP_Real*         varubs,             /**< array of ofdata->permvars variable UBs. If NULL, use SCIPvarGetUbLocal */
   int*               branchedvarindices, /**< array of given branching decisions, in branching order */
   SCIP_Bool*         inbranchedvarindices, /**< array stating whether variable with index in ofdata->permvars is
                                              *  contained in the branching decisions. */
   int                nbranchedvarindices /**< number of branching decisions */
)
{
   int i;
   int p;
   int* perm;
   int varid;
   SCIP_Real orbitlbmin;
   SCIP_Real orbitlbmax;
   SCIP_Real orbitubmin;
   int varidimage;
   PRDJSETWMIN* orbitmindjset;
   SCIP_Bool arraysgenerated;

   assert( scip != NULL );
   assert( ofdata != NULL );
   assert( chosenperms != NULL );
   assert( nchosenperms != NULL );
   assert( (varlbs == NULL) == (varubs == NULL) );
   assert( branchedvarindices != NULL );
   assert( inbranchedvarindices != NULL );
   assert( nbranchedvarindices >= 0 );

   /* compute the orbit of the branching variable of a stabilizing symmetry subgroup at this point.
    *
    * We are interested in propagating the lexicographic order constraint
    * \f$\sigma(x) \succeq \sigma(\gamma(x))$ for all \f$\gamma$ in the symmetry group.
    * Here, \f$\sigma$ restricts the solution vector x only to the branching variables, in the branching variable order.
    *
    * We thus want to consider the set (**) of permutations \gamma from the symmetry group
    * with \f$\sigma(x) = \sigma(\gamma(x))$ for all feasible solution vectors \f$x$ satisfying the constraint above.
    *
    * A subset of this set of permutations is the following set:
    * All permutations \f$\gamma$ in the symmetry group with \f$\sigma(x) \leq \sigma(\gamma(x))$ (elementwise).
    *
    * Note that these sets do not define a group per se.
    * We show a construction creating a subgroup of the symmetry group, which is a subset of the set (**).
    * We do this based on the (strong) generating set that is provided, and greedily build the symmetry subgroup.
    *
    * Start with S = {} as the empty generating set. Obviously, that is a subset of (**).
    * For each strong generator \f$\gamma$, Test if the group generated by the set S with \f$\gamma$ is a subset of **.
    * If so, add \f$\gamma$ to S. Otherwise, or if it cannot be determined easily, continue with the next generator.
    *
    * A simple test is the following: Compute the orbits of the group generated by the set S with \f$\gamma$,
    * and for every branching variable test if its upper bound is smaller or equal to its minimal orbit lower bound.
    * If that is the case, the test passes.
    *
    * It is easy to verify that this test accepts as least as many permutations as the readily implemented
    * orbital fixing method for binary variables. As such, a larger and still valid symmetry subgroup is detected.
    *
    * The steps above assume that all reductions found so far are compatible with this symmetry reduction method.
    * In other words, if a constraint/propagator that is not aware of the symmetry handling finds a reduction, then
    * this reduction can be applied to the whole orbit. This is only possible if condition (C4) holds, that is:
    *   (C4): If x is feasible and \f$\sigma(x) = \sigma(\gamma(x))$, then \f$\gamma(x)$ is also feasible.
    * Note that this condition can be violated. So, we additionally need to check this condition. We do this by making
    * sure that the intersection of the variable domains is non-empty.
    */
   *nchosenperms = 0;

   /* for creating the provisional disjoint set with minimum, variable lower bounds must be an array. Allocate. */
   arraysgenerated = varlbs == NULL;
   if ( arraysgenerated )
   {
      assert( varubs == NULL );
      SCIPallocBufferArray(scip, &varlbs, ofdata->npermvars);
      for (i = 0; i < ofdata->npermvars; ++i)
         varlbs[i] = SCIPvarGetLbLocal(ofdata->permvars[i]);

      SCIPallocBufferArray(scip, &varubs, ofdata->npermvars);
      for (i = 0; i < ofdata->npermvars; ++i)
         varlbs[i] = SCIPvarGetUbLocal(ofdata->permvars[i]);
   }

   /* create disjoint set that maintains orbit information of the group generated by chosenperms */
   SCIP_CALL( provisionalDisjointSetWithMinimumCreate(scip, &orbitmindjset, varlbs, varlbs, varubs,
      ofdata->npermvars, ofdata->npermvars * 2) );

   /* @todo investigate handling permutations in different order (algorithm is greedy) */
   for (p = 0; p < ofdata->nperms; ++p)
   {
      perm = ofdata->perms[p];

      /* provisionally, compute orbit information of perm added to chosenperms */
      for (i = 0; i < ofdata->npermvars; ++i)
      {
         if ( i < perm[i] )
         {
            SCIP_CALL( provisionalDisjointSetWithMinimumUnion(scip, orbitmindjset, i, perm[i]) );
         }
      }

      /* iterate over each branched variable and check */
      for (i = 0; i < nbranchedvarindices; ++i)
      {
         varid = branchedvarindices[i];
         varidimage = perm[varid];

         /* this variable branching does not affect this permutation */
         if ( varidimage == varid )
            continue;

         /* permit permutation only if the branching variable upper bound is smaller or equal to the minimum in the
          * orbit of the group generated by chosenperms added with perm.
          * Observe that the variable is part of its own orbit, so if it is a free variable, then this check always
          * disqualifies this permutation as it should.
          */
         SCIP_CALL( provisionalDisjointSetWithMinimumFind(scip, orbitmindjset, varid, NULL,
            &orbitlbmin, NULL, NULL) );
         if ( GT(scip, varubs[varid], orbitlbmin) )
         {
            /* If the upper bound of the branching variable is larger than the minimal lower bound of variables in the
             * orbit, disqualify this permutation. */
            break;
         }
      }
      /* if the above loop is broken, this permutation does not qualify for the stabilizer */
      if ( i < nbranchedvarindices )
      {
         /* undo provisional orbit information updates */
         provisionalDisjointSetWithMinimumUndo(orbitmindjset);
         continue;
      }

      /* check condition (C4) */
      for (i = 0; i < ofdata->npermvars; ++i)
      {
         /* @todo better is to iterate over every component in disjoint set, if possible */
         /* @todo instead of checking perm[i] == i everywhere, can also precompute the support of a permutation */
         if ( perm[i] == i )
            continue;

         SCIP_CALL( provisionalDisjointSetWithMinimumFind(scip, orbitmindjset, i, NULL,
            NULL, &orbitlbmax, &orbitubmin) );
         /* the intersection of the variable domains in the orbit of i is empty */
         if ( GT(scip, orbitlbmax, orbitubmin) )
            break;
      }
      /* if the above loop is broken, this permutation does not qualify for the stabilizer */
      if ( i < ofdata->npermvars )
      {
         /* undo provisional orbit information updates */
         provisionalDisjointSetWithMinimumUndo(orbitmindjset);
         continue;
      }

      /* permutation qualifies for the stabilizer. Add permutation */
      chosenperms[(*nchosenperms)++] = perm;

      /* make provisional updates to orbit information permanent */
      provisionalDisjointSetWithMinimumConfirm(orbitmindjset);
   }

   /* free djset */
   provisionalDisjointSetWithMinimumFree(scip, &orbitmindjset);

   if ( arraysgenerated )
   {
      SCIPfreeBufferArray(scip, &varubs);
      SCIPfreeBufferArray(scip, &varlbs);
      varlbs = NULL;
      varubs = NULL;
   }

   return SCIP_OKAY;
}

/** using bisection, finds the minimal entry of firstleq such that ids[idssort[*firstleq]] >= findid */
static
int bisectSortedArrayFindFirstLEQ(
   int*               ids,                /**< int array with entries */
   int*               idssort,            /**< array of indices of ids that sort ids */
   int                frameleft,          /**< search in idssort for index range [frameleft, frameright) */
   int                frameright,         /**< search in idssort for index range [frameleft, frameright) */
   int                findid              /**< entry value to find */
)
{
   int center;
   int id;

#ifndef NDEBUG
   int origframeleft;
   int origframeright;
   origframeleft = frameleft;
   origframeright = frameright;
#endif

   assert( ids != NULL );
   assert( idssort != NULL );
   assert( frameleft >= 0 );
   assert( frameright > frameleft );

   while (frameright - frameleft >= 2)
   {
      /* split [frameleft, frameright) in [frameleft, center) and [center, frameright) */
      center = frameleft + ((frameright - frameleft) / 2);
      assert( center > frameleft );
      assert( center < frameright );
      id = idssort[center];
      if ( ids[id] < findid )
      {
         /* first instance of orbitsetcomponentid is in [center, frameright) */
         frameleft = center;
      }
      else
      {
         /* first instance of orbitsetcomponentid is in [frameleft, center) */
         frameright = center;
      }
   }

   assert( frameright - frameleft == 1 );
   id = idssort[frameleft];
   if ( ids[id] < findid )
      ++frameleft;

   assert( frameleft >= origframeleft );
   assert( frameright <= origframeright );
   assert( frameleft >= origframeright || ids[idssort[frameleft]] >= findid );
   assert( frameleft - 1 < origframeleft || ids[idssort[frameleft - 1]] < findid );
   return frameleft;
}

/** dynamic orbital fixing, the orbital branching part */
static
SCIP_RETCODE orbitalFixingDynamicApplyOrbitalBranchingPropagations(
   SCIP*              scip,               /**< pointer to SCIP data structure */
   OFCDATA*           ofdata,             /**< pointer to data for orbital fixing data */
   SCIP_SHADOWTREE*   shadowtree,         /**< pointer to shadow tree */
   SCIP_Bool*         infeasible,         /**< pointer to store whether infeasibility is detected */
   int*               ngen                /**< pointer to store the number of determined variable domain reductions */
)
{
   SCIP_NODE* focusnode;
   SCIP_NODE* parentnode;
   SCIP_SHADOWNODE* shadowfocusnode;
   SCIP_SHADOWNODE* tmpshadownode;
   SCIP_SHADOWNODE** rootedshadowpath;
   int pathlength;
   int depth;
   int branchstep;
   int i;
   SCIP_Real* varlbs;
   SCIP_Real* varubs;
   SCIP_SHADOWBOUNDUPDATE* update;
   int* branchedvarindices;
   SCIP_Bool* inbranchedvarindices;
   int nbranchedvarindices;
   int varid;
   SCIP_SHADOWBOUNDUPDATE* branchingdecision;
   int branchingdecisionvarid;
   int** chosenperms;
   int* perm;
   int nchosenperms;
   int p;
   int* varorbitids;
   int* varorbitidssort;
   int orbitid;
   int orbitbegin;
   int orbitend;
   int idx;
   SCIP_Real orbitlb;
   SCIP_Real orbitub;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_DISJOINTSET* orbitset;
   int orbitsetcomponentid;

   assert( scip != NULL );
   assert( ofdata != NULL );
   assert( shadowtree != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   focusnode = SCIPgetFocusNode(scip);
   assert( focusnode == SCIPgetCurrentNode(scip) );
   assert( focusnode != NULL );

   /* do nothing if this method has already been called for this node */
   if ( ofdata->lastnode == focusnode )
      return SCIP_OKAY;

   ofdata->lastnode = focusnode;
   parentnode = SCIPnodeGetParent(focusnode);

   /* the root node has not been generated by branching decisions */
   if ( parentnode == NULL )
      return SCIP_OKAY;

   shadowfocusnode = SCIPshadowTreeGetShadowNode(shadowtree, focusnode);
   assert( shadowfocusnode != NULL );

   /* get the rooted path */
   /* @todo add depth field to shadow tree node to improve efficiency */
   pathlength = 0;
   tmpshadownode = shadowfocusnode;
   do
   {
      tmpshadownode = tmpshadownode->parent;
      ++pathlength;
   }
   while ( tmpshadownode != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &rootedshadowpath, pathlength) );
   i = pathlength;
   tmpshadownode = shadowfocusnode;
   while ( i > 0 )
   {
      rootedshadowpath[--i] = tmpshadownode;
      assert( tmpshadownode != NULL );
      tmpshadownode = tmpshadownode->parent;
   }
   assert( tmpshadownode == NULL );
   assert( i == 0 );

   /* replay fixings and propagations made until just before the focusnode */
   assert( ofdata->npermvars > 0 );  /* if it's 0, then we do not have to do anything at all */

   SCIP_CALL( SCIPallocBufferArray(scip, &varlbs, ofdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varubs, ofdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &branchedvarindices, ofdata->npermvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &inbranchedvarindices, ofdata->npermvars) );

   /* start with the bounds found after computing the symmetry group */
   for (i = 0; i < ofdata->npermvars; ++i)
      varlbs[i] = ofdata->globalvarlbs[i];
   for (i = 0; i < ofdata->npermvars; ++i)
      varubs[i] = ofdata->globalvarubs[i];

   nbranchedvarindices = 0;
   for (depth = 0; depth < pathlength - 1; ++depth)
   {
      tmpshadownode = rootedshadowpath[depth];

      /* receive propagations */
      for (i = 0; i < tmpshadownode->npropagations; ++i)
      {
         update = &(tmpshadownode->propagations[i]);
         varid = SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) update->var);
         assert( varid < ofdata->npermvars || varid == INT_MAX );
         assert( varid >= 0 );
         if ( varid < ofdata->npermvars )
         {
            assert( LE(scip, varlbs[varid], varubs[varid]) );
            switch (update->boundchgtype)
            {
               case SCIP_BOUNDTYPE_LOWER:
                  assert( GE(scip, update->newbound, varlbs[varid]) );
                  varlbs[varid] = update->newbound;
                  break;
               case SCIP_BOUNDTYPE_UPPER:
                  assert( LE(scip, update->newbound, varubs[varid]) );
                  varubs[varid] = update->newbound;
                  break;
               default:
                  assert( FALSE );
            }
            assert( LE(scip, varlbs[varid], varubs[varid]) );
         }
      }

      /* receive variable indices of branched variables */
      for (i = 0; i < tmpshadownode->nbranchingdecisions; ++i)
      {
         update = &(tmpshadownode->branchingdecisions[i]);
         varid = SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) update->var);
         assert( varid < ofdata->npermvars || varid == INT_MAX );
         assert( varid >= 0 );
         if ( varid < ofdata->npermvars )
         {
            if ( inbranchedvarindices[varid] )
               continue;
            branchedvarindices[nbranchedvarindices++] = varid;
            inbranchedvarindices[varid] = TRUE;
         }
      }
   }

   /* determine symmetry group at this point, apply branched variable, apply orbital branching for this
    * 
    * The branching variables are applied one-after-the-other.
    * So, the group before branching is determined, orbital branching to the branching variable, then the branching 
    * variable is applied, and possibly repeated for other branching variables.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &chosenperms, ofdata->nperms) );
   for (branchstep = 0; branchstep < shadowfocusnode->nbranchingdecisions; ++branchstep)
   {
      branchingdecision = &(shadowfocusnode->branchingdecisions[branchstep]);
      branchingdecisionvarid = SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) branchingdecision->var);
      assert( branchingdecisionvarid < ofdata->npermvars || branchingdecisionvarid == INT_MAX );
      assert( branchingdecisionvarid >= 0 );
      assert( branchingdecision->boundchgtype == SCIP_BOUNDTYPE_LOWER ?
         LE(scip, varlbs[branchingdecisionvarid], branchingdecision->newbound) :
         GE(scip, varubs[branchingdecisionvarid], branchingdecision->newbound) );
      assert( LE(scip, varlbs[branchingdecisionvarid], varubs[branchingdecisionvarid]) );

      /* branching decision will not have an effect on this */
      if ( branchingdecisionvarid >= ofdata->npermvars )
         continue;
      assert( branchingdecisionvarid >= 0 && branchingdecisionvarid < ofdata->npermvars );

      /* get the generating set of permutations of a subgroup of a stabilizing symmetry subgroup.
       *
       * Note: All information about branching decisions is kept in varlbs, varubs, and the branchedvarindices.
       */
      SCIP_CALL( orbitalFixingDynamicGetSymmetryStabilizerSubgroup(scip, ofdata, chosenperms, &nchosenperms,
         varlbs, varubs, branchedvarindices, inbranchedvarindices, nbranchedvarindices) );

      /* compute orbit containing branching var */
      SCIP_CALL( SCIPcreateDisjointset(scip, &orbitset, ofdata->npermvars) );

      /* put elements mapping to each other in same orbit */
      /* @todo a potential performance hazard; quadratic time */
      for (p = 0; p < nchosenperms; ++p)
      {
         perm = chosenperms[p];
         for (i = 0; i < ofdata->npermvars; ++i)
         {
            if ( i != perm[i] )
               SCIPdisjointsetUnion(orbitset, i, perm[i], FALSE);
         }
      }

      /* 1. ensure that the bounds are tightest possible just before the branching step (orbital fixing step) 
       *
       * If complete propagation was applied in the previous node,
       * then all variables in the same orbit have the same bounds just before branching,
       * so the bounds of the branching variable should be the tightest in its orbit by now.
       * It is possible that that is not the case. In that case, we do it here. 
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &varorbitids, ofdata->npermvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varorbitidssort, ofdata->npermvars) );
      for (i = 0; i < ofdata->npermvars; ++i)
         varorbitids[i] = SCIPdisjointsetFind(orbitset, i);
      SCIPsort(varorbitidssort, SCIPsortArgsortInt, varorbitids, ofdata->npermvars);

      orbitid = -1;
      orbitbegin = 0;
      for (orbitbegin = 0; orbitbegin < ofdata->npermvars; orbitbegin = orbitend)
      {
         /* get id of the orbit, and scan how large the orbit is */
         orbitid = varorbitids[varorbitidssort[orbitbegin]];
         for (orbitend = orbitbegin + 1; orbitend < ofdata->npermvars; ++orbitend)
         {
            if ( varorbitids[varorbitidssort[orbitend]] != orbitid )
               break;
         }

         /* orbits consisting of only one element cannot yield reductions */
         if ( orbitend - orbitbegin <= 1 )
            continue;

         /* get upper and lower bounds in orbit */
         orbitlb = -INFINITY;
         orbitub = INFINITY;
         for (i = orbitbegin; i < orbitend; ++i)
         {
            varid = varorbitidssort[i];
            assert( varid >= 0 );
            assert( varid < ofdata->npermvars );

            lb = varlbs[varid];
            if ( SCIPisInfinity(scip, -orbitlb) || GT(scip, lb, orbitlb) )
               orbitlb = lb;
            ub = varubs[varid];
            if ( SCIPisInfinity(scip, orbitub) || LT(scip, ub, orbitub) )
               orbitub = ub;
         }

         /* if bounds are incompatible, infeasibility is detected */
         if ( GT(scip, orbitlb, orbitub) )
         {
            *infeasible = TRUE;
            goto FREE;
         }
         assert( LE(scip, orbitlb, orbitub) );

         /* update variable bounds to be in this range */
         for (i = orbitbegin; i < orbitend; ++i)
         {
            varid = varorbitidssort[i];
            assert( varid >= 0 );
            assert( varid < ofdata->npermvars );

            assert( LE(scip, varlbs[varid], orbitlb) );
            varlbs[varid] = orbitlb;
            if ( !SCIPisInfinity(scip, -orbitlb) &&
               LT(scip, SCIPvarGetLbLocal(ofdata->permvars[varid]), orbitlb) )
            {
               SCIP_Bool tightened;
               SCIP_CALL( SCIPtightenVarLb(scip, ofdata->permvars[varid], orbitlb, TRUE, infeasible, &tightened) );

               /* propagator detected infeasibility in this node. Jump out of loop towards freeing everything. */
               if ( *infeasible )
                  goto FREE;
               assert( tightened );
               *ngen += 1;
            }

            assert( GE(scip, varubs[varid], orbitub) );
            varubs[varid] = orbitub;
            if ( !SCIPisInfinity(scip, orbitub) &&
               GT(scip, SCIPvarGetUbLocal(ofdata->permvars[varid]), orbitub) )
            {
               SCIP_Bool tightened;
               SCIP_CALL( SCIPtightenVarUb(scip, ofdata->permvars[varid], orbitub, TRUE, infeasible, &tightened) );

               /* propagator detected infeasibility in this node. Jump out of loop towards freeing everything. */
               if ( *infeasible )
                  goto FREE;
               assert( tightened );
               *ngen += 1;
            }
         }
      }
      assert( !*infeasible );

      /* 2. apply branching step to varlbs or varubs array
       *
       * Due to the steps above, it is possible that the branching step is redundant or infeasible. */
      assert( LE(scip, varlbs[branchingdecisionvarid], varubs[branchingdecisionvarid]) );
      switch (branchingdecision->boundchgtype)
      {
         case SCIP_BOUNDTYPE_LOWER:
            /* incompatible upper bound */
            if ( GT(scip, branchingdecision->newbound, varubs[branchingdecisionvarid]) )
            {
               *infeasible = TRUE;
               goto FREE;
            }

            assert( LE(scip, varlbs[branchingdecisionvarid], branchingdecision->newbound) );
            varlbs[branchingdecisionvarid] = branchingdecision->newbound;
            break;
         case SCIP_BOUNDTYPE_UPPER:
            /* incompatible lower bound */
            if ( LT(scip, branchingdecision->newbound, varlbs[branchingdecisionvarid]) )
            {
               *infeasible = TRUE;
               goto FREE;
            }

            assert( GE(scip, varubs[branchingdecisionvarid], branchingdecision->newbound) );
            varubs[branchingdecisionvarid] = branchingdecision->newbound;
            break;
         default:
            assert( FALSE );
      }

      /* 3. propagate that branching variable is >= the variables in its orbit
       *
       * Also apply the updates to the variable arrays 
       */

      /* get the orbit of the branching variable */
      orbitsetcomponentid = SCIPdisjointsetFind(orbitset, branchingdecisionvarid);

      /* find the orbit in the sorted array of orbits. npermvars can be huge, so use bisection. */
      orbitbegin = bisectSortedArrayFindFirstLEQ(varorbitids, varorbitidssort, 0, ofdata->npermvars, 
         orbitsetcomponentid);
      assert( orbitbegin >= 0 && orbitbegin < ofdata->npermvars );
      assert( varorbitids[varorbitidssort[orbitbegin]] == orbitsetcomponentid );
      assert( orbitbegin == 0 || varorbitids[varorbitidssort[orbitbegin - 1]] < orbitsetcomponentid );

      orbitend = bisectSortedArrayFindFirstLEQ(varorbitids, varorbitidssort, orbitbegin + 1, ofdata->npermvars, 
         orbitsetcomponentid + 1);
      assert( orbitend > 0 && orbitend <= ofdata->npermvars && orbitend > orbitbegin );
      assert( orbitend == ofdata->npermvars || varorbitids[varorbitidssort[orbitend]] > orbitsetcomponentid );
      assert( varorbitids[varorbitidssort[orbitend - 1]] == orbitsetcomponentid );

      /* propagate that branching variable is >= the variables in its orbit */
      for (idx = orbitbegin; idx < orbitend; ++idx)
      {
         varid = varorbitidssort[idx];
         assert( varorbitids[varid] == orbitsetcomponentid );

         /* ignore current branching variable */
         if ( varid == branchingdecisionvarid )
            continue;

         /* is variable varid in the orbit? */
         if ( SCIPdisjointsetFind(orbitset, varid) != orbitsetcomponentid )
            continue;

         /* all variables in the same orbit have the same bounds just before branching,
          * due to generalized orbital fixing. If that was not the case, these steps are applied just before applying
          * the branching step above. After the branching step, the branching variable bounds are most restricted.
          */
         assert( SCIPisInfinity(scip, -varlbs[branchingdecisionvarid])
            || GE(scip, varlbs[branchingdecisionvarid], varlbs[varid]) );
         assert( SCIPisInfinity(scip, varubs[branchingdecisionvarid])
            || LE(scip, varubs[branchingdecisionvarid], varubs[varid]) );
         /* bound changes already made could only have tightened the variable domains we are thinking about */
         assert( GE(scip, SCIPvarGetLbLocal(ofdata->permvars[varid]), varlbs[varid]) );
         assert( LE(scip, SCIPvarGetUbLocal(ofdata->permvars[varid]), varubs[varid]) );

         /* for branching variable x and variable y in its orbit, propagate x >= y. */
         /* modify UB of y-variables */
         assert( GE(scip, varubs[varid], varubs[branchingdecisionvarid]) );
         varubs[varid] = varubs[branchingdecisionvarid];
         if ( GT(scip, SCIPvarGetUbLocal(ofdata->permvars[varid]), varubs[branchingdecisionvarid]) )
         {
            SCIP_Bool tightened;
            SCIP_CALL( SCIPtightenVarUb(scip, ofdata->permvars[varid], varubs[branchingdecisionvarid], TRUE,
                  infeasible, &tightened) );

            /* propagator detected infeasibility in this node. Jump out of loop towards freeing everything. */
            if ( *infeasible )
               goto FREE;
            assert( tightened );
            *ngen += 1;
         }

         /* because variable domains are initially the same, the LB of the x-variables does not need to be modified. */
         assert( LE(scip, varlbs[varid], varlbs[branchingdecisionvarid]) );
      }

      FREE:
      SCIPfreeBufferArray(scip, &varorbitidssort);
      SCIPfreeBufferArray(scip, &varorbitids);
      SCIPfreeDisjointset(scip, &orbitset);

      if ( *infeasible )
         break;

      /* for the next branched variable at this node, if it's not already added,
       * mark the branching variable of this iteration as a branching variable. */
      if ( !inbranchedvarindices[branchingdecisionvarid] )
      {
         assert( nbranchedvarindices < ofdata->npermvars );
         branchedvarindices[nbranchedvarindices++] = branchingdecisionvarid;
         inbranchedvarindices[branchingdecisionvarid] = TRUE;
      }
   }
   SCIPfreeBufferArray(scip, &chosenperms);

   /* clean inbranchedvarindices array */
   for (i = 0; i < nbranchedvarindices; ++i)
   {
      varid = branchedvarindices[i];
      assert( varid >= 0 );
      assert( varid < ofdata->npermvars );
      assert( inbranchedvarindices[varid] );
      inbranchedvarindices[varid] = FALSE;
   }
#ifndef NDEBUG
   for (i = 0; i < ofdata->npermvars; ++i)
   {
      assert( inbranchedvarindices[i] == FALSE );
   }
#endif
   
   /* free everything */
   SCIPfreeCleanBufferArray(scip, &inbranchedvarindices);
   SCIPfreeBufferArray(scip, &branchedvarindices);
   SCIPfreeBufferArray(scip, &varubs);
   SCIPfreeBufferArray(scip, &varlbs);
   SCIPfreeBufferArray(scip, &rootedshadowpath);

   return SCIP_OKAY;
}

/** dynamic orbital fixing, the orbital fixing part */
static
SCIP_RETCODE orbitalFixingDynamicApplyOrbitalFixingPropagations(
   SCIP*              scip,               /**< pointer to SCIP data structure */
   OFCDATA*           ofdata,             /**< pointer to data for orbital fixing data */
   SCIP_SHADOWTREE*   shadowtree,         /**< pointer to shadow tree */
   SCIP_Bool*         infeasible,         /**< pointer to store whether infeasibility is detected */
   int*               ngen                /**< pointer to store the number of determined variable domain reductions */
)
{
   SCIP_NODE* focusnode;
   SCIP_SHADOWNODE* shadowfocusnode;
   SCIP_SHADOWNODE* tmpshadownode;
   int i;
   SCIP_SHADOWBOUNDUPDATE* update;
   int* branchedvarindices;
   SCIP_Bool* inbranchedvarindices;
   int nbranchedvarindices;
   int varid;
   int** chosenperms;
   int* perm;
   int nchosenperms;
   int p;
   SCIP_VAR* var;
   SCIP_DISJOINTSET* orbitset;
   int* varorbitids;
   int* varorbitidssort;
   int orbitid;
   int orbitbegin;
   int orbitend;
   SCIP_Real orbitlb;
   SCIP_Real orbitub;
   SCIP_Real lb;
   SCIP_Real ub;

   assert( scip != NULL );
   assert( ofdata != NULL );
   assert( shadowtree != NULL );
   assert( infeasible != NULL );
   assert( ngen != NULL );

   focusnode = SCIPgetFocusNode(scip);
   assert( focusnode == SCIPgetCurrentNode(scip) );
   assert( focusnode != NULL );

   shadowfocusnode = SCIPshadowTreeGetShadowNode(shadowtree, focusnode);
   assert( shadowfocusnode != NULL );

   /* get the branching variables until present, so including the branchings of the focusnode */
   assert( ofdata->npermvars > 0 );  /* if it's 0, then we do not have to do anything at all */

   SCIP_CALL( SCIPallocBufferArray(scip, &branchedvarindices, ofdata->npermvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &inbranchedvarindices, ofdata->npermvars) );

   nbranchedvarindices = 0;
   tmpshadownode = shadowfocusnode;
   while ( tmpshadownode != NULL )
   {
      /* receive variable indices of branched variables */
      for (i = 0; i < tmpshadownode->nbranchingdecisions; ++i)
      {
         update = &(tmpshadownode->branchingdecisions[i]);
         varid = SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) update->var);
         assert( varid < ofdata->npermvars || varid == INT_MAX );
         assert( varid >= 0 );
         if ( varid < ofdata->npermvars )
         {
            if ( inbranchedvarindices[varid] )
               continue;
            branchedvarindices[nbranchedvarindices++] = varid;
            inbranchedvarindices[varid] = TRUE;
         }
      }
      tmpshadownode = tmpshadownode->parent;
   }

   /* 1. compute the orbit of the branching variable of the stabilized symmetry subgroup at this point. */
   /* 1.1. identify the permutations of the symmetry group that are permitted */
   SCIP_CALL( SCIPallocBufferArray(scip, &chosenperms, ofdata->nperms) );
   SCIP_CALL( orbitalFixingDynamicGetSymmetryStabilizerSubgroup(scip, ofdata, chosenperms, &nchosenperms,
      NULL, NULL, branchedvarindices, inbranchedvarindices, nbranchedvarindices) );

   /* 1.2. compute orbits of this subgroup */
   SCIP_CALL( SCIPcreateDisjointset(scip, &orbitset, ofdata->npermvars) );

   /* put elements mapping to each other in same orbit */
   /* @todo this is O(nchosenperms * npermvars), which is a potential performance bottleneck.
      Alternative: precompute support per permutation at initialization, and iterate over these.*/
   for (p = 0; p < nchosenperms; ++p)
   {
      perm = chosenperms[p];
      for (i = 0; i < ofdata->npermvars; ++i)
      {
         if ( i != perm[i] )
            SCIPdisjointsetUnion(orbitset, i, perm[i], FALSE);
      }
   }

   /* 2. for each orbit, take the intersection of the domains */
   SCIP_CALL( SCIPallocBufferArray(scip, &varorbitids, ofdata->npermvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varorbitidssort, ofdata->npermvars) );
   for (i = 0; i < ofdata->npermvars; ++i)
      varorbitids[i] = SCIPdisjointsetFind(orbitset, i);
   SCIPsort(varorbitidssort, SCIPsortArgsortInt, varorbitids, ofdata->npermvars);

   orbitid = -1;
   orbitbegin = 0;
   for (orbitbegin = 0; orbitbegin < ofdata->npermvars; orbitbegin = orbitend)
   {
      /* get id of the orbit, and scan how large the orbit is */
      orbitid = varorbitids[varorbitidssort[orbitbegin]];
      for (orbitend = orbitbegin + 1; orbitend < ofdata->npermvars; ++orbitend)
      {
         if ( varorbitids[varorbitidssort[orbitend]] != orbitid )
            break;
      }

      /* orbits consisting of only one element cannot yield reductions */
      if ( orbitend - orbitbegin <= 1 )
         continue;

      /* get upper and lower bounds in orbit */
      orbitlb = -INFINITY;
      orbitub = INFINITY;
      for (i = orbitbegin; i < orbitend; ++i)
      {
         varid = varorbitidssort[i];
         assert( varid >= 0 );
         assert( varid < ofdata->npermvars );

         var = ofdata->permvars[varid];
         assert( var != NULL );

         lb = SCIPvarGetLbLocal(var);
         if ( GT(scip, lb, orbitlb) )
            orbitlb = lb;
         ub = SCIPvarGetUbLocal(var);
         if ( LT(scip, ub, orbitub) )
            orbitub = ub;
      }

      /* if bounds are incompatible, infeasibility is detected */
      if ( GT(scip, orbitlb, orbitub) )
      {
         *infeasible = TRUE;
         goto FREE;
      }
      assert( LE(scip, orbitlb, orbitub) );

      /* update variable bounds to be in this range */
      for (i = orbitbegin; i < orbitend; ++i)
      {
         varid = varorbitidssort[i];
         assert( varid >= 0 );
         assert( varid < ofdata->npermvars );

         var = ofdata->permvars[varid];
         assert( var != NULL );

         if ( LT(scip, SCIPvarGetLbLocal(var), orbitlb) )
         {
            SCIP_Bool tightened;
            SCIP_CALL( SCIPtightenVarLb(scip, var, orbitlb, TRUE, infeasible, &tightened) );

            /* propagator detected infeasibility in this node. Jump out of loop towards freeing everything. */
            if ( *infeasible )
               goto FREE;
            assert( tightened );
            *ngen += 1;
         }

         if ( GT(scip, SCIPvarGetUbLocal(var), orbitub) )
         {
            SCIP_Bool tightened;
            SCIP_CALL( SCIPtightenVarUb(scip, var, orbitub, TRUE, infeasible, &tightened) );

            /* propagator detected infeasibility in this node. Jump out of loop towards freeing everything. */
            if ( *infeasible )
               goto FREE;
            assert( tightened );
            *ngen += 1;
         }
      }
   }

   FREE:
   SCIPfreeBufferArray(scip, &varorbitidssort);
   SCIPfreeBufferArray(scip, &varorbitids);
   SCIPfreeDisjointset(scip, &orbitset);
   SCIPfreeBufferArray(scip, &chosenperms);

   /* clean inbranchedvarindices array */
   for (i = 0; i < nbranchedvarindices; ++i)
   {
      varid = branchedvarindices[i];
      assert( varid >= 0 );
      assert( varid < ofdata->npermvars );
      assert( inbranchedvarindices[varid] );
      inbranchedvarindices[varid] = FALSE;
   }
#ifndef NDEBUG
   for (i = 0; i < ofdata->npermvars; ++i)
   {
      assert( inbranchedvarindices[i] == FALSE );
   }
#endif

   SCIPfreeCleanBufferArray(scip, &inbranchedvarindices);
   SCIPfreeBufferArray(scip, &branchedvarindices);

   return SCIP_OKAY;
}


/** apply generalized orbital fixing on a symmetry group component using a two step mechanism
 *
 * 1. At the parent of our focus node (which is the current node, because we're not probing),
 *    compute the symmetry group just before branching. Then, for our branching variable x with variable y in its
 *    orbit, we mimic adding the constraint x >= y by variable bound propagations in this node.
 *
 *    In principle, this generalizes orbital branching in the binary case: propagation of x >= y yields
 *       1. In the 1-branch: 1 = x >= y is a tautology (since y is in {0, 1}). Nothing happens.
 *       0. In the 0-branch: 0 = x >= y implies y = 0. This is an exact description of orbital branching.
 *    REF: Ostrowski, James, et al. "Orbital branching." Mathematical Programming 126.1 (2011): 147-178.
 *
 *    (This only needs to be done once per node.)
 *
 * 2. At the focus node itself, compute the symmetry group.
 *    The symmetry group in this branch-and-bound tree node is a subgroup of the problem symmetry group
 *    as described in the function orbitalFixingDynamicGetSymmetryStabilizerSubgroup.
 *    For this symmetry subgroup, in each orbit, update the variable domains with the intersection of all variable
 *    domains in the orbit.
 *
 *    This generalizes orbital fixing in the binary case.
 *    REF: Margot 2002, Margot 2003, Orbital Branching, Ostrowski's PhD thesis.
 *
 */
static
SCIP_RETCODE orbitalFixingPropagateComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA*          orbifixdata,        /**< orbital fixing data structure */
   OFCDATA*              ofdata,             /**< orbital fixing component data */
   SCIP_SHADOWTREE*      shadowtree,         /**< pointer to shadow tree */
   SCIP_Bool*            infeasible,         /**< whether infeasibility is found */
   int*                  nred                /**< number of domain reductions */
   )
{
   /* step 1 */
   SCIP_CALL( orbitalFixingDynamicApplyOrbitalBranchingPropagations(scip, ofdata, shadowtree, infeasible, nred) );
   if ( *infeasible )
      return SCIP_OKAY;

   /* step 2 */
   SCIP_CALL( orbitalFixingDynamicApplyOrbitalFixingPropagations(scip, ofdata, shadowtree, infeasible, nred) );
   if ( *infeasible )
      return SCIP_OKAY;

   return SCIP_OKAY;
}


/** adds component */
static
SCIP_RETCODE addComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA*          orbifixdata,        /**< pointer to the dynamic orbital fixing data */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int**                 perms,              /**< permutations in the component */
   int                   nperms              /**< number of permutations in the component */
   )
{
   OFCDATA* ofdata;
   int i;
   int j;
   int p;
   int* origperm;
   int* newperm;
   int newidx;
   int newpermidx;

   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( perms != NULL );
   assert( nperms > 0 );

   SCIP_CALL( SCIPallocBlockMemory(scip, &ofdata) );

   /* correct indices by removing fixed points */

   /* determine the number of vars that are moved by the component, assign to ofdata->npermvars */
   ofdata->npermvars = 0;
   for (i = 0; i < npermvars; ++i)
   {
      /* is index i moved by any of the permutations in the component? */
      for (p = 0; p < nperms; ++p)
      {
         if ( perms[p][i] != i )
         {
            ++ofdata->npermvars;
            break;
         }
      }
   }

   /* do not support the setting where the component is empty */
   if ( ofdata->npermvars <= 0 )
      return SCIP_ERROR;

   /* create index-corrected permvars array and the inverse */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ofdata->permvars, ofdata->npermvars) );
   SCIP_CALL( SCIPhashmapCreate(&ofdata->permvarmap, SCIPblkmem(scip), ofdata->npermvars) );

   j = 0;
   for (i = 0; i < npermvars; ++i)
   {
      /* permvars array must be unique */
      assert( SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) permvars[i]) == INT_MAX );

      /* is index i moved by any of the permutations in the component? */
      for (p = 0; p < nperms; ++p)
      {
         if ( perms[p][i] != i )
         {
            ofdata->permvars[j] = permvars[i];
            SCIP_CALL( SCIPhashmapInsertInt(ofdata->permvarmap, (void*) permvars[i], j) );
            ++j;
            break;
         }
      }
   }
   assert( j == ofdata->npermvars );

   /* allocate permutations */
   ofdata->nperms = nperms;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ofdata->perms, nperms) );
   for (p = 0; p < nperms; ++p)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ofdata->perms[p], ofdata->npermvars) );
      origperm = perms[p];
      newperm = ofdata->perms[p];

      for (i = 0; i < npermvars; ++i)
      {
         newidx = SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) permvars[i]);
         if ( newidx >= ofdata->npermvars )
            continue;
         assert( newidx >= 0 );
         assert( newidx < ofdata->npermvars );
         assert( ofdata->permvars[newidx] == permvars[i] );
         newpermidx = SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) permvars[origperm[i]]);
         assert( newpermidx >= 0 );
         assert( newidx < ofdata->npermvars ); /* this is in the orbit of any permutation, so cannot be INT_MAX */
         assert( ofdata->permvars[newpermidx] == permvars[origperm[i]] );

         newperm[newidx] = newpermidx;
      }
   }

   /* global variable bounds */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ofdata->globalvarlbs, ofdata->npermvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ofdata->globalvarubs, ofdata->npermvars) );
   for (i = 0; i < ofdata->npermvars; ++i)
   {
      ofdata->globalvarlbs[i] = SCIPvarGetLbGlobal(ofdata->permvars[i]);
      ofdata->globalvarubs[i] = SCIPvarGetUbGlobal(ofdata->permvars[i]);
   }

   /* catch global variable bound change event */
   for (i = 0; i < ofdata->npermvars; ++i)
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, ofdata->permvars[i], SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
         orbifixdata->globalfixeventhdlr, (SCIP_EVENTDATA*) ofdata, NULL) );
   }

   /* resize component array if needed */
   assert( orbifixdata->ncomponents >= 0 );
   assert( (orbifixdata->ncomponents == 0) == (orbifixdata->componentdatas == NULL) );
   assert( orbifixdata->ncomponents <= orbifixdata->maxncomponents );
   if ( orbifixdata->ncomponents == orbifixdata->maxncomponents )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, orbifixdata->ncomponents + 1);
      assert( newsize >= 0 );

      if ( orbifixdata->ncomponents == 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &orbifixdata->componentdatas, newsize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &orbifixdata->componentdatas,
            orbifixdata->ncomponents, newsize) );
      }

      orbifixdata->maxncomponents = newsize;
   }

   /* add component */
   assert( orbifixdata->ncomponents < orbifixdata->maxncomponents );
   orbifixdata->componentdatas[orbifixdata->ncomponents++] = ofdata;

   return SCIP_OKAY;
}


/** frees component */
static
SCIP_RETCODE freeComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA*          orbifixdata,        /**< pointer to the dynamic orbital fixing data */
   OFCDATA**             ofdata              /**< pointer to component data */
   )
{
   int i;
   int p;

   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( ofdata != NULL );
   assert( *ofdata != NULL );
   assert( (*ofdata)->globalvarlbs != NULL );
   assert( (*ofdata)->globalvarubs != NULL );
   assert( (*ofdata)->nperms > 0 );
   assert( (*ofdata)->npermvars > 0 );
   assert( (*ofdata)->perms != NULL );
   assert( (*ofdata)->permvarmap != NULL );
   assert( (*ofdata)->permvars != NULL );
   assert( (*ofdata)->npermvars > 0 );

   assert( SCIPisTransformed(scip) );

   /* free global variable bound change event */
   if ( SCIPgetStage(scip) != SCIP_STAGE_FREE )
   {
      /* events at the freeing stage may not be dropped, because they are already getting dropped */
      for (i = (*ofdata)->npermvars - 1; i >= 0; --i)
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*ofdata)->permvars[i],
            SCIP_EVENTTYPE_GLBCHANGED | SCIP_EVENTTYPE_GUBCHANGED,
            orbifixdata->globalfixeventhdlr, (SCIP_EVENTDATA*) (*ofdata), -1) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &(*ofdata)->globalvarubs, (*ofdata)->npermvars);
   SCIPfreeBlockMemoryArray(scip, &(*ofdata)->globalvarlbs, (*ofdata)->npermvars);

   for (p = (*ofdata)->nperms -1; p >= 0; --p)
   {
      SCIPfreeBlockMemoryArray(scip, &(*ofdata)->perms[p], (*ofdata)->npermvars);
   }
   SCIPfreeBlockMemoryArray(scip, &(*ofdata)->perms, (*ofdata)->nperms);

   SCIPhashmapFree(&(*ofdata)->permvarmap);
   SCIPfreeBlockMemoryArray(scip, &(*ofdata)->permvars, (*ofdata)->npermvars);

   SCIPfreeBlockMemory(scip, ofdata);

   return SCIP_OKAY;
}


/*
 * Event handler callback methods
 */

/** maintain global variable fixings for reductions found during presolving or at the root node */
static
SCIP_DECL_EVENTEXEC(eventExecGlobalBoundChange)
{
   OFCDATA* ofdata;
   SCIP_VAR* var;
   int varidx;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_SYMMETRY_NAME) == 0 );
   assert( event != NULL );

   /* only update the global bounds during presolving or at the root node */
   if ( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPgetNNodes(scip) > 1 )
      return SCIP_OKAY;

   ofdata = (OFCDATA*) eventdata;
   var = SCIPeventGetVar(event);
   assert( var != NULL );
   assert( SCIPvarIsTransformed(var) );

   assert( ofdata->permvarmap != NULL );
   varidx = SCIPhashmapGetImageInt(ofdata->permvarmap, (void*) var);

   switch ( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_GUBCHANGED:
      assert( ofdata->globalvarubs[varidx] == SCIPeventGetOldbound(event) );
      ofdata->globalvarubs[varidx] = SCIPeventGetNewbound(event);
      break;
   case SCIP_EVENTTYPE_GLBCHANGED:
      assert( ofdata->globalvarlbs[varidx] == SCIPeventGetOldbound(event) );
      ofdata->globalvarlbs[varidx] = SCIPeventGetNewbound(event);
      break;
   default:
      SCIPABORT();
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/*
 * Interface methods
 */

/** propagate orbital fixing */
SCIP_RETCODE SCIPorbitalFixingPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA*          orbifixdata,        /**< orbital fixing data structure */
   SCIP_Bool*            infeasible,         /**< whether infeasibility is found */
   int*                  nred,               /**< number of domain reductions */
   SCIP_Bool*            didrun              /**< whether any propagator actually ran */
   )
{
   OFCDATA* ofdata;
   SCIP_SHADOWTREE* shadowtree;
   int c;

   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( infeasible != NULL );
   assert( nred != NULL );
   assert( didrun != NULL );

   *infeasible = FALSE;
   *nred = 0;

   /* no components, no orbital fixing */
   assert( orbifixdata->ncomponents >= 0 );
   if ( orbifixdata->ncomponents == 0 )
      return SCIP_OKAY;

   /* do nothing if we are in a probing node */
   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* do not run again in repropagation, since the path to the root might have changed */
   if ( SCIPinRepropagation(scip) )
      return SCIP_OKAY;

   assert( orbifixdata->shadowtreeeventhdlr != NULL );
   shadowtree = SCIPgetShadowTree(orbifixdata->shadowtreeeventhdlr);
   assert( shadowtree != NULL );

   for (c = 0; c < orbifixdata->ncomponents; ++c)
   {
      ofdata = orbifixdata->componentdatas[c];
      assert( ofdata != NULL );
      assert( ofdata->nperms > 0 );
      SCIP_CALL( orbitalFixingPropagateComponent(scip, orbifixdata, ofdata, shadowtree, infeasible, nred) );
      *didrun = TRUE;

      if ( *infeasible )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** adds component for orbital fixing */
SCIP_RETCODE SCIPorbitalFixingAddComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA*          orbifixdata,        /**< orbital fixing data structure */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int**                 perms,              /**< permutations in the component */
   int                   nperms              /**< number of permutations in the component */
   )
{
   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( permvars != NULL );
   assert( npermvars > 0 );
   assert( perms != NULL );
   assert( nperms > 0 );

   /* dynamic symmetry reductions cannot be performed on original problem */
   assert( SCIPisTransformed(scip) );

   SCIP_CALL( addComponent(scip, orbifixdata, permvars, npermvars, perms, nperms) );

   return SCIP_OKAY;
}


/** resets orbital fixing data structure (clears all components) */
SCIP_RETCODE SCIPorbitalFixingReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA*          orbifixdata         /**< orbital fixing data structure */
   )
{
   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( orbifixdata->ncomponents >= 0 );
   assert( (orbifixdata->ncomponents == 0) == (orbifixdata->componentdatas == NULL) );
   assert( orbifixdata->ncomponents <= orbifixdata->maxncomponents );
   assert( orbifixdata->shadowtreeeventhdlr != NULL );

   while ( orbifixdata->ncomponents > 0 )
   {
      SCIP_CALL( freeComponent(scip, orbifixdata, &(orbifixdata->componentdatas[--orbifixdata->ncomponents])) );
   }

   assert( orbifixdata->ncomponents == 0 );
   SCIPfreeBlockMemoryArray(scip, &orbifixdata->componentdatas, orbifixdata->maxncomponents);
   orbifixdata->componentdatas = NULL;
   orbifixdata->maxncomponents = 0;

   return SCIP_OKAY;
}


/** free orbital fixing data */
SCIP_RETCODE SCIPorbitalFixingFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA**         orbifixdata         /**< orbital fixing data structure */
   )
{
   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( *orbifixdata != NULL );

   SCIP_CALL( SCIPorbitalFixingReset(scip, *orbifixdata) );

   SCIPfreeBlockMemory(scip, orbifixdata);
   return SCIP_OKAY;
}


/** initializes structures needed for orbital fixing
 * 
 * This is only done exactly once.
 */
SCIP_RETCODE SCIPorbitalFixingInclude(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OFDATA**         orbifixdata,        /**< pointer to orbital fixing data structure to populate */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr /**< pointer to the shadow tree eventhdlr */
   )
{
   assert( scip != NULL );
   assert( orbifixdata != NULL );
   assert( shadowtreeeventhdlr != NULL );

   SCIP_CALL( SCIPcheckStage(scip, "SCIPorbitalFixingInclude", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPallocBlockMemory(scip, orbifixdata) );

   (*orbifixdata)->componentdatas = NULL;
   (*orbifixdata)->ncomponents = 0;
   (*orbifixdata)->maxncomponents = 0;
   (*orbifixdata)->shadowtreeeventhdlr = shadowtreeeventhdlr;

   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &(*orbifixdata)->globalfixeventhdlr,
      EVENTHDLR_SYMMETRY_NAME, EVENTHDLR_SYMMETRY_DESC, eventExecGlobalBoundChange,
      (SCIP_EVENTHDLRDATA*) (*orbifixdata)) );

   return SCIP_OKAY;
}
