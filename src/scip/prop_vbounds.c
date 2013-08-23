/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_vbounds.c
 * @brief  variable upper and lower bound propagator
 * @author Stefan Heinz
 * @author Jens Schulz
 * @author Gerald Gamrath
 *
 * This propagator uses global bound information provided by SCIP to deduce global and local bound changes.
 * It can take into account
 * - implications (bound change following from specific value of a binary variable)
 * - cliques (set of binary variables, each with a corresponding value, of which at most one variable can get the value)
 * - variable lower/upper bounds (bounds of arbitrary variables that depend linearly on the value of another variable)
 *
 * The propagator does not look at a variable in whole, but at one point in time only handles one specific bound (lower
 * or upper) of a variable and deduces changes for lower or upper bounds of other variables. The concept is as follows:
 *
 * 1) Extract variable bound data
 *
 *    Implications and cliques are stored in a way such that given a variable and its new value, we can access all bound
 *    changes that can be deduced from setting the variable to that value. However, for variable bounds, this currently
 *    does not hold, they are only stored in the other direction, i.e. for a bound of a given variable, we have a list
 *    of all other bounds of variables that directly influence the bound of the given variable and a linear function
 *    describing how they do this.
 *    For the propagation, we need the other direction, thus we store it in the propagator data when the branch-and-bound
 *    solving process is about to begin.
 *
 * 2) Topological sorting of bounds of variable
 *
 *    We compute a topological order of the bounds of variables. This is needed to define an order in which we will
 *    regard bounds of variables in the propagation process in order to avoid unneccessarily regarding the same variable
 *    bound multiple times because it was changed in the meantime when propagating another bound of a variable.
 *    Therefore, we implictly regard a directed graph, in which each node corresponds to a bound of a variable and there
 *    exists a directed edge from one node to another, if the bound corresponding to the former node influences the
 *    bound corresponding to the latter node. This is done by iteratively running a DFS until all nodes were visited.
 *    Note that there might be cycles in the graph, which are randomly broken, so the order is only almost topological.
 *
 * 3) Collecting bound changes
 *
 *    For each bound of a variable, which can trigger bound changes of other variables, the propagator catches all
 *    events informing about a global change of the bound or a local toghtening of the bound. The event handler
 *    then adds the bound of the variable to a priority queue, with the key in the priority queue corresponding
 *    to the position of the bound in the topological sort.
 *
 * 4) Propagating Bounds
 *
 *    As long as there are bounds contained in the priority queue, the propagator pops one bound from the queue, which
 *    is the one most at the beginning of the topological sort, so it should not be influenced by propagating other
 *    bounds currently contained in the queue. Starting at this bound, all implication, clique, and variable bound
 *    information is used to deduce tigther bounds for other variables and change the bounds, if a tighter one is found.
 *    These bound changes trigger an event that will lead to adding the corresponding bound to the priority queue,
 *    if it is not contained, yet. The process is iterated until the priority queue contains no more bounds.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_vbounds.h"

/**@name Propagator properties
 *
 * @{
 */

#define PROP_NAME              "vbounds"
#define PROP_DESC              "propagates variable upper and lower bounds"
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP | SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY           3000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

/**@} */

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "vbounds"
#define EVENTHDLR_DESC         "bound change event handler for for vbounds propagator"

/**@} */

/**@name Default parameter values
 *
 * @{
 */

#define DEFAULT_USEBDWIDENING      TRUE      /**< should bound widening be used to initialize conflict analysis? */
#define DEFAULT_USEIMPLICS         FALSE     /**< should implications be propagated? */
#define DEFAULT_USECLIQUES         FALSE     /**< should cliques be propagated? */
#define DEFAULT_USEVBOUNDS         TRUE      /**< should variable bounds be propagated? */
#define DEFAULT_DOTOPOSORT         TRUE      /**< should the bounds be topologically sorted in advance? */
#define DEFAULT_SORTCLIQUES        FALSE     /**< should cliques be regarded for the topological sort? */

/**@} */

/**@name Propagator defines
 *
 * @{
 *
 * The propagator works on indices representing a bound of a variable. This index will be called bound index in the
 * following. For a given active variable with problem index i (note that active variables have problem indices
 * between 0 and nactivevariable - 1), the bound index of its lower bound is 2*i, the bound index of its upper
 * bound is 2*i + 1. The other way around, a given bound index i corresponds to the variable with problem index
 * i/2 (rounded down), and to the lower bound, if i is even, to the upper bound if i is odd.
 * The following macros can be used to convert bound index into variable problem index and boundtype and vice versa.
 */
#define getLbIndex(idx) (2*(idx))
#define getUbIndex(idx) (2*(idx)+1)
#define getVarIndex(idx) ((idx)/2)
#define getBoundtype(idx) (((idx) % 2 == 0) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER)
#define isIndexLowerbound(idx) ((idx) % 2 == 0)
#define getBoundString(lower) ((lower) ? "lb" : "ub")
#define getBoundtypeString(type) ((type) == SCIP_BOUNDTYPE_LOWER ? "lower" : "upper")
#define indexGetBoundString(idx) (getBoundString(isIndexLowerbound(idx)))

/**@} */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for catching bound changes */
   SCIP_VAR**            vars;               /**< array containing all variable which are considered within the propagator */
   SCIP_HASHMAP*         varhashmap;         /**< hashmap mapping from variable to index in the vars array */
   int*                  topoorder;          /**< array mapping on the bounds of variables in topological order;
                                              *   or -1, if the bound that should be at that position has no outgoing
                                              *   implications, cliques, or vbounds;
                                              *   i.e., for i < j and topoorder[i] != -1 != topoorder[j], the variable
                                              *   and boundtype represented by index topoorder[i] are earlier in the
                                              *   topological order than those represented by index topoorder[j]
                                              */
   int**                 vboundboundedidx;   /**< array storing for each bound index the bound indices of all bounds
                                              *   influenced by this bound through variable bounds */
   SCIP_Real**           vboundcoefs;        /**< array storing for each bound index the coefficients in the variable
                                              *   bounds influencing the corresponding bound index stored in
                                              *   vboundboundedidx */
   SCIP_Real**           vboundconstants;    /**< array storing for each bound index the constants in the variable
                                              *   bounds influencing the corresponding bound index stored in
                                              *   vboundboundedidx */
   int*                  nvbounds;           /**< array storing for each bound index the number of vbounds stored */
   int*                  vboundsize;         /**< array with sizes of vbound arrays for the nodes */
   int                   nbounds;            /**< number of bounds of variables regarded (two times number of active variables) */
   SCIP_PQUEUE*          propqueue;          /**< priority queue to handle the bounds of variables that were changed and have to be propagated */
   SCIP_Bool*            inqueue;            /**< boolean array to store whether a bound of a variable is already contained in propqueue */
   SCIP_Bool             initialized;        /**< was the data for propagation already initialized? */
   SCIP_Bool             usebdwidening;      /**< should bound widening be used to initialize conflict analysis? */
   SCIP_Bool             useimplics;         /**< should implications be propagated? */
   SCIP_Bool             usecliques;         /**< should cliques be propagated? */
   SCIP_Bool             usevbounds;         /**< should variable bounds be propagated? */
   SCIP_Bool             dotoposort;         /**< should the bounds be topologically sorted in advance? */
   SCIP_Bool             sortcliques;        /**< should cliques be regarded for the topological sort? */
};

/** inference information */
struct InferInfo
{
   union
   {
      struct
      {
         unsigned int    pos:31;             /**< position of the variable which forced that propagation */
         unsigned int    boundtype:1;        /**< bound type which was the reason (0: lower, 1: upper) */
      } asbits;
      int                asint;              /**< inference information as a single int value */
   } val;
};
typedef struct InferInfo INFERINFO;

/** converts an integer into an inference information */
static
INFERINFO intToInferInfo(
   int                   i                   /**< integer to convert */
   )
{
   INFERINFO inferinfo;

   inferinfo.val.asint = i;

   return inferinfo;
}

/** converts an inference information into an int */
static
int inferInfoToInt(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return inferinfo.val.asint;
}

/** returns the propagation rule stored in the inference information */
static
SCIP_BOUNDTYPE inferInfoGetBoundtype(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   assert((SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype == SCIP_BOUNDTYPE_LOWER
      || (SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype == SCIP_BOUNDTYPE_UPPER);
   return (SCIP_BOUNDTYPE)inferinfo.val.asbits.boundtype;
}

/** returns the position stored in the inference information */
static
int inferInfoGetPos(
   INFERINFO             inferinfo           /**< inference information to convert */
   )
{
   return (int) inferinfo.val.asbits.pos;
}

/** constructs an inference information out of a position of a variable and a boundtype */
static
INFERINFO getInferInfo(
   int                   pos,                /**< position of the variable which forced that propagation */
   SCIP_BOUNDTYPE        boundtype           /**< propagation rule that deduced the value */
   )
{
   INFERINFO inferinfo;

   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);
   assert((int)boundtype >= 0 && (int)boundtype <= 1); /*lint !e685 !e568q*/
   assert(pos >= 0);

   inferinfo.val.asbits.pos = (unsigned int) pos; /*lint !e732*/
   inferinfo.val.asbits.boundtype = (unsigned int) boundtype; /*lint !e641*/

   return inferinfo;
}

/*
 * Local methods
 */

/* returns the lower bound index of a variable */
static
int varGetLbIndex(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var                 /**< variable to get the index for */
   )
{
   assert(SCIPhashmapExists(propdata->varhashmap, var) == ((size_t)SCIPhashmapGetImage(propdata->varhashmap, var) > 0));

   return getLbIndex((int)(size_t)SCIPhashmapGetImage(propdata->varhashmap, var) - 1);
}

/* returns the upper bound index of a variable */
static
int varGetUbIndex(
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var                 /**< variable to get the index for */
   )
{
   assert(SCIPhashmapExists(propdata->varhashmap, var) == ((size_t)SCIPhashmapGetImage(propdata->varhashmap, var) > 0));

   return getUbIndex((int)(size_t)SCIPhashmapGetImage(propdata->varhashmap, var) - 1);
}

/** reset propagation data */
static
void resetPropdata(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   propdata->vars = NULL;
   propdata->varhashmap = NULL;
   propdata->topoorder = NULL;
   propdata->vboundboundedidx = NULL;
   propdata->vboundcoefs = NULL;
   propdata->vboundconstants = NULL;
   propdata->nvbounds = NULL;
   propdata->vboundsize = NULL;
   propdata->nbounds = 0;
   propdata->initialized = FALSE;
}

/** catches events for variables */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Bool lower;
   int nbounds;
   int v;
   int idx;

   assert(scip != NULL);
   assert(propdata != NULL);
   assert(propdata->vars != NULL);
   assert(propdata->topoorder != NULL);

   /* catch variable events according to computed eventtypes */
   eventhdlr = propdata->eventhdlr;
   assert(eventhdlr != NULL);

   vars = propdata->vars;
   nbounds = propdata->nbounds;

   /* setup events */
   for( v = 0; v < nbounds; ++v )
   {
      idx = propdata->topoorder[v];
      assert(idx >= 0 && idx < nbounds);

      var = vars[getVarIndex(idx)];
      lower = isIndexLowerbound(idx);

      /* if the bound does not influence another bound by implications, cliques, or vbounds,
       * we do not create an event and do not catch changes of the bound;
       * we mark this by setting the value in topoorder to -1
       */
      if( propdata->nvbounds[idx] == 0 && SCIPvarGetNImpls(var, lower) == 0 && SCIPvarGetNCliques(var, lower) == 0 )
      {
         propdata->topoorder[v] = -1;
         continue;
      }

      /* determine eventtype that we want to catch depending on boundtype of variable */
      if( lower )
         eventtype = SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_GLBCHANGED;
      else
         eventtype = SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_GUBCHANGED;

      SCIP_CALL( SCIPcatchVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*) (size_t) v, NULL) );
   }

   return SCIP_OKAY;
}


/** drops events for variables */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTTYPE eventtype;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Bool lower;
   int nbounds;
   int v;
   int idx;

   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;
   assert(eventhdlr != NULL);

   vars = propdata->vars;
   nbounds = propdata->nbounds;

   for( v = 0; v < nbounds; ++v )
   {
      idx = propdata->topoorder[v];

      if( idx == -1 )
         continue;

      assert(idx >= 0 && idx < nbounds);

      var = vars[getVarIndex(idx)];
      lower = isIndexLowerbound(idx);

      /* determine eventtype that we catch and now want to drop depending on boundtype of variable */
      if( lower )
         eventtype = SCIP_EVENTTYPE_LBTIGHTENED | SCIP_EVENTTYPE_GLBCHANGED;
      else
         eventtype = SCIP_EVENTTYPE_UBTIGHTENED | SCIP_EVENTTYPE_GUBCHANGED;

      SCIP_CALL( SCIPdropVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*) (size_t) v, -1) );
   }

   return SCIP_OKAY;
}

#define INITMEMSIZE 5

/* adds a vbound to the propagator data to store it internally and allow forward propagation */
static
SCIP_RETCODE addVbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   startidx,           /**< index of bound of variable influencing the other variable */
   int                   endidx,             /**< index of bound of variable which is influenced */
   SCIP_Real             coef,               /**< coefficient in the variable bound */
   SCIP_Real             constant            /**< constant in the variable bound */
   )
{
   int nvbounds;

   assert(scip != NULL);
   assert(propdata != NULL);

   if( propdata->vboundsize[startidx] == 0 )
   {
      /* allocate memory for storing vbounds */
      propdata->vboundsize[startidx] = INITMEMSIZE;

      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundboundedidx[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundcoefs[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundconstants[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
   }
   else if( propdata->nvbounds[startidx] >= propdata->vboundsize[startidx] )
   {
      /* reallocate memory for storing vbounds */
      propdata->vboundsize[startidx] = SCIPcalcMemGrowSize(scip, propdata->nvbounds[startidx] + 1);
      assert(propdata->nvbounds[startidx] < propdata->vboundsize[startidx]);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundboundedidx[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundcoefs[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
      SCIP_CALL( SCIPreallocMemoryArray(scip, &propdata->vboundconstants[startidx], propdata->vboundsize[startidx]) ); /*lint !e866*/
   }

   nvbounds = propdata->nvbounds[startidx];
   propdata->vboundboundedidx[startidx][nvbounds] = endidx;
   propdata->vboundcoefs[startidx][nvbounds] = coef;
   propdata->vboundconstants[startidx][nvbounds] = constant;
   (propdata->nvbounds[startidx])++;

   return SCIP_OKAY;
}


/** comparison method for two indices in the topoorder array, preferring higher indices because the order is reverse
 *  topological
 */
static
SCIP_DECL_SORTPTRCOMP(compVarboundIndices)
{
   int idx1 = (int)(size_t)elem1;
   int idx2 = (int)(size_t)elem2;

   return idx2 - idx1;
}

/** performs depth-first-search in the implicitly given directed graph from the given start index */
static
SCIP_RETCODE dfs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   int                   startnode,          /**< node to start the depth-first-search */
   SCIP_Bool*            visited,            /**< array to store for each node, whether it was already visited */
   int*                  dfsstack,           /**< array of size number of nodes to store the stack;
                                              *   only needed for performance reasons */
   int*                  stacknextedge,      /**< array of size number of nodes to store the number of adjacent nodes
                                              *   already visited for each node on the stack; only needed for
                                              *   performance reasons */
   int*                  dfsnodes,           /**< array of nodes that can be reached starting at startnode, in reverse
                                              *   dfs order */
   int*                  ndfsnodes           /**< pointer to store number of nodes that can be reached starting at
                                              *   startnode */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_Bool lower;
   int stacksize;
   int curridx;
   int nimpls;
   int idx;

   assert(startnode >= 0);
   assert(startnode < propdata->nbounds);
   assert(visited != NULL);
   assert(visited[startnode] == FALSE);
   assert(dfsstack != NULL);
   assert(dfsnodes != NULL);
   assert(ndfsnodes != NULL);

   vars = propdata->vars;

   /* put start node on the stack */
   dfsstack[0] = startnode;
   stacknextedge[0] = 0;
   stacksize = 1;
   idx = -1;

   /* we run until no more bounds indices are on the stack, i.e. all changed bounds were propagated */
   while( stacksize > 0 )
   {
      /* get next node from stack */
      curridx = dfsstack[stacksize - 1];

      /* mark current node as visited */
      assert(visited[curridx] == (stacknextedge[stacksize - 1] != 0));
      visited[curridx] = TRUE;

      startvar = vars[getVarIndex(curridx)];
      lower = isIndexLowerbound(curridx);

      nimpls = 0;

      if( propdata->sortcliques && propdata->usecliques && stacknextedge[stacksize - 1] == 0 )
         stacknextedge[stacksize - 1] = -1;

      /* stacknextedge is negative, if the last visited edge from the current node belongs to a clique;
       * the index of the clique in the variable's clique list equals abs(stacknextedge) - 1
       */
      if( propdata->sortcliques && propdata->usecliques && stacknextedge[stacksize - 1] < 0 )
      {
         SCIP_CLIQUE** cliques;
         int ncliques;
         int j;
         int i;
         SCIP_Bool found;

         ncliques = SCIPvarGetNCliques(startvar, lower);
         cliques = SCIPvarGetCliques(startvar, lower);
         found = FALSE;

         assert(stacknextedge[stacksize - 1] == -1 || -stacknextedge[stacksize - 1] - 1 < ncliques);

         /* iterate over all not yet handled cliques and search for an unvisited node */
         for( j = -stacknextedge[stacksize - 1] - 1; j < ncliques; ++j )
         {
            SCIP_VAR** cliquevars;
            SCIP_Bool* cliquevals;
            int ncliquevars;

            cliquevars = SCIPcliqueGetVars(cliques[j]);
            cliquevals = SCIPcliqueGetValues(cliques[j]);
            ncliquevars = SCIPcliqueGetNVars(cliques[j]);

            for( i = 0; i < ncliquevars; ++i )
            {
               if( cliquevars[i] == startvar )
                  continue;

               if( cliquevals[i] )
                  idx = varGetUbIndex(propdata, cliquevars[i]);
               else
                  idx = varGetLbIndex(propdata, cliquevars[i]);

               /* break when the first unvisited node is reached */
               if( idx >= 0 && !visited[idx] )
               {
                  found = TRUE;
                  break;
               }
            }
            if( found )
               break;
         }

         /* we stopped because we found an unhandled node and not because we reached the end of the list */
         if( found )
         {
            assert(idx >= 0);
            assert(!visited[idx]);
            assert(j < ncliques);

            SCIPdebugMessage("clique: %s(%s) -> %s(%s)\n", getBoundString(lower), SCIPvarGetName(startvar),
               indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));

            /* put the adjacent node onto the stack */
            dfsstack[stacksize] = idx;
            stacknextedge[stacksize] = 0;
            stacknextedge[stacksize - 1] = -j - 1;
            stacksize++;
            assert(stacksize <= propdata->nbounds);

            /* restart while loop, get next index from stack */
            continue;
         }
         else
         {
            /* we did not find an edge to an unhandled node given by a clique */
            stacknextedge[stacksize - 1] = 0;
         }
      }
      assert(stacknextedge[stacksize - 1] >= 0);

      /* go over edges given by implications */
      if( propdata->useimplics )
      {
         nimpls = SCIPvarGetNImpls(startvar, lower);

         if( stacknextedge[stacksize - 1] < nimpls )
         {
            SCIP_VAR** implvars;
            SCIP_BOUNDTYPE* impltypes;
            int* implids;
            int i;

            implvars = SCIPvarGetImplVars(startvar, lower);
            impltypes = SCIPvarGetImplTypes(startvar, lower);
            implids = SCIPvarGetImplIds(startvar, lower);

            for( i = stacknextedge[stacksize - 1]; i < nimpls; ++i )
            {
               /* it might happen that implications point to inactive variables (normally, those are removed when a
                * variable becomes inactive, but in some cases, it cannot be done), we have to ignore these variables
                */
               if( !SCIPvarIsActive(implvars[i]) )
                  continue;

               /* implication is just a shortcut, so we dont regard it now, because will later go the long way, anyway;
                * however, if we do regard cliques for the topological order, we use them to get a better order
                */
               if( propdata->usecliques && !propdata->sortcliques && implids[i] < 0 )
                  continue;

               idx = (impltypes[i] == SCIP_BOUNDTYPE_LOWER ? varGetLbIndex(propdata, implvars[i]) : varGetUbIndex(propdata, implvars[i]));

               /* break when the first unvisited node is reached */
               if( idx >= 0 && !visited[idx] )
                  break;
            }

            /* we stopped because we found an unhandled node and not because we reached the end of the list */
            if( i < nimpls )
            {
               assert(!visited[idx]);

               SCIPdebugMessage("impl: %s(%s) -> %s(%s)\n", getBoundString(lower), SCIPvarGetName(startvar),
                  indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));


               /* put the adjacent node onto the stack */
               dfsstack[stacksize] = idx;
               stacknextedge[stacksize] = 0;
               stacknextedge[stacksize - 1] = i + 1;
               stacksize++;
               assert(stacksize <= propdata->nbounds);

               /* restart while loop, get next index from stack */
               continue;
            }
            else
            {
               stacknextedge[stacksize - 1] = nimpls;
            }
         }
      }
      assert(stacknextedge[stacksize - 1] >= nimpls);

      /* go over edges corresponding by varbounds */
      if( propdata->usevbounds )
      {
         int nvbounds;
         int* vboundidx;
         int i;

         nvbounds = propdata->nvbounds[curridx];
         vboundidx = propdata->vboundboundedidx[curridx];

         /* iterate over all vbounds for the given bound */
         for( i = 0; i < nvbounds; ++i )
         {
            idx = vboundidx[i];
            assert(idx >= 0);

            /* break when the first unvisited node is reached */
            if( !visited[idx] )
               break;
         }

         /* we stopped because we found an unhandled node and not because we reached the end of the list */
         if( i < nvbounds )
         {
            assert(!visited[idx]);

            SCIPdebugMessage("vbound: %s(%s) -> %s(%s)\n", getBoundString(lower), SCIPvarGetName(startvar),
               indexGetBoundString(idx), SCIPvarGetName(vars[getVarIndex(idx)]));

            /* put the adjacent node onto the stack */
            dfsstack[stacksize] = idx;
            stacknextedge[stacksize] = 0;
            stacknextedge[stacksize - 1] = nimpls + i + 1;
            stacksize++;
            assert(stacksize <= propdata->nbounds);

            /* restart while loop, get next index from stack */
            continue;
         }

      }

      /* the current node was completely handled, remove it from stack */
      stacksize--;

      SCIPdebugMessage("topoorder[%d] = %s(%s)\n", *ndfsnodes, getBoundString(lower), SCIPvarGetName(startvar));

      /* store node in the sorted nodes array */
      dfsnodes[(*ndfsnodes)] = curridx;
      (*ndfsnodes)++;
   }

   return SCIP_OKAY;
}


/** sort the bounds of variables topologically */
static
SCIP_RETCODE topologicalSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   int* dfsstack;
   int* stacknextedge;
   int nsortednodes;
   int nbounds;
   int i;

   assert(scip != NULL);
   assert(propdata != NULL);

   nbounds = propdata->nbounds;

   SCIP_CALL( SCIPallocBufferArray(scip, &dfsstack, nbounds) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stacknextedge, nbounds) );

   nsortednodes = 0;

#ifndef NDEBUG
   for( i = 0; i < nbounds; ++i )
      assert(!propdata->inqueue[i]);
#endif

   /* while there are unvisited nodes, run dfs starting from one of these nodes; the dfs orders are stored in the
    * topoorder array, later dfs calls are just appended after the stacks of previous dfs calls, which gives us a
    * reverse topological order
    */
   for( i = 0; i < nbounds; ++i )
   {
      if( !propdata->inqueue[i] )
      {
         SCIP_CALL( dfs(scip, propdata, i, propdata->inqueue, dfsstack, stacknextedge, propdata->topoorder, &nsortednodes) );
      }
   }
   assert(nsortednodes == nbounds);

   BMSclearMemoryArray(propdata->inqueue, nbounds);

   SCIPfreeBufferArray(scip, &stacknextedge);
   SCIPfreeBufferArray(scip, &dfsstack);

   return SCIP_OKAY;
}

/** initializes the internal data for the variable bounds propagator */
static
SCIP_RETCODE initData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop                /**< vbounds propagator */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   int nvars;
   int nbounds;
   int startidx;
   int v;
   int n;

   assert(scip != NULL);
   assert(prop != NULL);

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(!propdata->initialized);

   SCIPdebugMessage("initialize vbounds propagator for problem <%s>\n", SCIPgetProbName(scip));

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nbounds = 2 * nvars;

   /* store size of the bounds of variables array */
   propdata->nbounds = nbounds;

   if( nbounds == 0 )
      return SCIP_OKAY;

   propdata->initialized = TRUE;

   /* prepare priority queue structure */
   SCIP_CALL( SCIPpqueueCreate(&propdata->propqueue, nvars, 2.0, compVarboundIndices) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->inqueue, nbounds) );
   BMSclearMemoryArray(propdata->inqueue, nbounds);

   /* we need to copy the variable since this array is the basis of the propagator and the corresponding variable array
    * within SCIP might change during the search
    */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &propdata->vars, vars, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&propdata->varhashmap, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * nvars)) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPhashmapInsert(propdata->varhashmap, propdata->vars[v], (void*)(size_t)(v + 1)) );
   }

   /* allocate memory for the arrays of the propdata */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->topoorder, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundboundedidx, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundcoefs, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundconstants, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->nvbounds, nbounds) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &propdata->vboundsize, nbounds) );
   BMSclearMemoryArray(propdata->vboundboundedidx, nbounds);
   BMSclearMemoryArray(propdata->vboundcoefs, nbounds);
   BMSclearMemoryArray(propdata->vboundconstants, nbounds);
   BMSclearMemoryArray(propdata->nvbounds, nbounds);
   BMSclearMemoryArray(propdata->vboundsize, nbounds);

   for( v = 0; v < nbounds; ++v )
   {
      propdata->topoorder[v] = v;
      propdata->vboundboundedidx[v] = NULL;
      propdata->vboundcoefs[v] = NULL;
      propdata->vboundconstants[v] = NULL;
      propdata->nvbounds[v] = 0;
      propdata->vboundsize[v] = 0;
   }

   /* collect information about varbounds */
   for( v = 0; v < nbounds; ++v )
   {
      SCIP_VAR** vbvars;
      SCIP_VAR* var;
      SCIP_Real* coefs;
      SCIP_Real* constants;
      SCIP_Bool lower;
      int nvbvars;

      var = vars[getVarIndex(v)];
      lower = isIndexLowerbound(v);

      /* get the variable bound informations for the current variable */
      if( lower )
      {
         vbvars = SCIPvarGetVlbVars(var);
         coefs = SCIPvarGetVlbCoefs(var);
         constants = SCIPvarGetVlbConstants(var);
         nvbvars = SCIPvarGetNVlbs(var);
      }
      else
      {
         vbvars = SCIPvarGetVubVars(var);
         coefs = SCIPvarGetVubCoefs(var);
         constants = SCIPvarGetVubConstants(var);
         nvbvars = SCIPvarGetNVubs(var);
      }

      /* loop over all variable lower bounds; a variable lower bound has the form: x >= b*y + d,
       * a variable upper bound the form x <= b*y + d */
      for( n = 0; n < nvbvars; ++n )
      {
         SCIP_VAR* vbvar;
         SCIP_Real coef;
         SCIP_Real constant;

         vbvar = vbvars[n];
         coef = coefs[n];
         constant = constants[n];
         assert(vbvar != NULL);

         /* transform variable bound variable to an active variable, if possible */
         SCIP_CALL( SCIPgetProbvarSum(scip, &vbvar, &coef, &constant) );
         assert(vbvar != NULL);

         if( !SCIPvarIsActive(vbvar) )
            continue;

         /* if the coefficient is positive, the type of bound is the same for the bounded and the bounding variable */
         if( SCIPisPositive(scip, coef) )
            startidx = (lower ? varGetLbIndex(propdata, vbvar) : varGetUbIndex(propdata, vbvar));
         else
            startidx = (lower ? varGetUbIndex(propdata, vbvar) : varGetLbIndex(propdata, vbvar));
         assert(startidx >= 0);

         /* If the vbvar is binary, the vbound should be stored as an implication already.
          * However, it might happen that vbvar was integer when the variable bound was added, but was converted
          * to a binary variable later during presolving when its upper bound was changed to 1. In this case,
          */
         if( SCIPvarGetType(vbvar) == SCIP_VARTYPE_BINARY
            && SCIPvarHasImplic(vbvar, isIndexLowerbound(startidx), var, getBoundtype(v)) )
         {
#if 0
            SCIP_VAR** implvars;
            SCIP_Real* implbounds;
            SCIP_BOUNDTYPE* impltypes;
            int nimplvars;
            int j;
            SCIP_Bool startlower;

            startlower = isIndexLowerbound(startidx);

            implvars = SCIPvarGetImplVars(vbvar, startlower);
            impltypes = SCIPvarGetImplTypes(vbvar, startlower);
            implbounds = SCIPvarGetImplBounds(vbvar, startlower);
            nimplvars = SCIPvarGetNImpls(vbvar, startlower);

            for( j = 0; j < nimplvars; ++j )
            {
               if( (implvars[j] == var) && (lower == (impltypes[j] == SCIP_BOUNDTYPE_LOWER)) )
               {
                  if( lower )
                     assert(SCIPisGE(scip, implbounds[j], (startlower ? coef : 0.0) + constant));
                  else
                     assert(SCIPisLE(scip, implbounds[j], (startlower ? coef : 0.0) + constant));
                  break;
               }
            }
            assert(j  < nimplvars);
#endif
            SCIPdebugMessage("varbound <%s> %s %g * <%s> + %g not added to propagator data due to reverse implication\n",
               SCIPvarGetName(var), (lower ? ">=" : "<="), coef,
               SCIPvarGetName(vbvar), constant);
         }
         else
         {
            SCIP_CALL( addVbound(scip, propdata, startidx, v, coef, constant) );

            SCIPdebugMessage("varbound <%s> %s %g * <%s> + %g added to propagator data\n",
               SCIPvarGetName(var), (lower ? ">=" : "<="), coef,
               SCIPvarGetName(vbvar), constant);

         }
      }
   }

   /* sort the bounds topologically */
   if( propdata->dotoposort )
   {
      SCIP_CALL( topologicalSort(scip, propdata) );
   }

   /* catch variable events */
   SCIP_CALL( catchEvents(scip, propdata) );

   return SCIP_OKAY;
}


/** resolves a propagation by adding the variable which implied that bound change */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable to be reported */
   SCIP_BOUNDTYPE        boundtype,          /**< bound to be reported */
   SCIP_BDCHGIDX*        bdchgidx            /**< the index of the bound change, representing the point of time where the change took place, or NULL for the current local bounds */
   )
{
   assert(propdata != NULL);
   assert(boundtype == SCIP_BOUNDTYPE_LOWER || boundtype == SCIP_BOUNDTYPE_UPPER);

   SCIPdebugMessage(" -> add %s bound of variable <%s> as reason\n",
      getBoundtypeString(boundtype), SCIPvarGetName(var));

   switch( boundtype )
   {
   case SCIP_BOUNDTYPE_LOWER:
      SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
      break;
   case SCIP_BOUNDTYPE_UPPER:
      SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
      break;
   default:
      SCIPerrorMessage("invalid bound type <%d>\n", boundtype);
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** relaxes bound of give variable as long as the given inference bound still leads to a cutoff and add that bound
 *  change to the conflict set
 */
static
SCIP_RETCODE relaxVbdvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which the upper bound should be relaxed */
   SCIP_BOUNDTYPE        boundtype,          /**< boundtype used for the variable bound variable */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place, or NULL for the current local bounds */
   SCIP_Real             relaxedbd           /**< relaxed bound */
   )
{
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPaddConflictRelaxedLb(scip, var, bdchgidx, relaxedbd) );
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      SCIP_CALL( SCIPaddConflictRelaxedUb(scip, var, bdchgidx, relaxedbd) );
   }

   return SCIP_OKAY;
}

/** compute the relaxed bound which is sufficient to propagate the inference lower bound of given variable */
static
SCIP_Real computeRelaxedLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which was propagated */
   SCIP_Real             inferlb,            /**< inference lower bound */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant            /**< inference variable bound constant used */
   )
{
   SCIP_Real relaxedbd;

   if( SCIPvarIsIntegral(var) )
      relaxedbd = (inferlb - 1.0 + 2*SCIPfeastol(scip) - constant) / coef;
   else
      relaxedbd = (inferlb - constant) / coef;

   /* check the computed relaxed lower/upper bound is a proper reason for the inference bound which has to be explained */
   assert(SCIPisEQ(scip, inferlb, SCIPadjustedVarLb(scip, var, relaxedbd * coef + constant)));

   return relaxedbd;
}


/** analyzes an infeasibility which was reached by updating the lower bound of the inference variable above its upper
 *  bound
 */
static
SCIP_RETCODE analyzeConflictLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             infervar,           /**< variable which led to a cutoff */
   SCIP_Real             inferlb,            /**< lower bound which led to infeasibility */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the lower bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the lower bound change */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant,           /**< inference variable bound constant used */
   SCIP_Bool             canwide             /**< can bound widening be used (for vbounds) or not (for inplications or cliques) */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(infervar != NULL);
   assert(SCIPisEQ(scip, SCIPvarGetUbLocal(infervar), SCIPvarGetUbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisEQ(scip, SCIPvarGetUbAtIndex(infervar, NULL, TRUE), SCIPvarGetUbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisGT(scip, inferlb, SCIPvarGetUbLocal(infervar)));
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* check if conflict analysis is applicable */
   if( !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   if( canwide && propdata->usebdwidening )
   {
      SCIP_Real relaxedbd;
      SCIP_Real relaxedub;

      SCIPdebugMessage("try to create conflict using bound widening order: inference variable, variable bound variable\n");

      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );

      /* adjust lower bound */
      inferlb = SCIPadjustedVarLb(scip, infervar, inferlb);

      /* compute a relaxed upper bound which would be sufficient to be still infeasible */
      if( SCIPvarIsIntegral(infervar) )
         relaxedub = inferlb - 1.0;
      else
         relaxedub = inferlb - 2*SCIPfeastol(scip);

      /* try to relax inference variable upper bound such that the infeasibility is still given */
      SCIP_CALL( SCIPaddConflictRelaxedUb(scip, infervar, NULL, relaxedub) );

      /* collect the upper bound which is reported to the conflict analysis */
      relaxedub = SCIPgetConflictVarUb(scip, infervar);

      /* adjust inference bound with respect to the upper bound reported to the conflict analysis */
      if( SCIPvarIsIntegral(infervar) )
         relaxedub = relaxedub + 1.0;
      else
         relaxedub = relaxedub + 2*SCIPfeastol(scip);

      /* compute the relaxed bound which is sufficient to propagate the inference lower bound of given variable */
      relaxedbd = computeRelaxedLowerbound(scip, infervar, relaxedub, coef, constant);

      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, vbdvar, boundtype, NULL, relaxedbd) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }
   else
   {
      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );

      /* add upper bound of the variable for which we tried to change the lower bound */
      SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL) );

      /* add (correct) bound of the variable which let to the new lower bound */
      SCIP_CALL( resolvePropagation(scip, propdata, vbdvar, boundtype, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }

   return SCIP_OKAY;
}

/** compute the relaxed bound which is sufficient to propagate the inference upper bound of given variable */
static
SCIP_Real computeRelaxedUpperbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which was propagated */
   SCIP_Real             inferub,            /**< inference upper bound */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant            /**< inference variable bound constant used */
   )
{
   SCIP_Real relaxedbd;

   if( SCIPvarIsIntegral(var) )
      relaxedbd = (inferub + 1.0 - 2*SCIPfeastol(scip) - constant) / coef;
   else
      relaxedbd = (inferub - constant) / coef;

   /* check the computed relaxed lower/upper bound is a proper reason for the inference bound which has to be explained */
   assert(SCIPisEQ(scip, inferub, SCIPadjustedVarUb(scip, var, relaxedbd * coef + constant)));

   return relaxedbd;
}

/** analyzes an infeasibility which was reached by updating the upper bound of the inference variable below its lower
 *  bound
 */
static
SCIP_RETCODE analyzeConflictUpperbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             infervar,           /**< variable which led to a cutoff */
   SCIP_Real             inferub,            /**< upper bound which led to infeasibility */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the upper bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the upper bound change */
   SCIP_Real             coef,               /**< inference variable bound coefficient used */
   SCIP_Real             constant,           /**< inference variable bound constant used */
   SCIP_Bool             canwide             /**< can bound widening be used (for vbounds) or not (for inplications or cliques) */
   )
{
   assert(scip != NULL);
   assert(propdata != NULL);
   assert(infervar != NULL);
   assert(SCIPisEQ(scip, SCIPvarGetLbLocal(infervar), SCIPvarGetLbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(infervar, NULL, TRUE), SCIPvarGetLbAtIndex(infervar, NULL, FALSE)));
   assert(SCIPisLT(scip, inferub, SCIPvarGetLbLocal(infervar)));
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* check if conflict analysis is applicable */
   if( !SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   if( canwide && propdata->usebdwidening )
   {
      SCIP_Real relaxedbd;
      SCIP_Real relaxedlb;

      SCIPdebugMessage("try to create conflict using bound widening order: inference variable, variable bound variable\n");

      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );

      /* adjust upper bound */
      inferub = SCIPadjustedVarUb(scip, infervar, inferub);

      /* compute a relaxed lower bound which would be sufficient to be still infeasible */
      if( SCIPvarIsIntegral(infervar) )
         relaxedlb = inferub + 1.0;
      else
         relaxedlb = inferub + 2*SCIPfeastol(scip);

      /* try to relax inference variable lower bound such that the infeasibility is still given */
      SCIP_CALL( SCIPaddConflictRelaxedLb(scip, infervar, NULL, relaxedlb) );

      /* collect the lower bound which is reported to the conflict analysis */
      relaxedlb = SCIPgetConflictVarLb(scip, infervar);

      /* adjust inference bound with respect to the upper bound reported to the conflict analysis */
      if( SCIPvarIsIntegral(infervar) )
         relaxedlb = relaxedlb - 1.0;
      else
         relaxedlb = relaxedlb - 2*SCIPfeastol(scip);

      /* compute the relaxed bound which is sufficient to propagate the inference upper bound of given variable */
      relaxedbd = computeRelaxedUpperbound(scip, infervar, relaxedlb, coef, constant);

      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, vbdvar, boundtype, NULL, relaxedbd) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }
   else
   {
      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );

      /* add lower bound of the variable for which we tried to change the upper bound */
      SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL) );

      /* add (correct) bound of the variable which let to the new upper bound */
      SCIP_CALL( resolvePropagation(scip, propdata, vbdvar, boundtype, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
   }

   return SCIP_OKAY;
}


/* tries to tighten the (global) lower bound of the given variable to the given new bound */
static
SCIP_RETCODE tightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable whose lower bound should be tightened */
   SCIP_Real             newlb,              /**< new lower bound for the variable */
   SCIP_Bool             global,             /**< is the bound globally valid? */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the lower bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the lower bound change */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_Real             coef,               /**< coefficient in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Real             constant,           /**< constant in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Bool             canwide,            /**< can bound widening be used (for vbounds) or not (for inplications or cliques) */
   int*                  nchgbds,            /**< pointer to increase, if a bound was changed */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation */
   )
{
   INFERINFO inferinfo;
   SCIP_Real lb;
   SCIP_Bool tightened;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(propdata != NULL);
   assert(var != NULL);
   assert(nchgbds != NULL);
   assert(result != NULL);

   lb = SCIPvarGetLbLocal(var);

   /* check that the new upper bound is better */
   if( (SCIPvarIsIntegral(var) && newlb - lb > 0.5) || (force && SCIPisGT(scip, newlb, lb)) )
      force = TRUE;
   else
      force = FALSE;

   /* try to tighten the lower bound */
   if( global )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newlb, force, &infeasible, &tightened) );
   }
   else
   {
      inferinfo = getInferInfo(boundtype == SCIP_BOUNDTYPE_LOWER ? varGetLbIndex(propdata, vbdvar) : varGetUbIndex(propdata, vbdvar), boundtype);

      SCIP_CALL( SCIPinferVarLbProp(scip, var, newlb, prop, inferInfoToInt(inferinfo), force, &infeasible, &tightened) );
   }

   if( infeasible )
   {
      /* the infeasible results comes from the fact that the new lower bound lies above the current upper bound */
      assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(var)));
      assert(!global || SCIPisGT(scip, newlb, SCIPvarGetUbGlobal(var)));

      SCIPdebugMessage("tightening%s lower bound of variable <%s> to %g due the %s bound of variable <%s> led to infeasibility\n",
         (global ? " global" : ""), SCIPvarGetName(var), newlb, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));

      if( global )
      {
         /* cutoff the root node */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      }
      else
      {
         /* analyzes a infeasibility via conflict analysis */
         SCIP_CALL( analyzeConflictLowerbound(scip, propdata, var, newlb, vbdvar, boundtype, coef, constant, canwide) );
      }
      *result = SCIP_CUTOFF;
   }
   else if( tightened )
   {
      SCIPdebugMessage("tightened%s lower bound of variable <%s> to %g due the %s bound of variable <%s>\n",
         (global ? " global" : ""), SCIPvarGetName(var), newlb, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));
      (*nchgbds)++;
   }

   return SCIP_OKAY;
}

/* tries to tighten the (global) upper bound of the given variable to the given new bound */
static
SCIP_RETCODE tightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var,                /**< variable whose upper bound should be tightened */
   SCIP_Real             newub,              /**< new upper bound of the variable */
   SCIP_Bool             global,             /**< is the bound globally valid? */
   SCIP_VAR*             vbdvar,             /**< variable which is the reason for the upper bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< bound which is the reason for the upper bound change */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_Real             coef,               /**< coefficient in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Real             constant,           /**< constant in vbound constraint causing the propagation;
                                              *   or 0.0 if propagation is caused by clique or implication */
   SCIP_Bool             canwide,            /**< can bound widening be used (for vbounds) or not (for inplications or cliques) */
   int*                  nchgbds,            /**< pointer to increase, if a bound was changed */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation */
   )
{
   INFERINFO inferinfo;
   SCIP_Real ub;
   SCIP_Bool tightened;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(propdata != NULL);
   assert(var != NULL);
   assert(nchgbds != NULL);
   assert(result != NULL);

   ub = SCIPvarGetUbLocal(var);

   /* check that the new upper bound is better */
   if( (SCIPvarIsIntegral(var) && ub - newub > 0.5) || (force && SCIPisLT(scip, newub, ub)) )
      force = TRUE;
   else
      force = FALSE;

   /* try to tighten the upper bound */
   if( global )
   {
      SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newub, force, &infeasible, &tightened) );
   }
   else
   {
      inferinfo = getInferInfo(boundtype == SCIP_BOUNDTYPE_LOWER ? varGetLbIndex(propdata, vbdvar) : varGetUbIndex(propdata, vbdvar), boundtype);

      SCIP_CALL( SCIPinferVarUbProp(scip, var, newub, prop, inferInfoToInt(inferinfo), force, &infeasible, &tightened) );
   }

   if( infeasible )
   {
      /* the infeasible results comes from the fact that the new upper bound lies below the current lower bound */
      assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(var)));
      assert(!global || SCIPisLT(scip, newub, SCIPvarGetLbGlobal(var)));

      SCIPdebugMessage("tightening%s upper bound of variable <%s> to %g due the %s bound of variable <%s> led to infeasibility\n",
         (global ? " global" : ""), SCIPvarGetName(var), newub, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));

      if( global )
      {
         /* cutoff the root node */
         SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
      }
      else
      {
         /* analyzes a infeasibility via conflict analysis */
         SCIP_CALL( analyzeConflictUpperbound(scip, propdata, var, newub, vbdvar, boundtype, coef, constant, canwide) );
      }
      *result = SCIP_CUTOFF;
   }
   else if( tightened )
   {
      SCIPdebugMessage("tightened%s upper bound of variable <%s> to %g due the %s bound of variable <%s>\n",
         (global ? " global" : ""), SCIPvarGetName(var), newub, getBoundtypeString(boundtype), SCIPvarGetName(vbdvar));
      (*nchgbds)++;
   }

   return SCIP_OKAY;
}

/** performs propagation of variables lower and upper bounds, implications, and cliques */
static
SCIP_RETCODE propagateVbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< vbounds propagator */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_BOUNDTYPE starttype;
   SCIP_Real startbound;
   SCIP_Real globalbound;
   int startpos;
   int topopos;
   int v;
   int n;
   int nchgbds;
   int nbounds;
   SCIP_Bool lower;
   SCIP_Bool global;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(result != NULL);

   (*result) = SCIP_DIDNOTRUN;

   /* we do not run the propagator in presolving, because we want to avoid doing the expensive creation of the graph twice */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      return SCIP_OKAY;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* initialize propagator data needed for propagation, if not done yet */
   if( !propdata->initialized )
   {
      SCIP_CALL( initData(scip, prop) );
   }
   assert(propdata->nbounds == 0 || propdata->propqueue != NULL);

   vars = propdata->vars;
   nbounds = propdata->nbounds;

   if( nbounds == 0 )
      return SCIP_OKAY;

   /* propagate all variables if we are in repropagation */
   if( SCIPinRepropagation(scip) )
   {
      SCIP_VAR* var;
      int idx;

      for( v = nbounds - 1; v >= 0; --v )
      {
         idx = propdata->topoorder[v];
         if( idx != -1 && !propdata->inqueue[v] )
         {
            var = vars[getVarIndex(idx)];
            lower = isIndexLowerbound(idx);
            if( !SCIPvarIsBinary(var) || (lower && SCIPvarGetLbLocal(var) > 0.5)
                  || (!lower && SCIPvarGetUbLocal(var) < 0.5) )
            {
               SCIP_CALL( SCIPpqueueInsert(propdata->propqueue, (void*)(size_t)(v + 1)) );
               propdata->inqueue[v] = TRUE;
            }
         }
      }
   }

   /* return if no bound changes are in the priority queue (no changed bounds to handle since last propagation) */
   if( SCIPpqueueNElems(propdata->propqueue) == 0 )
   {
      (*result) = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   nchgbds = 0;

   SCIPdebugMessage("varbound propagator: %d elements in the propagation queue\n", SCIPpqueueNElems(propdata->propqueue));

   /* get variable bound of highest priority from priority queue and try to deduce bound changes for other variables;
    * the priority queue is ordered w.r.t the topological sort of the varbound graph
    */
   while( SCIPpqueueNElems(propdata->propqueue) > 0 )
   {
      topopos = ((int)(size_t)SCIPpqueueRemove(propdata->propqueue)) - 1;
      assert(propdata->inqueue[topopos]);
      startpos = propdata->topoorder[topopos];
      assert(startpos >= 0);
      propdata->inqueue[topopos] = FALSE;

      startvar = vars[getVarIndex(startpos)];
      starttype = getBoundtype(startpos);
      lower = (starttype == SCIP_BOUNDTYPE_LOWER);
      startbound = ( lower ? SCIPvarGetLbLocal(startvar) : SCIPvarGetUbLocal(startvar) );
      globalbound = ( lower ? SCIPvarGetLbGlobal(startvar) : SCIPvarGetUbGlobal(startvar));
      global = SCIPisEQ(scip, startbound, globalbound);

      SCIPdebugMessage("propagate new %s bound of %g of variable <%s>:\n",
         getBoundtypeString(starttype), startbound, SCIPvarGetName(startvar));

      /* there should be neither implications nor cliques for non-binary variables */
      assert(SCIPvarIsBinary(startvar) || SCIPvarGetNImpls(startvar, lower) == 0);
      assert(SCIPvarIsBinary(startvar) || SCIPvarGetNCliques(startvar, lower) == 0);

      if( SCIPvarIsBinary(startvar) )
      {
         /* we only propagate binary variables if the lower bound changed to 1.0 or the upper bound changed to 0.0 */
         if( lower != (startbound > 0.5) )
            continue;

         /* propagate implications */
         if( propdata->useimplics )
         {
            int nimplvars;

            /* if the lower bound of the startvar was changed, it was fixed to 1.0, otherwise it was fixed to 0.0;
             * get all implications for this varfixing
             */
            nimplvars = SCIPvarGetNImpls(startvar, lower);

            /* if there are implications for the varfixing, propagate them */
            if( nimplvars > 0 )
            {
               SCIP_VAR** implvars;
               SCIP_BOUNDTYPE* impltypes;
               SCIP_Real* implbounds;
               int* implids;

               implvars = SCIPvarGetImplVars(startvar, lower);
               impltypes = SCIPvarGetImplTypes(startvar, lower);
               implbounds = SCIPvarGetImplBounds(startvar, lower);
               implids = SCIPvarGetImplIds(startvar, lower);

               for( n = 0; n < nimplvars; ++n )
               {
                  /* implication is just a shortcut, so we do not propagate it now,
                   * because we will propagate the longer way, anyway
                   */
                  if( implids[n] < 0 )
                     continue;

                  /* it might happen that implications point to inactive variables (normally, those are removed when a
                   * variable becomes inactive, but in some cases, it cannot be done), we have to ignore these variables
                   */
                  if( !SCIPvarIsActive(implvars[n]) )
                     continue;

                  if( impltypes[n] == SCIP_BOUNDTYPE_LOWER )
                  {
                     SCIP_CALL( tightenVarLb(scip, prop, propdata, implvars[n], implbounds[n], global, startvar,
                           starttype, force, 0.0, 0.0, FALSE, &nchgbds, result) );
                  }
                  else
                  {
                     SCIP_CALL( tightenVarUb(scip, prop, propdata, implvars[n], implbounds[n], global, startvar,
                           starttype, force, 0.0, 0.0, FALSE, &nchgbds, result) );
                  }

                  if( *result == SCIP_CUTOFF )
                     return SCIP_OKAY;
               }
            }
         }

         /* propagate cliques */
         if( propdata->usecliques )
         {
            int ncliques;

            /* if the lower bound of the startvar was changed, it was fixed to 1.0, otherwise it was fixed to 0.0;
             * get all cliques for this varfixing
             */
            ncliques = SCIPvarGetNCliques(startvar, lower);

            /* if there are cliques for the varfixing, propagate them */
            if( ncliques > 0 )
            {
               SCIP_CLIQUE** cliques;
               int j;

               cliques = SCIPvarGetCliques(startvar, lower);

               for( j = 0; j < ncliques; ++j )
               {
                  SCIP_VAR** cliquevars;
                  SCIP_Bool* cliquevals;
                  int ncliquevars;

                  cliquevars = SCIPcliqueGetVars(cliques[j]);
                  cliquevals = SCIPcliqueGetValues(cliques[j]);
                  ncliquevars = SCIPcliqueGetNVars(cliques[j]);

                  /* fix all variables except for the startvar to the value which is not in the clique */
                  for( n = 0; n < ncliquevars; ++n )
                  {
                     if( cliquevars[n] == startvar )
                        continue;

                     /* try to tighten the bound */
                     if( cliquevals[n] )
                     {
                        /* unnegated variable is in clique, so it has to be fixed to 0.0 */
                        SCIP_CALL( tightenVarUb(scip, prop, propdata, cliquevars[n], 0.0, global, startvar, starttype,
                              force, 0.0, 0.0, FALSE, &nchgbds, result) );
                     }
                     else
                     {
                        /* negated variable is in clique, so it has to be fixed to 1.0 */
                        SCIP_CALL( tightenVarLb(scip, prop, propdata, cliquevars[n], 1.0, global, startvar, starttype,
                              force, 0.0, 0.0, FALSE, &nchgbds, result) );
                     }
                     if( *result == SCIP_CUTOFF )
                        return SCIP_OKAY;
                  }
               }
            }
         }
      }

      /* propagate vbounds */
      if( propdata->usevbounds )
      {
         SCIP_VAR* boundedvar;
         SCIP_Real newbound;
         SCIP_Real coef;
         SCIP_Real constant;

         /* iterate over all vbounds for the given bound */
         for( n = 0; n < propdata->nvbounds[startpos]; ++n )
         {
            boundedvar = vars[getVarIndex(propdata->vboundboundedidx[startpos][n])];
            coef = propdata->vboundcoefs[startpos][n];
            constant = propdata->vboundconstants[startpos][n];

            /* compute new bound */
            newbound = startbound * coef + constant;

            /* try to tighten the bound */
            if( isIndexLowerbound(propdata->vboundboundedidx[startpos][n]) )
            {
               SCIP_CALL( tightenVarLb(scip, prop, propdata, boundedvar, newbound, global, startvar, starttype, force,
                     coef, constant, TRUE, &nchgbds, result) );
            }
            else
            {
               SCIP_CALL( tightenVarUb(scip, prop, propdata, boundedvar, newbound, global, startvar, starttype, force,
                     coef, constant, TRUE, &nchgbds, result) );
            }

            if( *result == SCIP_CUTOFF )
               return SCIP_OKAY;
         }
      }
   }

   SCIPdebugMessage("tightened %d variable bounds\n", nchgbds);

   /* set the result depending on whether bound changes were found or not */
   if( nchgbds > 0 )
      (*result) = SCIP_REDUCEDDOM;
   else
      (*result) = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/**@name Callback methods of propagator
 *
 * @{
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyVbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropVbounds(scip) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);

   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int v;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* free data stored for propagation */
   if( propdata->initialized )
   {
      /* drop all variable events */
      SCIP_CALL( dropEvents(scip, propdata) );

      /* release all variables */
      for( v = 0; v < propdata->nbounds; ++v )
      {
         /* free vbound data */
         if( propdata->vboundsize[v] > 0 )
         {
            SCIPfreeMemoryArray(scip, &propdata->vboundboundedidx[v]);
            SCIPfreeMemoryArray(scip, &propdata->vboundcoefs[v]);
            SCIPfreeMemoryArray(scip, &propdata->vboundconstants[v]);
         }
      }

      /* free priority queue */
      SCIPpqueueFree(&propdata->propqueue);

      /* free arrays */
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundsize, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->nvbounds, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundconstants, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundcoefs, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->vboundboundedidx, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->inqueue, propdata->nbounds);
      SCIPfreeBlockMemoryArray(scip, &propdata->topoorder, propdata->nbounds);

      /* free variable array and hashmap */
      SCIPhashmapFree(&propdata->varhashmap);
      SCIPfreeBlockMemoryArray(scip, &propdata->vars, propdata->nbounds / 2);
   }

   /* reset propagation data */
   resetPropdata(propdata);

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecVbounds)
{  /*lint --e{715}*/

   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, FALSE, result) );

   assert((*result) == SCIP_CUTOFF || (*result) == SCIP_DIDNOTRUN
      || (*result) == SCIP_DIDNOTFIND || (*result) == SCIP_REDUCEDDOM);

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropVbounds)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR* startvar;
   SCIP_BOUNDTYPE starttype;
   int pos;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   starttype = inferInfoGetBoundtype(intToInferInfo(inferinfo));
   pos = inferInfoGetPos(intToInferInfo(inferinfo));
   assert(pos >= 0);
   assert(pos < propdata->nbounds);

   vars = propdata->vars;
   assert(vars != NULL);
   startvar = vars[getVarIndex(pos)];
   assert(startvar != NULL);
   assert(startvar != infervar);

   SCIPdebugMessage("explain %s bound change of variable <%s>\n",
      getBoundtypeString(boundtype), SCIPvarGetName(infervar));

   if( !SCIPvarIsBinary(startvar) && propdata->usebdwidening )
   {
      int* vboundboundedidx;
      SCIP_Real constant;
      SCIP_Real coef;
      int inferidx;
      int nvbounds;
      int b;

      nvbounds = propdata->nvbounds[pos];
      vboundboundedidx = propdata->vboundboundedidx[pos];

      inferidx = boundtype == SCIP_BOUNDTYPE_LOWER ? varGetLbIndex(propdata, infervar) : varGetUbIndex(propdata, infervar);
      assert(inferidx >= 0);

      for( b = 0; b < nvbounds; ++b )
      {
         if( vboundboundedidx[b] == inferidx )
            break;
      }
      assert(b < nvbounds);

      coef = propdata->vboundcoefs[pos][b];
      constant = propdata->vboundconstants[pos][b];
      assert(!SCIPisZero(scip, coef));

      /* compute the relaxed bound which is sufficient to propagate the inference bound of given variable */
      if( boundtype == SCIP_BOUNDTYPE_LOWER )
         relaxedbd = computeRelaxedLowerbound(scip, infervar, relaxedbd, coef, constant);
      else
         relaxedbd = computeRelaxedUpperbound(scip, infervar, relaxedbd, coef, constant);

      /* try to relax variable bound variable */
      SCIP_CALL( relaxVbdvar(scip, startvar, starttype, bdchgidx, relaxedbd) );
   }
   else
   {
      SCIP_CALL( resolvePropagation(scip, propdata, startvar, starttype, bdchgidx) );
   }

   (*result) = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVbound)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int idx;

   assert(eventhdlr != NULL);

   propdata = (SCIP_PROPDATA*)SCIPeventhdlrGetData(eventhdlr);
   assert(propdata != NULL);

   idx = (int) (size_t) eventdata;
   assert(idx >= 0);

   SCIPdebugMessage("eventexec (type=%u): try to add sort index %d: %s(%s) to priority queue\n", SCIPeventGetType(event),
      idx, indexGetBoundString(propdata->topoorder[idx]),
      SCIPvarGetName(propdata->vars[getVarIndex(propdata->topoorder[idx])]));

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_GUBCHANGED && SCIPvarIsBinary(SCIPeventGetVar(event))
      && SCIPeventGetNewbound(event) > 0.5 )
      return SCIP_OKAY;

   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_GLBCHANGED && SCIPvarIsBinary(SCIPeventGetVar(event))
      && SCIPeventGetNewbound(event) < 0.5 )
      return SCIP_OKAY;

   assert(getVarIndex(propdata->topoorder[idx]) < SCIPgetNVars(scip));
   assert(SCIPvarGetType(propdata->vars[getVarIndex(propdata->topoorder[idx])]) != SCIP_VARTYPE_BINARY
      || (isIndexLowerbound(propdata->topoorder[idx]) == (SCIPeventGetNewbound(event) > 0.5)));

   /* add the bound change to the propagation queue, if it is not already contained */
   if( !propdata->inqueue[idx] )
   {
      SCIP_CALL( SCIPpqueueInsert(propdata->propqueue, (void*)(size_t)(idx + 1)) );
      propdata->inqueue[idx] = TRUE;
   }
   assert(SCIPpqueueNElems(propdata->propqueue) > 0);

   return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the vbounds propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create pseudoobj propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );

   /*  reset propagation data */
   resetPropdata(propdata);

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecVbounds, propdata) );
   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyVbounds) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeVbounds) );
   SCIP_CALL( SCIPsetPropExitsol(scip, prop, propExitsolVbounds) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropVbounds) );

   /* include event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &propdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecVbound, (SCIP_EVENTHDLRDATA*)propdata) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/usebdwidening", "should bound widening be used to initialize conflict analysis?",
         &propdata->usebdwidening, FALSE, DEFAULT_USEBDWIDENING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/useimplics", "should implications be propagated?",
         &propdata->useimplics, FALSE, DEFAULT_USEIMPLICS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/usecliques", "should cliques be propagated?",
         &propdata->usecliques, FALSE, DEFAULT_USECLIQUES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/usevbounds", "should vbounds be propagated?",
         &propdata->usevbounds, FALSE, DEFAULT_USEVBOUNDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/dotoposort", "should the bounds be topologically sorted in advance?",
         &propdata->dotoposort, FALSE, DEFAULT_DOTOPOSORT, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/"PROP_NAME"/sortcliques", "should cliques be regarded for the topological sort?",
         &propdata->sortcliques, FALSE, DEFAULT_SORTCLIQUES, NULL, NULL) );

   return SCIP_OKAY;
}

/** returns TRUE if the propagator has the status that all variable lower and upper bounds are propgated */
SCIP_Bool SCIPisPropagatedVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   return (SCIPpqueueNElems(propdata->propqueue) == 0);
}

/** performs propagation of variables lower and upper bounds */
SCIP_RETCODE SCIPexecPropVbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             force,              /**< should domain changes for continuous variables be forced */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   SCIP_PROP* prop;

   prop = SCIPfindProp(scip, PROP_NAME);
   assert(prop != NULL);

   /* perform variable lower and upper bound propagation */
   SCIP_CALL( propagateVbounds(scip, prop, force, result) );

   assert((*result) == SCIP_CUTOFF || (*result) == SCIP_DIDNOTRUN
      || (*result) == SCIP_DIDNOTFIND || (*result) == SCIP_REDUCEDDOM);

   return SCIP_OKAY;
}

/**@} */
