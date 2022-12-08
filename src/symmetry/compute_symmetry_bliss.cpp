/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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

/**@file   compute_symmetry_bliss.cpp
 * @brief  interface for symmetry computations to bliss
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* include bliss graph */
#include <bliss/defs.hh>
#include <bliss/graph.hh>

#include <string.h>
#include <vector>
#include <list>
#include <math.h>

#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_pow.h"
#include "scip/expr.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"

using std::vector;


/** struct for bliss callback */
struct BLISS_Data
{
   SCIP*                 scip;               /**< SCIP pointer */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nperms;             /**< number of permutations */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
   int                   nmaxperms;          /**< maximal number of permutations */
   int                   maxgenerators;      /**< maximal number of generators constructed (= 0 if unlimited) */
};

/* ------------------- map for operator types ------------------- */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyOptype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare the types of two operators according to their name, level and, in case of power, exponent.
 */
static
SCIP_DECL_HASHKEYEQ(SYMhashKeyEQOptype)
{
   SYM_OPTYPE* k1;
   SYM_OPTYPE* k2;

   k1 = (SYM_OPTYPE*) key1;
   k2 = (SYM_OPTYPE*) key2;

   /* first check operator name */
   if ( SCIPexprGetHdlr(k1->expr) != SCIPexprGetHdlr(k2->expr) )
      return FALSE;

   /* for pow expressions, also check exponent (TODO should that happen for signpow as well?) */
   if ( SCIPisExprPower((SCIP*)userptr, k1->expr )
      && SCIPgetExponentExprPow(k1->expr) != SCIPgetExponentExprPow(k2->expr) )  /*lint !e777*/
      return FALSE;

   /* if still undecided, take level */
   if ( k1->level != k2->level )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValOptype)
{  /*lint --e{715}*/
   SYM_OPTYPE* k;
   SCIP_Real exponent;

   k = (SYM_OPTYPE*) key;

   if ( SCIPisExprPower((SCIP*)userptr, k->expr) )
      exponent = SCIPgetExponentExprPow(k->expr);
   else
      exponent = 1.0;

   return SCIPhashThree(SCIPrealHashCode(exponent), k->level,
      SCIPhashKeyValString(NULL, static_cast<void*>(const_cast<char*>(SCIPexprhdlrGetName(SCIPexprGetHdlr(k->expr))))));
}

/* ------------------- map for constant types ------------------- */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyConsttype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare two constants according to their values.
 */
static
SCIP_DECL_HASHKEYEQ(SYMhashKeyEQConsttype)
{
   SYM_CONSTTYPE* k1;
   SYM_CONSTTYPE* k2;

   k1 = (SYM_CONSTTYPE*) key1;
   k2 = (SYM_CONSTTYPE*) key2;

   return (SCIP_Bool)(k1->value == k2->value);  /*lint !e777*/
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValConsttype)
{  /*lint --e{715}*/
   SYM_CONSTTYPE* k;

   k = (SYM_CONSTTYPE*) key;

   return SCIPrealHashCode(k->value);
}

/* ------------------- map for constraint side types ------------------- */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(SYMhashGetKeyRhstype)
{  /*lint --e{715}*/
   return elem;
}

/** returns TRUE iff both keys are equal
 *
 *  Compare two constraint sides according to lhs and rhs.
 */
static
SCIP_DECL_HASHKEYEQ(SYMhashKeyEQRhstype)
{
   SYM_RHSTYPE* k1;
   SYM_RHSTYPE* k2;

   k1 = (SYM_RHSTYPE*) key1;
   k2 = (SYM_RHSTYPE*) key2;

   if ( k1->lhs != k2->lhs )  /*lint !e777*/
      return FALSE;

   return (SCIP_Bool)(k1->rhs == k2->rhs);  /*lint !e777*/
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(SYMhashKeyValRhstype)
{  /*lint --e{715}*/
   SYM_RHSTYPE* k;

   k = (SYM_RHSTYPE*) key;

   return SCIPhashTwo(SCIPrealHashCode(k->lhs), SCIPrealHashCode(k->rhs));
}

/** callback function for bliss */
static
void blisshook(
   void*                 user_param,         /**< parameter supplied at call to bliss */
   unsigned int          n,                  /**< size of aut vector */
   const unsigned int*   aut                 /**< automorphism */
   )
{
   assert( aut != NULL );
   assert( user_param != NULL );

   BLISS_Data* data = static_cast<BLISS_Data*>(user_param);
   assert( data->scip != NULL );
   assert( data->npermvars < (int) n );
   assert( data->maxgenerators >= 0);

   /* make sure we do not generate more that maxgenerators many permutations, if the limit in bliss is not available */
   if ( data->maxgenerators != 0 && data->nperms >= data->maxgenerators )
      return;

   /* copy first part of automorphism */
   bool isIdentity = true;
   int* p = 0;
   if ( SCIPallocBlockMemoryArray(data->scip, &p, data->npermvars) != SCIP_OKAY )
      return;

   for (int j = 0; j < data->npermvars; ++j)
   {
      /* convert index of variable-level 0-nodes to variable indices */
      p[j] = (int) aut[j];
      if ( p[j] != j )
         isIdentity = false;
   }

   /* ignore trivial generators, i.e. generators that only permute the constraints */
   if ( isIdentity )
   {
      SCIPfreeBlockMemoryArray(data->scip, &p, data->npermvars);
      return;
   }

   /* check whether we should allocate space for perms */
   if ( data->nmaxperms <= 0 )
   {
      if ( data->maxgenerators == 0 )
         data->nmaxperms = 100;   /* seems to cover many cases */
      else
         data->nmaxperms = data->maxgenerators;

      if ( SCIPallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms) != SCIP_OKAY )
         return;
   }
   else if ( data->nperms >= data->nmaxperms )    /* check whether we need to resize */
   {
      int newsize = SCIPcalcMemGrowSize(data->scip, data->nperms + 1);
      assert( newsize >= data->nperms );
      assert( data->maxgenerators == 0 );

      if ( SCIPreallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms, newsize) != SCIP_OKAY )
         return;

      data->nmaxperms = newsize;
   }

   data->perms[data->nperms++] = p;
}

/** callback function for bliss in case of reflection symmetries */
static
void blisshookReflSym(
   void*                 user_param,         /**< parameter supplied at call to bliss */
   unsigned int          n,                  /**< size of aut vector */
   const unsigned int*   aut                 /**< automorphism */
   )
{
   assert( aut != NULL );
   assert( user_param != NULL );

   BLISS_Data* data = static_cast<BLISS_Data*>(user_param);
   assert( data->scip != NULL );
   assert( 2 * data->npermvars < (int) n );
   assert( data->maxgenerators >= 0);

   /* make sure we do not generate more that maxgenerators many permutations, if the limit in bliss is not available */
   if ( data->maxgenerators != 0 && data->nperms >= data->maxgenerators )
      return;

   /* copy first part of automorphism */
   bool isIdentity = true;
   int* p = 0;
   if ( SCIPallocBlockMemoryArray(data->scip, &p, 2 * data->npermvars) != SCIP_OKAY )
      return;

   for (int j = 0; j < 2 * data->npermvars; ++j)
   {
      /* convert index of variable-level 0-nodes to variable indices */
      p[j] = (int) aut[j];
      if ( p[j] != j )
         isIdentity = false;
   }

   /* ignore trivial generators, i.e. generators that only permute the constraints */
   if ( isIdentity )
   {
      SCIPfreeBlockMemoryArray(data->scip, &p, 2 * data->npermvars);
      return;
   }

   /* check whether we should allocate space for perms */
   if ( data->nmaxperms <= 0 )
   {
      if ( data->maxgenerators == 0 )
         data->nmaxperms = 100;   /* seems to cover many cases */
      else
         data->nmaxperms = data->maxgenerators;

      if ( SCIPallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms) != SCIP_OKAY )
         return;
   }
   else if ( data->nperms >= data->nmaxperms )    /* check whether we need to resize */
   {
      int newsize = SCIPcalcMemGrowSize(data->scip, data->nperms + 1);
      assert( newsize >= data->nperms );
      assert( data->maxgenerators == 0 );

      if ( SCIPreallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms, newsize) != SCIP_OKAY )
         return;

      data->nmaxperms = newsize;
   }

   data->perms[data->nperms++] = p;
}

/** Creates the nodes in the graph that correspond to variables. Each variable type gets a unique color
 *
 *  @pre graph should be empty when this is called
 */
static
SCIP_RETCODE createVariableNodes(
   SCIP*                 scip,               /**< SCIP instance */
   bliss::Graph*         G,                  /**< Graph to be constructed */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix (also contains the relevant variables) */
   int&                  nnodes,             /**< buffer to store number of nodes in graph */
   const int&            nedges,             /**< buffer to store number of edges in graph */
   int&                  nusedcolors         /**< buffer to store number of used colors */
   )
{
   assert( scip != NULL );
   assert( G != NULL );
   assert( nnodes == 0 );
   assert( nedges == 0 );
   assert( nusedcolors == 0 );
   SCIPdebugMsg(scip, "Creating graph with colored nodes for variables.\n");

   /* add nodes for variables */
   for (int v = 0; v < matrixdata->npermvars; ++v)
   {
      const int color = matrixdata->permvarcolors[v];
      assert( 0 <= color && color < matrixdata->nuniquevars );

#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) color);
      assert( node == v );
#else
      (void) G->add_vertex((unsigned) color);
#endif

      ++nnodes;
   }

   /* this is not exactly true, since we skip auxvars, but it doesn't matter if some colors are not used at all */
   nusedcolors = matrixdata->nuniquevars;

   return SCIP_OKAY;
}

/** Construct linear part of colored graph for symmetry computations
 *
 *  Construct graph:
 *  - Each variable gets a different node.
 *  - Each constraint gets a different node.
 *  - Each matrix coefficient gets a different node that is connected to the two nodes
 *    corresponding to the respective constraint and variable.
 *
 *  Each different variable, rhs, matrix coefficient gets a different color that is attached to the corresponding entries.
 *
 *  @pre This method assumes that the nodes corresponding to permutation variables are already in the graph and that
 *  their node number is equal to their index.
 */
static
SCIP_RETCODE fillGraphByLinearConss(
   SCIP*                 scip,               /**< SCIP instance */
   bliss::Graph*         G,                  /**< Graph to be constructed */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   int&                  nnodes,             /**< buffer to store number of nodes in graph */
   int&                  nedges,             /**< buffer to store number of edges in graph */
   int&                  nusedcolors,        /**< buffer to store number of used colors */
   SCIP_Bool&            success             /**< whether the construction was successful */
   )
{
   assert( nnodes == (int) G->get_nof_vertices() );
   assert( nusedcolors <= nnodes );

   SCIPdebugMsg(scip, "Filling graph with colored coefficient nodes for linear part.\n");

   success = TRUE;

   /* add nodes for rhs of constraints */
   for (int c = 0; c < matrixdata->nrhscoef; ++c)
   {
      const int color = matrixdata->rhscoefcolors[c];
      assert( 0 <= color && color < matrixdata->nuniquerhs );

#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) (nusedcolors + color));
      assert( node == matrixdata->npermvars + c );
#else
      (void) G->add_vertex((unsigned) (nusedcolors + color));
#endif

      ++nnodes;
   }
   assert( (int) G->get_nof_vertices() == matrixdata->npermvars + matrixdata->nrhscoef );
   nusedcolors += matrixdata->nuniquerhs;

   /* Grouping of nodes depends on the number of nodes in the bipartite graph class.
    * If there are more variables than constraints, we group by constraints.
    * That is, given several variable nodes which are incident to one constraint node by the same color,
    * we join these variable nodes to the constraint node by only one intermediate node.
    */
   const bool groupByConstraints = matrixdata->nrhscoef < matrixdata->npermvars;
   if ( groupByConstraints )
      SCIPdebugMsg(scip, "Group intermediate nodes by constraints.\n");
   else
      SCIPdebugMsg(scip, "Group intermediate nodes by variables.\n");

   /* "colored" edges based on all matrix coefficients - loop through ordered matrix coefficients */
   int ninternodes;
   if ( groupByConstraints )
      ninternodes = matrixdata->nrhscoef;
   else
      ninternodes = matrixdata->npermvars;

   int* internodes = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &internodes, ninternodes) ); /*lint !e530*/
   for (int l = 0; l < ninternodes; ++l)
      internodes[l] = -1;

   /* We pass through the matrix coeficients, grouped by color, i.e., different coefficients. If the coeffients appear
    * in the same row or column, it suffices to only generate a single node (depending on groupByConstraints). We store
    * this node in the array internodes. In order to avoid reinitialization, we store the node number with increasing
    * numbers for each color. The smallest number for the current color is stored in firstcolornodenumber. */
   int oldcolor = -1;
#ifndef NDEBUG
   SCIP_Real oldcoef = SCIP_INVALID;
#endif
   int firstcolornodenumber = -1;
   for (int j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int idx = matrixdata->matidx[j];
      assert( 0 <= idx && idx < matrixdata->nmatcoef );

      /* find color corresponding to matrix coefficient */
      const int color = matrixdata->matcoefcolors[idx];
      assert( 0 <= color && color < matrixdata->nuniquemat );

      assert( 0 <= matrixdata->matrhsidx[idx] && matrixdata->matrhsidx[idx] < matrixdata->nrhscoef );
      assert( 0 <= matrixdata->matvaridx[idx] && matrixdata->matvaridx[idx] < matrixdata->npermvars );

      const int rhsnode = matrixdata->npermvars + matrixdata->matrhsidx[idx];
      const int varnode = matrixdata->matvaridx[idx];
      assert( matrixdata->npermvars <= rhsnode && rhsnode < matrixdata->npermvars + matrixdata->nrhscoef );
      assert( rhsnode < (int) G->get_nof_vertices() );
      assert( varnode < (int) G->get_nof_vertices() );

      /* if we have only one color, we do not need intermediate nodes */
      if ( matrixdata->nuniquemat == 1 )
      {
         G->add_edge((unsigned) varnode, (unsigned) rhsnode);
         ++nedges;
      }
      else
      {
         /* if new group of coefficients has been reached */
         if ( color != oldcolor )
         {
            assert( ! SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );
            oldcolor = color;
            firstcolornodenumber = nnodes;
#ifndef NDEBUG
            oldcoef = matrixdata->matcoef[idx];
#endif
         }
         else
            assert( SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );

         int varrhsidx;
         if ( groupByConstraints )
            varrhsidx = matrixdata->matrhsidx[idx];
         else
            varrhsidx = matrixdata->matvaridx[idx];
         assert( 0 <= varrhsidx && varrhsidx < ninternodes );

         if ( internodes[varrhsidx] < firstcolornodenumber )
         {
            internodes[varrhsidx] = (int) G->add_vertex((unsigned) (nusedcolors + color));
            ++nnodes;
         }
         assert( internodes[varrhsidx] >= matrixdata->npermvars + matrixdata->nrhscoef );
         assert( internodes[varrhsidx] >= firstcolornodenumber );

         /* determine whether graph would be too large for bliss (can only handle int) */
         if ( nnodes >= INT_MAX/2 )
         {
            success = FALSE;
            break;
         }

         G->add_edge((unsigned) varnode, (unsigned) internodes[varrhsidx]);
         G->add_edge((unsigned) rhsnode, (unsigned) internodes[varrhsidx]);
         nedges += 2;
      }
   }

   nusedcolors += matrixdata->nuniquemat;

   SCIPfreeBufferArray(scip, &internodes);
   return SCIP_OKAY;
}

/** Construct non-linear part of colored graph for symmetry computations
 *
 *  Construct graph:
 *  - Each node of the expression trees gets a different node.
 *  - Each coefficient of a sum expression gets its own node connected to the node of the corresponding child.
 *  - Each constraint (with lhs and (!) rhs) gets its own node connected to the corresponding node of the root expression.
 *
 *  @note: In contrast to the linear part, lhs and rhs are treated together here, so that each unique combination of lhs
 *  and rhs gets its own node. This makes the implementation a lot simpler with the small downside, that different
 *  formulations of the same constraints would not be detected as equivalent, e.g. for
 *      0 <= x1 + x2 <= 1
 *      0 <= x3 + x4
 *           x3 + x4 <= 1
 *  there would be no symmetry between (x1,x2) and (x3,x4) detected.
 *
 *  Each different constraint (sides), sum-expression coefficient, constant and operator type gets a
 *  different color that is attached to the corresponding entries.
 *
 *  @pre This method assumes that the nodes corresponding to permutation variables are already in the graph and that
 *  their node number is equal to their index.
 */
static
SCIP_RETCODE fillGraphByNonlinearConss(
   SCIP*                 scip,               /**< SCIP instance */
   bliss::Graph*         G,                  /**< Graph to be constructed */
   SYM_EXPRDATA*         exprdata,           /**< data for nonlinear constraints */
   int&                  nnodes,             /**< buffer to store number of nodes in graph */
   int&                  nedges,             /**< buffer to store number of edges in graph */
   int&                  nusedcolors,        /**< number of used colors ind the graph so far */
   SCIP_Bool&            success             /**< whether the construction was successful */
   )
{
   SCIP_HASHTABLE* optypemap;
   SCIP_HASHTABLE* consttypemap;
   SCIP_HASHTABLE* sumcoefmap;
   SCIP_HASHTABLE* rhstypemap;
   SYM_OPTYPE* uniqueoparray = NULL;
   SYM_CONSTTYPE* uniqueconstarray = NULL;
   SYM_CONSTTYPE* sumcoefarray = NULL;
   SYM_RHSTYPE* uniquerhsarray = NULL;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   int nuniqueops = 0;
   int nuniqueconsts = 0;
   int nuniquecoefs = 0;
   int nuniquerhs = 0;
   int oparraysize = exprdata->nuniqueoperators;
   int constarraysize = exprdata->nuniqueconstants;
   int coefarraysize = exprdata->nuniquecoefs;
   int rhsarraysize;

   assert( scip != NULL );
   assert( G != NULL );
   assert( exprdata != NULL );
   assert( nnodes == (int) G->get_nof_vertices() );
   assert( nnodes >= nusedcolors );

   success = TRUE; /*lint !e838*/

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   nconss = conshdlr != NULL ? SCIPconshdlrGetNConss(conshdlr) : 0;
   if ( nconss == 0 )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   rhsarraysize = nconss;

   SCIPdebugMsg(scip, "Filling graph with colored coefficient nodes for non-linear part.\n");

   /* create maps for optypes, constants, sum coefficients and rhs to indices */
   SCIP_CALL( SCIPhashtableCreate(&optypemap, SCIPblkmem(scip), oparraysize, SYMhashGetKeyOptype,
         SYMhashKeyEQOptype, SYMhashKeyValOptype, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&consttypemap, SCIPblkmem(scip), constarraysize, SYMhashGetKeyConsttype,
         SYMhashKeyEQConsttype, SYMhashKeyValConsttype, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&sumcoefmap, SCIPblkmem(scip), coefarraysize, SYMhashGetKeyConsttype,
         SYMhashKeyEQConsttype, SYMhashKeyValConsttype, (void*) scip) );
   SCIP_CALL( SCIPhashtableCreate(&rhstypemap, SCIPblkmem(scip), rhsarraysize, SYMhashGetKeyRhstype,
         SYMhashKeyEQRhstype, SYMhashKeyValRhstype, (void*) scip) );

   assert( optypemap != NULL );
   assert( consttypemap != NULL );
   assert( sumcoefmap != NULL );
   assert( rhstypemap != NULL );

   /* allocate space for mappings from optypes, constants, sum coefficients and rhs to colors */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniqueoparray, oparraysize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniqueconstarray, constarraysize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sumcoefarray, coefarraysize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &uniquerhsarray, rhsarraysize) );

   SCIP_EXPRITER* it;
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );

   /* iterate over all expressions and add the corresponding nodes to the graph */
   for (int i = 0; i < nconss; ++i)
   {
      SCIP_EXPR* rootexpr;
      vector<int> visitednodes(0);
      vector<SCIP_Bool> ischildofsum(0);
      int currentlevel = 0;

      rootexpr = SCIPgetExprNonlinear(conss[i]);

      SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
      SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_LEAVEEXPR);

      for (SCIP_EXPR* expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it)) /*lint !e441*/ /*lint !e440*/
      {
         switch( SCIPexpriterGetStageDFS(it) )
         {
            /* upon entering an expression, check its type and add nodes and edges if neccessary */
            case SCIP_EXPRITER_ENTEREXPR:
            {
               int node = -1;
               int parentnode = -1;
               int color = -1;

               /* for variable expressions, get the corresponding node that is already in the graph */
               if ( SCIPisExprVar(scip, expr) )
               {
                  SCIP_VAR* var = SCIPgetVarExprVar(expr);

                  /* check whether the variable is active; if not, then replace the inactive variable by its aggregation
                   * or its fixed value; note that this step is equivalent as representing an inactive variable as sum
                   * expression
                   */
                  if ( SCIPvarIsActive(var) )
                  {
                     node = SCIPvarGetProbindex(var);
                     assert( node < (int) G->get_nof_vertices() );
                  }
                  else
                  {
                     SCIP_VAR** vars = NULL;
                     SCIP_Real* vals = NULL;
                     SCIP_Real constant = 0;
                     int varsize = 1;
                     int requiredsize;
                     int k;

                     SCIP_CALL( SCIPallocBufferArray(scip, &vars, varsize) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &vals, varsize) );

                     vars[0] = var;
                     vals[0] = 1.0;

                     SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, vals, &varsize, varsize, &constant, &requiredsize, TRUE) );

                     if ( requiredsize > varsize )
                     {
                        SCIP_CALL( SCIPreallocBufferArray(scip, &vars, requiredsize) );
                        SCIP_CALL( SCIPreallocBufferArray(scip, &vals, requiredsize) );

                        SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, vals, &varsize, requiredsize, &constant, &requiredsize, TRUE) );
                        assert( requiredsize <= varsize );
                     }

                     parentnode = visitednodes[visitednodes.size() - 1];
                     assert( parentnode < (int) G->get_nof_vertices() );

                     /* create nodes for all aggregation variables and coefficients and connect them to the parent node */
                     for ( k = 0; k < requiredsize; ++k )
                     {
                        SYM_CONSTTYPE* ct;
                        int internode;

                        assert( vars[k] != NULL );
                        assert( vals[k] != 0.0 );
                        assert( nuniquecoefs < coefarraysize );

                        ct = &sumcoefarray[nuniquecoefs];
                        ct->value = vals[k];

                        if ( !SCIPhashtableExists(sumcoefmap, (void *) ct) )
                        {
                           SCIP_CALL( SCIPhashtableInsert(sumcoefmap, (void *) ct) );
                           ct->color = nusedcolors++;
                           color = ct->color;
                           nuniquecoefs++;
                        }
                        else
                        {
                           color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(sumcoefmap, (void *) ct))->color;
                        }

                        /* add the intermediate node with the corresponding color */
                        internode = (int) G->add_vertex((unsigned) color);
                        ++nnodes;

                        assert( internode < (int) G->get_nof_vertices() );

                        G->add_edge((unsigned) internode, (unsigned) parentnode);
                        ++nedges;

                        /* connect the intermediate node to its corresponding variable node */
                        node = SCIPvarGetProbindex(vars[k]);
                        assert( node < (int) G->get_nof_vertices() );

                        G->add_edge((unsigned) node, (unsigned) internode);
                        ++nedges;
                     }

                     /* add the node for the constant */
                     if ( constant != 0.0 )
                     {
                        SYM_CONSTTYPE* ct;

                        /* check whether we have to resize */
                        if ( nuniqueconsts >= constarraysize )
                        {
                           int newsize = SCIPcalcMemGrowSize(scip, nuniqueconsts+1);
                           assert( newsize >= 0 );
                           SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &uniqueconstarray, constarraysize, newsize) );
                           constarraysize = newsize;
                        }

                        assert( nuniqueconsts < constarraysize );

                        ct = &uniqueconstarray[nuniqueconsts];
                        ct->value = constant;

                        if ( !SCIPhashtableExists(consttypemap, (void *) ct) )
                        {
                           SCIP_CALL( SCIPhashtableInsert(consttypemap, (void *) ct) );
                           ct->color = nusedcolors++;
                           color = ct->color;
                           nuniqueconsts++;
                        }
                        else
                        {
                           color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(consttypemap, (void *) ct))->color;
                        }

                        /* add the node with a new color */
                        node = (int) G->add_vertex((unsigned) color);
                        ++nnodes;

                        assert( node < (int) G->get_nof_vertices() );

                        G->add_edge((unsigned) node, (unsigned) parentnode);
                        ++nedges;
                     }

                     SCIPfreeBufferArray(scip, &vals);
                     SCIPfreeBufferArray(scip, &vars);

                     /* add a filler node since it will be removed in the next iteration anyway */
                     visitednodes.push_back(nnodes);
                     ischildofsum.push_back(FALSE);
                     ++currentlevel;

                     break;
                  }
               }
               /* for constant expressions, get the color of its type (value) or assign a new one */
               else if ( SCIPisExprValue(scip, expr) )
               {
                  SYM_CONSTTYPE* ct;

                  assert( nuniqueconsts < constarraysize );

                  ct = &uniqueconstarray[nuniqueconsts];
                  ct->value = SCIPgetValueExprValue(expr);

                  if ( !SCIPhashtableExists(consttypemap, (void *) ct) )
                  {
                     SCIP_CALL( SCIPhashtableInsert(consttypemap, (void *) ct) );
                     ct->color = nusedcolors++;
                     color = ct->color;
                     nuniqueconsts++;
                  }
                  else
                  {
                     color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(consttypemap, (void *) ct))->color;
                  }
               }
               /* for all other expressions, get the color of its operator type or assign a new one */
               else
               {
                  SYM_OPTYPE* ot;

                  assert( nuniqueops < oparraysize );

                  ot = &uniqueoparray[nuniqueops];

                  ot->expr = expr;
                  ot->level = currentlevel;

                  if ( !SCIPhashtableExists(optypemap, (void *) ot) )
                  {
                     SCIP_CALL( SCIPhashtableInsert(optypemap, (void *) ot) );
                     ot->color = nusedcolors++;
                     color = ot->color;
                     nuniqueops++;
                  }
                  else
                  {
                     color = ((SYM_OPTYPE*) SCIPhashtableRetrieve(optypemap, (void *) ot))->color;
                  }
               }

               /* if this is the root expression, add the constraint side node (will be parent of expression node) */
               if ( SCIPexpriterGetParentDFS(it) == NULL )
               {
                  /* add the node corresponding to the constraint */
                  SYM_RHSTYPE* rt;
                  int parentcolor;

                  assert( nuniquerhs < rhsarraysize );

                  rt = &uniquerhsarray[nuniquerhs];
                  rt->lhs = SCIPgetLhsNonlinear(conss[i]);
                  rt->rhs = SCIPgetRhsNonlinear(conss[i]);

                  if ( !SCIPhashtableExists(rhstypemap, (void *) rt) )
                  {
                     SCIP_CALL( SCIPhashtableInsert(rhstypemap, (void *) rt) );
                     rt->color = nusedcolors++;
                     parentcolor = rt->color;
                     nuniquerhs++;
                  }
                  else
                  {
                     parentcolor = ((SYM_RHSTYPE*) SCIPhashtableRetrieve(rhstypemap, (void *) rt))->color;
                  }

                  /* add the constraint side node with the corresponding color */
                  parentnode = (int) G->add_vertex((unsigned) parentcolor);
                  ++nnodes;

                  assert( parentnode < (int) G->get_nof_vertices() );
               }
               /* otherwise, get the parentnode stored in visitednodes */
               else
               {
                  parentnode = visitednodes[visitednodes.size() - 1];
                  assert( parentnode < (int) G->get_nof_vertices() );
               }

               /* in all cases apart from variable expressions, the new node is added with the corresponding color */
               if ( color != -1 )
               {
                  node = (int) G->add_vertex((unsigned) color);
                  ++nnodes;

                  assert( node < (int) G->get_nof_vertices() );
               }

               /* store the new node so that it can be used as parentnode later */
               assert( node != -1 );
               visitednodes.push_back(node);
               ischildofsum.push_back(FALSE);

               /* connect the current node with its parent */
               assert( parentnode != -1 );
               G->add_edge((unsigned) node, (unsigned) parentnode);
               ++nedges;

               /* for sum expression, also add intermediate nodes for the coefficients */
               if ( SCIPisExprSum(scip, expr) )
               {
                  SCIP_Real* coefs = SCIPgetCoefsExprSum(expr);
                  int internode;

                  /* iterate over children from last to first, such that visitednodes array is in correct order */
                  for (int j = SCIPexprGetNChildren(expr) - 1; j >= 0; --j)
                  {
                     SYM_CONSTTYPE* ct;

                     assert( nuniquecoefs < coefarraysize );

                     ct = &sumcoefarray[nuniquecoefs];
                     ct->value = coefs[j];

                     if ( !SCIPhashtableExists(sumcoefmap, (void *) ct) )
                     {
                        SCIP_CALL( SCIPhashtableInsert(sumcoefmap, (void *) ct) );
                        ct->color = nusedcolors++;
                        color = ct->color;
                        nuniquecoefs++;
                     }
                     else
                     {
                        color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(sumcoefmap, (void *) ct))->color;
                     }

                     /* add the intermediate node with the corresponding color */
                     internode = (int) G->add_vertex((unsigned) color);
                     ++nnodes;
                     visitednodes.push_back(internode);
                     ischildofsum.push_back(TRUE);

                     assert( internode < (int) G->get_nof_vertices() );

                     G->add_edge((unsigned) internode, (unsigned) node);
                     ++nedges;
                  }

                  /* add node for the constant term of the sum expression */
                  SCIP_Real constval = SCIPgetConstantExprSum(expr);
                  if ( constval != 0.0 )
                  {
                     SYM_CONSTTYPE* ct;

                     /* check whether we have to resize */
                     if ( nuniqueconsts >= constarraysize )
                     {
                        int newsize = SCIPcalcMemGrowSize(scip, nuniqueconsts+1);
                        assert( newsize >= 0 );
                        SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &uniqueconstarray, constarraysize, newsize) );
                        constarraysize = newsize;
                     }

                     assert( nuniqueconsts < constarraysize );

                     ct = &uniqueconstarray[nuniqueconsts];
                     ct->value = constval;

                     if ( !SCIPhashtableExists(consttypemap, (void *) ct) )
                     {
                        SCIP_CALL( SCIPhashtableInsert(consttypemap, (void *) ct) );
                        ct->color = nusedcolors++;
                        color = ct->color;
                        nuniqueconsts++;
                     }
                     else
                     {
                        color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(consttypemap, (void *) ct))->color;
                     }

                     /* add the node with a new color */
                     internode = (int) G->add_vertex((unsigned) color);
                     ++nnodes;

                     assert( node < (int) G->get_nof_vertices() );

                     G->add_edge((unsigned) internode, (unsigned) node);
                     ++nedges;
                  }
               }

               currentlevel++;
               break;
            }

            /* when leaving an expression, the nodes that are not needed anymore are erased from the respective arrays */
            case SCIP_EXPRITER_LEAVEEXPR:
            {
               visitednodes.pop_back();
               ischildofsum.pop_back();
               currentlevel--;

               /* When leaving the child of a sum expression, we have to pop again to get rid of the intermediate nodes
                * used for the coefficients of summands
                */
               if ( !ischildofsum.empty() && ischildofsum[ischildofsum.size() - 1] )
               {
                  visitednodes.pop_back();
                  ischildofsum.pop_back();
               }

               break;
            }

            default:
               SCIPABORT(); /* we should never be called in this stage */
               break;
         }
      }

      assert( currentlevel == 0 );
      assert( visitednodes.empty() );
      assert( ischildofsum.empty() );

      /* determine whether graph would be too large for bliss (can only handle int) */
      if ( nnodes >= INT_MAX/2 )
      {
         success = FALSE; /*lint !e838*/
         break;
      }
   }

   /* free everything */
   SCIPfreeExpriter(&it);
   SCIPfreeBlockMemoryArrayNull(scip, &uniquerhsarray, rhsarraysize);
   SCIPfreeBlockMemoryArrayNull(scip, &sumcoefarray, coefarraysize);
   SCIPfreeBlockMemoryArrayNull(scip, &uniqueconstarray, constarraysize);
   SCIPfreeBlockMemoryArrayNull(scip, &uniqueoparray, oparraysize);
   SCIPhashtableFree(&rhstypemap);
   SCIPhashtableFree(&sumcoefmap);
   SCIPhashtableFree(&consttypemap);
   SCIPhashtableFree(&optypemap);

   return SCIP_OKAY;
}

/** returns the color of a tree node in its respective class */
int colorInClass(
   int                   nodeidx,            /**< index of node in tree */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data for reflection symmetry detection */
   SCIP_Bool             inverse             /**< whether the color of the negated input shall be returned */
   )
{
   int color;

   assert( reflsymdata != NULL );

   switch ( reflsymdata->trees[nodeidx] )
   {
   case SYM_NODETYPE_OPERATOR :
      color = reflsymdata->opscolors[reflsymdata->treemap[nodeidx]];
      break;
   case SYM_NODETYPE_COEF :
      if ( inverse )
         color = reflsymdata->invcoefcolors[reflsymdata->treemap[nodeidx]];
      else
         color = reflsymdata->coefcolors[reflsymdata->treemap[nodeidx]];
      break;
   case SYM_NODETYPE_VAR :
      if ( inverse )
         color = reflsymdata->invvarcolors[reflsymdata->treemap[nodeidx]];
      else
         color = reflsymdata->varcolors[reflsymdata->treemap[nodeidx]];
      break;
   default :
      assert( reflsymdata->trees[nodeidx] == SYM_NODETYPE_VAL );
      color = reflsymdata->valcolors[reflsymdata->treemap[nodeidx]];
   }

   return color;
}

/** returns the types of flips allowed for variable arguments of an operator */
static
SYM_FLIPTYPE getOperatorFliptype(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data of CIP */
   int                   opidx               /**< index of operator */
   )
{
   SYM_NODETYPE* trees;
   int* treemap;
   SCIP_Real* treevals;
   SCIP_EXPRHDLR** treeops;
   SCIP_EXPRHDLR* op;

   assert( reflsymdata != NULL );
   assert( 0 <= opidx && opidx < reflsymdata->ntrees );
   assert( reflsymdata->trees[opidx] == SYM_NODETYPE_OPERATOR );

   trees = reflsymdata->trees;
   treemap = reflsymdata->treemap;
   treevals = reflsymdata->treevals;
   treeops = reflsymdata->treeops;

   /* operators without proper arguments do not allow any flips */
   if ( opidx >= reflsymdata->ntrees - 2 )
      return SYM_FLIPTYPE_NONE;

   /* check which types of flips are allowed by an operator */
   op = treeops[treemap[opidx]];
   if ( op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_OR || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_XOR ||
      op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_BDDISJ || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_GEQ ||
      op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_OBJ || op == SCIPfindExprhdlr(scip, "sum") )
      return SYM_FLIPTYPE_SHIFT_ODD;
   else if ( treeops[treemap[opidx]] == SCIPfindExprhdlr(scip, "prod") ||
      treeops[treemap[opidx]] == SCIPfindExprhdlr(scip, "sin") )
      return SYM_FLIPTYPE_ODD;
   else if ( op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_SOS1 || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_SOS2 ||
      op ==  SCIPfindExprhdlr(scip, "abs") || op ==  SCIPfindExprhdlr(scip, "cos") )
      return SYM_FLIPTYPE_EVEN;
   else if ( op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_AND || op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_INDICATOR ||
      op == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_EQ )
      return SYM_FLIPTYPE_NONE;
   else if ( treeops[treemap[opidx]] == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_SIGNPOWER )
   {
      if ( trees[opidx + 1] == SYM_NODETYPE_VAL )
      {
         if ( trees[opidx + 2] != SYM_NODETYPE_COEF )
            return SYM_FLIPTYPE_NONE;
         else
            return SYM_FLIPTYPE_EVEN;
      }
      else
      {
         if ( trees[opidx + 1] != SYM_NODETYPE_COEF )
            return SYM_FLIPTYPE_NONE;
         else
            return SYM_FLIPTYPE_EVEN;
      }
   }
   else if ( treeops[treemap[opidx]] == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_POWER )
   {
      /* power 1 yields odd function */
      if ( trees[opidx + 1] != SYM_NODETYPE_VAL )
      {
         /* has not a single variable as argument */
         if ( trees[opidx + 2] != SYM_NODETYPE_COEF )
            return SYM_FLIPTYPE_NONE;
         else
            return SYM_FLIPTYPE_ODD;
      }

      /* check whether exponent leads to even/odd function */
      if ( ! SCIPisIntegral(scip, treevals[treemap[opidx + 1]]) )
         return SYM_FLIPTYPE_NONE;
      else if ( (int) treevals[treemap[opidx + 1]] % 2 == 0 )
         return SYM_FLIPTYPE_EVEN;
      else
         return SYM_FLIPTYPE_ODD;
   }
   else if ( treeops[treemap[opidx]] == (SCIP_EXPRHDLR*) SYM_CONSOPTYPE_BIPROD )
      return SYM_FLIPTYPE_BIPROD;

   return SYM_FLIPTYPE_NONE;
}

/** creates the graph for detecting reflection symmetries */
SCIP_RETCODE SCIPcreateReflectionSymmetryDetectionGraph(
   SCIP*                 scip,               /**< SCIP pointer */
   bliss::Graph*         G,                  /**< graph to be constructed */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data of CIP */
   SCIP_Bool&            success             /**< whether the construction was successful */
   )
{
   int nnodes = 0;
   int nusedcolors = 0;

   assert( scip != NULL );
   assert( G != NULL );
   assert( reflsymdata != NULL );

   success = TRUE;

   SCIPdebugMsg(scip, "Add variable nodes to symmetry detection graph.\n");

   /* add nodes for variables */
   for (int v = 0; v < reflsymdata->ntreevars; ++v)
   {
      const int color = reflsymdata->varcolors[v];
      assert( 0 <= color && color < reflsymdata->nuniquevars );

#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) color);
      assert( node == v );
#else
      (void) G->add_vertex((unsigned) color);
#endif

      ++nnodes;
   }

   /* add nodes for negated variables and add edges to the original variable*/
   for (int v = 0; v < reflsymdata->ntreevars; ++v)
   {
      const int color = reflsymdata->invvarcolors[v];
      assert( 0 <= color && color < reflsymdata->nuniquevars );

      int node = (int) G->add_vertex((unsigned) color);
      assert( node == reflsymdata->ntreevars + v );

      G->add_edge(node, node - reflsymdata->ntreevars);

      ++nnodes;
   }
   nusedcolors = reflsymdata->nuniquevars;

   /* add nodes for rhs coefficients of constraints */
   for (int c = 0; c < reflsymdata->nuniquerhs; ++c)
   {
#ifndef NDEBUG
      int node = (int) G->add_vertex((unsigned) (nusedcolors + c));
      assert( node == 2 * reflsymdata->ntreevars + c );
#else
      (void) G->add_vertex((unsigned) (nusedcolors + c));
#endif
   }
   assert( (int) G->get_nof_vertices() == 2 * reflsymdata->ntreevars + reflsymdata->nuniquerhs );
   nusedcolors += reflsymdata->nuniquerhs;

   /* starting positions of different types of nodes */
   int rhsstart = 2 * reflsymdata->ntreevars;

   /* starting positions of colors of different type of nodes */
   int rhscolstart = reflsymdata->nuniquevars;
   int opcolstart = rhscolstart + reflsymdata->nuniquerhs;
   int coefcolstart = opcolstart + reflsymdata->nuniqueops;
   int valcolstart = coefcolstart + reflsymdata->nuniquecoefs;

   /* add nodes and edges for constraints */
   for (int c = 0; c < reflsymdata->ntreerhs; ++c)
   {
      /* keep track of path to root node, possible since tree is stored in DFS fashion */
      std::vector<unsigned> pathtoroot;
      std::vector<unsigned> pathtorootidx;
      std::vector<int> fliptypes;

      /* iterate through tree and create nodes and edges for constraint */
      for (int i = reflsymdata->treebegins[c]; i < reflsymdata->treebegins[c + 1]; ++i)
      {
         /* remove elements from the path that are not predecessors of the current node */
         while ( pathtoroot.size() > 0 && (int) pathtorootidx.back() > reflsymdata->treeparentidx[i] )
         {
            pathtorootidx.pop_back();
            pathtoroot.pop_back();
            fliptypes.pop_back();
         }

         if ( reflsymdata->treeparentidx[i] == -1 )
         {
            /* add root node and connect it with corresponding rhs node */
            assert( reflsymdata->trees[i] == SYM_NODETYPE_OPERATOR );

            unsigned node = G->add_vertex((unsigned) opcolstart + colorInClass(i, reflsymdata, FALSE));
            G->add_edge(node, (unsigned) (rhsstart + reflsymdata->rhscolors[c]));

            pathtoroot.push_back(node);
            pathtorootidx.push_back(i);
            fliptypes.push_back(getOperatorFliptype(scip, reflsymdata, i));
         }
         else if ( reflsymdata->trees[i] == SYM_NODETYPE_OPERATOR )
         {
            /* create node for the operator and connect it with its parent node */
            unsigned node = G->add_vertex((unsigned) opcolstart + colorInClass(i, reflsymdata, FALSE));
            G->add_edge(node, pathtoroot.back());

            SYM_FLIPTYPE fliptype = getOperatorFliptype(scip, reflsymdata, i);
            if ( fliptype != SYM_FLIPTYPE_BIPROD )
            {
               pathtoroot.push_back(node);
               pathtorootidx.push_back(i);
               fliptypes.push_back(fliptype);
            }
            else
            {
               /* already create nodes/edges for values, coefficients, and variables */
               assert( i + 1 < reflsymdata->treebegins[c + 1] - 1 );
               ++i;
               if ( reflsymdata->trees[i] == SYM_NODETYPE_VAL )
               {
                  unsigned valnode = G->add_vertex((unsigned) valcolstart + colorInClass(i, reflsymdata, FALSE));
                  G->add_edge(valnode, node);
                  ++i;
               }

               assert( i < reflsymdata->treebegins[c + 1] - 3 );
               assert( reflsymdata->trees[i] == SYM_NODETYPE_COEF );
               assert( reflsymdata->trees[i + 1] == SYM_NODETYPE_VAR );
               assert( reflsymdata->trees[i + 2] == SYM_NODETYPE_COEF );
               assert( reflsymdata->trees[i + 3] == SYM_NODETYPE_VAR );
               assert( colorInClass(i, reflsymdata, FALSE) == colorInClass(i + 2, reflsymdata, FALSE) );

               unsigned coef1node = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, FALSE));
               unsigned var1node = reflsymdata->treevaridx[reflsymdata->treemap[i + 1]];
               unsigned invvar1node = var1node + reflsymdata->ntreevars;
               unsigned coef2node = G->add_vertex((unsigned) coefcolstart + colorInClass(i + 2, reflsymdata, FALSE));
               unsigned var2node = reflsymdata->treevaridx[reflsymdata->treemap[i + 3]];
               unsigned invvar2node = var2node + reflsymdata->ntreevars;
               G->add_edge(coef1node, node);
               G->add_edge(coef2node, node);
               G->add_edge(var1node, coef1node);
               G->add_edge(var2node, coef1node);
               G->add_edge(invvar1node, coef2node);
               G->add_edge(invvar2node, coef2node);
               i += 3;
            }
         }
         else if ( reflsymdata->trees[i] == SYM_NODETYPE_COEF )
         {
            assert( i < reflsymdata->treebegins[c + 1] - 1 );
            assert( reflsymdata->trees[i + 1] == SYM_NODETYPE_VAR );

            /* create node for coefficient (and its negation if flips are allowed) */
            unsigned coefnode = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, FALSE));
            unsigned varnode = reflsymdata->treevaridx[reflsymdata->treemap[i + 1]];

            G->add_edge(coefnode, pathtoroot.back());
            G->add_edge(coefnode, varnode);

            switch ( fliptypes.back() )
            {
            case SYM_FLIPTYPE_EVEN :
            {
               unsigned invvarnode = varnode + reflsymdata->ntreevars;
               unsigned invcoefnode = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, TRUE));
               G->add_edge(coefnode, invcoefnode);
               G->add_edge(invcoefnode, pathtoroot.back());
               G->add_edge(invcoefnode, invvarnode);
               G->add_edge(coefnode, invvarnode);
               G->add_edge(invcoefnode, varnode);
               break;
            }
            case SYM_FLIPTYPE_SHIFT_ODD :
            case SYM_FLIPTYPE_ODD :
            {
               unsigned invcoefnode = G->add_vertex((unsigned) coefcolstart + colorInClass(i, reflsymdata, TRUE));
               unsigned invvarnode = varnode + reflsymdata->ntreevars;
               G->add_edge(coefnode, invcoefnode);
               G->add_edge(invcoefnode, pathtoroot.back());
               G->add_edge(invcoefnode, invvarnode);
               break;
            }
            default :
               assert ( fliptypes.back() == SYM_FLIPTYPE_NONE );
            }

            /* we have already taken the variable into account */
            ++i;
         }
         else if ( reflsymdata->trees[i] == SYM_NODETYPE_VAR )
         {
            SCIPinfoMessage(scip, NULL, "ERROR: terminate construction of symmetry detection graph, each variable should have a coefficient.\n");

            success = FALSE;
         }
         else
         {
            assert( reflsymdata->trees[i] == SYM_NODETYPE_VAL );

            /* create node for value and connect it with its parent node */
            unsigned node = G->add_vertex((unsigned) valcolstart + colorInClass(i, reflsymdata, FALSE));

            G->add_edge(node, pathtoroot.back());
         }
      }
   }

   return SCIP_OKAY;
}

/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}

char*
initStaticBlissName( );

static char* blissname = initStaticBlissName();

char*
initStaticBlissName( )
{
   blissname = new char[100];
#ifdef BLISS_PATCH_PRESENT
   (void) snprintf(blissname, 100, "bliss %sp", bliss::version);
#else
   (void) snprintf(blissname, 100, "bliss %s", bliss::version);
#endif
   return blissname;
}


/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
   return blissname;
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "Computing Graph Automorphism Groups by T. Junttila and P. Kaski (www.tcs.hut.fi/Software/bliss/)";
}

/** compute generators of symmetry group */
SCIP_RETCODE SYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   SYM_EXPRDATA*         exprdata,           /**< data for nonlinear constraints */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize      /**< pointer to store size of group */
   )
{
   assert( scip != NULL );
   assert( matrixdata != NULL );
   assert( exprdata != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( maxgenerators >= 0 );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;

   int nnodes = 0;
   int nedges = 0;
   int nusedcolors = 0;
   SCIP_Bool success = FALSE;

   /* create bliss graph */
   bliss::Graph G(0);

   /* create nodes corresponding to variables */
   SCIP_CALL( createVariableNodes(scip, &G, matrixdata, nnodes, nedges, nusedcolors) );

   assert( nnodes == matrixdata->npermvars );
   assert( nusedcolors == matrixdata->nuniquevars );

   /* fill graph with nodes for variables and linear constraints */
   SCIP_CALL( fillGraphByLinearConss(scip, &G, matrixdata, nnodes, nedges, nusedcolors, success) );

   if ( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Graph construction failed during linear part.\n");
      return SCIP_OKAY;
   }

   /* add the nodes for nonlinear constraints to the graph */
   SCIP_CALL( fillGraphByNonlinearConss(scip, &G, exprdata, nnodes, nedges, nusedcolors, success) );

   if ( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Graph construction failed during non-linear part.\n");
      return SCIP_OKAY;
   }

#ifdef SCIP_OUTPUT
   G.write_dot("debug.dot");
#endif

   SCIPdebugMsg(scip, "Symmetry detection graph has %u nodes.\n", G.get_nof_vertices());

   /* compute automorphisms */
   bliss::Stats stats;
   BLISS_Data data;
   data.scip = scip;
   data.npermvars = matrixdata->npermvars;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   G.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   G.set_component_recursion(false);

   /* do not use a node limit, but set generator limit */
#ifdef BLISS_PATCH_PRESENT
   G.set_search_limits(0, (unsigned) maxgenerators);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and pass it to the blisshook above */
   auto reportglue = [&](unsigned int n, const unsigned int* aut) {
      blisshook((void*)&data, n, aut);
   };

   /* lambda function to have access to stats and terminate the search if maxgenerators are reached */
   auto term = [&]() {
      return (stats.get_nof_generators() >= (long unsigned int) maxgenerators);
   };

   /* start search */
   G.find_automorphisms(stats, reportglue, term);
#else
   /* start search */
   G.find_automorphisms(stats, blisshook, (void*) &data);
#endif


#ifdef SCIP_OUTPUT
   (void) stats.print(stdout);
#endif

   /* prepare return values */
   if ( data.nperms > 0 )
   {
      *perms = data.perms;
      *nperms = data.nperms;
      *nmaxperms = data.nmaxperms;
   }
   else
   {
      assert( data.perms == NULL );
      assert( data.nmaxperms == 0 );
   }

   /* determine log10 of symmetry group size */
   *log10groupsize = (SCIP_Real) log10l(stats.get_group_size_approx());

   return SCIP_OKAY;
}

/** compute generators of reflection symmetry group */
SCIP_RETCODE SYMcomputeReflectionSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_REFLSYMDATA*      reflsymdata,        /**< data for computing reflection symmetries */
   int                   maxgenerators,      /**< maximum number of generators to be computed */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize      /**< pointer to store size of group */
   )
{
   SCIP_Bool success;

   assert( scip != NULL );
   assert( reflsymdata != NULL );

   /* create bliss graph */
   bliss::Graph G(0);

   /* create nodes corresponding to variables */
   SCIP_CALL( SCIPcreateReflectionSymmetryDetectionGraph(scip, &G, reflsymdata, success) );

   if ( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Graph construction failed.\n");
      return SCIP_OKAY;
   }

#ifdef SCIP_OUTPUT
   G.write_dot("debug.dot");
#endif

   SCIPdebugMsg(scip, "Symmetry detection graph has %u nodes.\n", G.get_nof_vertices());

   /* compute automorphisms */
   bliss::Stats stats;
   BLISS_Data data;
   data.scip = scip;
   data.npermvars = reflsymdata->ntreevars;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   G.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   G.set_component_recursion(false);

   /* do not use a node limit, but set generator limit */
#ifdef BLISS_PATCH_PRESENT
   G.set_search_limits(0, (unsigned) maxgenerators);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to data and pass it to the blisshook above */
   auto reportglue = [&](unsigned int n, const unsigned int* aut) {
      blisshookReflSym((void*)&data, n, aut);
   };

   /* lambda function to have access to stats and terminate the search if maxgenerators are reached */
   auto term = [&]() {
      return (stats.get_nof_generators() >= (long unsigned int) maxgenerators);
   };

   /* start search */
   G.find_automorphisms(stats, reportglue, term);
#else
   /* start search */
   G.find_automorphisms(stats, blisshook, (void*) &data);
#endif

#ifdef SCIP_OUTPUT
   (void) stats.print(stdout);
#endif

   /* prepare return values */
   if ( data.nperms > 0 )
   {
      *perms = data.perms;
      *nperms = data.nperms;
      *nmaxperms = data.nmaxperms;
   }
   else
   {
      assert( data.perms == NULL );
      assert( data.nmaxperms == 0 );
   }

   /* determine log10 of symmetry group size */
   *log10groupsize = (SCIP_Real) log10l(stats.get_group_size_approx());


   printf("Found %d generators, log10 group size %.1f\n", data.nperms, (SCIP_Real) log10l(stats.get_group_size_approx()));
   for (int p = 0; p < data.nperms; ++p)
   {
      SCIP_Bool hasreflection = FALSE;

      for (int v = 0; v < data.npermvars; ++v)
      {
         if ( data.perms[p][v] >= data.npermvars )
         {
            hasreflection = TRUE;
            break;
         }
      }

      if ( hasreflection )
      {
         // print permutation in cycle representation
         SCIP_Bool* covered;

         SCIP_CALL( SCIPallocClearBufferArray(scip, &covered, data.npermvars) );

         printf("perm %d: ", p);

         for (int v = 0; v < data.npermvars; ++v)
         {
            if ( covered[v] )
               continue;

            if ( data.perms[p][v] == v )
               continue;

            printf("(%d", v);
            covered[v] = TRUE;

            int w = data.perms[p][v];
            while ( w != v )
            {
               if ( w >= data.npermvars )
               {
                  printf(" -%d", w);
                  covered[w - data.npermvars] = TRUE;
               }
               else
               {
                  printf(" %d", w);
                  covered[w] = TRUE;
               }

               w = data.perms[p][w];
            }
            printf(")");
         }

         printf("\n");
         for (int v = 0; v < data.npermvars; ++v)
            covered[v] = FALSE;

         for (int v = 0; v < data.npermvars; ++v)
         {
            if ( covered[v] )
               continue;

            if ( data.perms[p][v] == v )
               continue;

            printf("(%s", SCIPvarGetName(SCIPgetVars(scip)[v]));
            covered[v] = TRUE;

            int w = data.perms[p][v];
            while ( w != v )
            {
               if ( w >= data.npermvars )
               {
                  printf(" -%s", SCIPvarGetName(SCIPgetVars(scip)[w - data.npermvars]));
                  covered[w - data.npermvars] = TRUE;
               }
               else
               {
                  printf(" %s", SCIPvarGetName(SCIPgetVars(scip)[w]));
                  covered[w] = TRUE;
               }

               w = data.perms[p][w];
            }
            printf(")");
         }

         SCIPfreeBufferArray(scip, &covered);
         printf("\n\n");
      }
   }

   return SCIP_OKAY;
}
