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

/**@file   compute_symmetry_sassy.c
 * @brief  interface for symmetry computations to sassy as a preprocessor to bliss
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/* include bliss */
#include <bliss/defs.hh>
#include <bliss/graph.hh>

/* include sassy */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4189)  // local variable is initialized but not referenced
# pragma warning(disable: 4267)  // conversion of size_t to int (at sassy/preprocessor.h:2897)
# pragma warning(disable: 4388)  // compare signed and unsigned expression
# pragma warning(disable: 4456)  // shadowed variable
# pragma warning(disable: 4430)  // missing type specifier
#endif

/* the actual include */
#include <sassy/preprocessor.h>

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wunused-but-set-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wshadow"
#endif

#ifdef _MSC_VER
# pragma warning(pop)
#endif

#include <sassy/tools/bliss_converter.h>

#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_pow.h"
#include "scip/expr.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"
#include "scip/scip_mem.h"


/** struct for symmetry callback */
struct SYMMETRY_Data
{
   SCIP*                 scip;               /**< SCIP pointer */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nperms;             /**< number of permutations */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
   int                   nmaxperms;          /**< maximal number of permutations */
   int                   maxgenerators;      /**< maximal number of generators to be constructed (= 0 if unlimited) */
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
   if ( SCIPisExprPower((SCIP*) userptr, k1->expr )
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
   uint64_t result;

   k = (SYM_OPTYPE*) key;

   if ( SCIPisExprPower((SCIP*) userptr, k->expr) )
      exponent = SCIPgetExponentExprPow(k->expr);
   else
      exponent = 1.0;

   result = SCIPhashThree(SCIPrealHashCode(exponent), k->level,
      SCIPhashKeyValString(NULL, static_cast<void*>(const_cast<char*>(SCIPexprhdlrGetName(SCIPexprGetHdlr(k->expr))))));

   return result;
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


/* ------------------- hook functions ------------------- */

/** callback function for sassy */  /*lint -e{715}*/
static
void sassyhook(
   void*                 user_param,         /**< parameter supplied at call to sassy */
   int                   n,                  /**< dimension of permutations */
   const int*            aut,                /**< permutation */
   int                   nsupp,              /**< support size */
   const int*            suppa               /**< support list */
   )
{
   assert( aut != NULL );
   assert( user_param != NULL );

   SYMMETRY_Data* data = static_cast<SYMMETRY_Data*>(user_param);
   assert( data->scip != NULL );
   assert( data->npermvars < (int) n );
   assert( data->maxgenerators >= 0 );

   /* make sure we do not generate more that maxgenerators many permutations */
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


/* ------------------- other functions ------------------- */

/** determine number of nodes and edges */
static
SCIP_RETCODE determineGraphSize(
   SCIP*                 scip,               /**< SCIP instance */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   SYM_EXPRDATA*         exprdata,           /**< data for nonlinear constraints */
   int*                  nnodes,             /**< pointer to store the total number of nodes in graph */
   int*                  nedges,             /**< pointer to store the total number of edges in graph */
   int*                  nlinearnodes,       /**< pointer to store the number of internal nodes for linear constraints */
   int*                  nnonlinearnodes,    /**< pointer to store the number of internal nodes for nonlinear constraints */
   int*                  nlinearedges,       /**< pointer to store the number of edges for linear constraints */
   int*                  nnonlinearedges,    /**< pointer to store the number of edges for nonlinear constraints */
   int**                 degrees,            /**< pointer to store the degrees of the nodes */
   int*                  maxdegrees,         /**< pointer to store the maximal size of the degree array */
   SCIP_Bool*            success             /**< pointer to store whether the construction was successful */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool groupByConstraints;
   int* internodes = NULL;
   int nmaxinternodes;
   int oldcolor = -1;
#ifndef NDEBUG
   SCIP_Real oldcoef = SCIP_INVALID;
#endif
   int firstcolornodenumber = -1;
   int nconss;
   int j;

   assert( scip != NULL );
   assert( matrixdata != NULL );
   assert( exprdata != NULL );
   assert( nnodes != NULL );
   assert( nedges != NULL );
   assert( nlinearnodes != NULL );
   assert( nnonlinearnodes != NULL );
   assert( nlinearedges != NULL );
   assert( nnonlinearedges != NULL );
   assert( degrees != NULL );
   assert( maxdegrees != NULL );
   assert( success != NULL );

   *success = TRUE;

   /* count nodes for variables */
   *nnodes = matrixdata->npermvars;

   /* count nodes for rhs of constraints */
   *nnodes += matrixdata->nrhscoef;

   /* allocate memory for degrees (will grow dynamically) */
   *degrees = NULL;
   *maxdegrees = 0;
   SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes + 100) );
   for (j = 0; j < *nnodes; ++j)
      (*degrees)[j] = 0;

   /* initialize counters */
   *nedges = 0;
   *nlinearnodes = 0;
   *nnonlinearnodes = 0;
   *nlinearedges = 0;
   *nnonlinearedges = 0;

   /* Determine grouping depending on the number of rhs vs. variables; see fillGraphByLinearConss(). */
   if ( matrixdata->nrhscoef < matrixdata->npermvars )
      groupByConstraints = TRUE;
   else
      groupByConstraints = FALSE;

   /* determine size of intermediate nodes */
   if ( groupByConstraints )
      nmaxinternodes = matrixdata->nrhscoef;
   else
      nmaxinternodes = matrixdata->npermvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &internodes, nmaxinternodes) ); /*lint !e530*/
   for (j = 0; j < nmaxinternodes; ++j)
      internodes[j] = -1;

   /* loop through all matrix coefficients */
   for (j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int varrhsidx;
      int rhsnode;
      int varnode;
      int color;
      int idx;

      idx = matrixdata->matidx[j];
      assert( 0 <= idx && idx < matrixdata->nmatcoef );

      /* find color corresponding to matrix coefficient */
      color = matrixdata->matcoefcolors[idx];
      assert( 0 <= color && color < matrixdata->nuniquemat );

      assert( 0 <= matrixdata->matrhsidx[idx] && matrixdata->matrhsidx[idx] < matrixdata->nrhscoef );
      assert( 0 <= matrixdata->matvaridx[idx] && matrixdata->matvaridx[idx] < matrixdata->npermvars );

      rhsnode = matrixdata->npermvars + matrixdata->matrhsidx[idx];
      varnode = matrixdata->matvaridx[idx];
      assert( matrixdata->npermvars <= rhsnode && rhsnode < matrixdata->npermvars + matrixdata->nrhscoef );
      assert( rhsnode < *nnodes );
      assert( varnode < *nnodes );

      if ( matrixdata->nuniquemat == 1 )
      {
         /* We do not need intermediate nodes if we have only one coefficient class; just add edges. */
         ++(*degrees)[varnode];
         ++(*degrees)[rhsnode];
         ++(*nlinearedges);
      }
      else
      {
         SCIP_Bool newinternode = FALSE;
         int internode;

         /* if new group of coefficients has been reached */
         if ( color != oldcolor )
         {
            assert( ! SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );
            oldcolor = color;
            firstcolornodenumber = *nnodes;
#ifndef NDEBUG
            oldcoef = matrixdata->matcoef[idx];
#endif
         }
         else
            assert( SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );

         if ( groupByConstraints )
            varrhsidx = matrixdata->matrhsidx[idx];
         else
            varrhsidx = matrixdata->matvaridx[idx];
         assert( 0 <= varrhsidx && varrhsidx < nmaxinternodes );

         if ( internodes[varrhsidx] < firstcolornodenumber )
         {
            internodes[varrhsidx] = (*nnodes)++;
            ++(*nlinearnodes);

            /* ensure memory for degrees */
            SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes) );
            (*degrees)[internodes[varrhsidx]] = 0;
            newinternode = TRUE;
         }
         internode = internodes[varrhsidx];
         assert( internode >= matrixdata->npermvars + matrixdata->nrhscoef );
         assert( internode >= firstcolornodenumber );

         /* determine whether graph would be too large */
         if ( *nnodes >= INT_MAX/2 )
         {
            *success = FALSE;
            break;
         }

         if ( groupByConstraints )
         {
            if ( newinternode )
            {
               ++(*degrees)[rhsnode];
               ++(*degrees)[internode];
               ++(*nlinearedges);
            }
            ++(*degrees)[varnode];
            ++(*degrees)[internode];
            ++(*nlinearedges);
         }
         else
         {
            if ( newinternode )
            {
               ++(*degrees)[varnode];
               ++(*degrees)[internode];
               ++(*nlinearedges);
            }
            ++(*degrees)[rhsnode];
            ++(*degrees)[internode];
            ++(*nlinearedges);
         }
      }
   }
   SCIPfreeBufferArray(scip, &internodes);

   /* now treat nonlinear constraints */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   nconss = conshdlr != NULL ? SCIPconshdlrGetNConss(conshdlr) : 0;
   if ( nconss > 0 && *success )
   {
      SCIP_CONS** conss;
      SCIP_EXPRITER* it;
      SCIP_VAR** vars = NULL;
      SCIP_Real* vals = NULL;
      int* visitednodes = NULL;
      int* ischildofsum = NULL;
      int maxvisitednodes;
      int maxischildofsum;
      int numvisitednodes = 0;
      int numischildofsum = 0;
      int varssize;
      int i;

      conss = SCIPconshdlrGetConss(conshdlr);

      /* prepare iterator */
      SCIP_CALL( SCIPcreateExpriter(scip, &it) );

      /* prepare stacks */
      maxvisitednodes = exprdata->nuniqueoperators + exprdata->nuniqueconstants + exprdata->nuniquecoefs;
      maxischildofsum = maxvisitednodes;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &visitednodes, maxvisitednodes) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ischildofsum, maxischildofsum) );

      /* get number of variables */
      varssize = SCIPgetNVars(scip);

      /* iterate over all expressions and add the corresponding nodes to the graph */
      for (i = 0; i < nconss; ++i)
      {
         SCIP_EXPR* rootexpr;
         SCIP_EXPR* expr;
#ifndef NDEBUG
         int currentlevel = 0;
#endif

         rootexpr = SCIPgetExprNonlinear(conss[i]);

         SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
         SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_LEAVEEXPR);

         for (expr = SCIPexpriterGetCurrent(it); ! SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it)) /*lint !e441*/ /*lint !e440*/
         {
            /* upon entering an expression, check its type and add nodes and edges if neccessary */
            switch ( SCIPexpriterGetStageDFS(it) )
            {
            case SCIP_EXPRITER_ENTEREXPR:
            {
               int node = -1;
               int parentnode = -1;
               SCIP_Bool isVarExpr = FALSE;

               /* for variable expressions, get the corresponding node that is already in the graph */
               if ( SCIPisExprVar(scip, expr) )
               {
                  SCIP_VAR* var;

                  var = SCIPgetVarExprVar(expr);
                  isVarExpr = TRUE;

                  /* Check whether the variable is active; if not, then replace the inactive variable by its aggregation
                   * or its fixed value; note that this step is equivalent as representing an inactive variable as sum
                   * expression.
                   */
                  if ( SCIPvarIsActive(var) )
                  {
                     node = SCIPvarGetProbindex(var);
                     assert( node < *nnodes );
                  }
                  else
                  {
                     SCIP_Real constant = 0.0;
                     int nvars;
                     int requiredsize;
                     int k;

                     if ( vars == NULL )
                     {
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, varssize) );
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vals, varssize) );
                     }
                     assert( vars != NULL && vals != NULL );

                     vars[0] = var;
                     vals[0] = 1.0;
                     nvars = 1;

                     SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, vals, &nvars, varssize, &constant, &requiredsize, TRUE) );
                     assert( requiredsize <= varssize );

                     assert( numvisitednodes > 0 );
                     parentnode = visitednodes[numvisitednodes-1];
                     assert( parentnode < *nnodes );

                     /* create nodes for all aggregation variables and coefficients and connect them to the parent node */
                     for (k = 0; k < nvars; ++k)
                     {
                        int internode;

                        assert( vars[k] != NULL );
                        assert( vals[k] != 0.0 );

                        /* add node */
                        internode = (*nnodes)++;
                        ++(*nnonlinearnodes);

                        /* ensure size of degrees */
                        SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes) );
                        (*degrees)[internode] = 0;

                        ++(*degrees)[internode];
                        ++(*degrees)[parentnode];
                        ++(*nnonlinearedges);

                        /* connect the intermediate node to its corresponding variable node */
                        node = SCIPvarGetProbindex(vars[k]);
                        assert( node < *nnodes );

                        ++(*degrees)[internode];
                        ++(*degrees)[node];
                        ++(*nnonlinearedges);
                     }

                     /* add the node for the constant */
                     if ( constant != 0.0 )
                     {
                        /* add node */
                        node = (*nnodes)++;
                        ++(*nnonlinearnodes);

                        /* ensure size of degrees */
                        SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes) );
                        (*degrees)[node] = 0;

                        ++(*degrees)[node];
                        ++(*degrees)[parentnode];
                        ++(*nnonlinearedges);
                     }

                     /* add a filler node since it will be removed in the next iteration anyway */
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &visitednodes, &maxvisitednodes, numvisitednodes+1) );
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ischildofsum, &maxischildofsum, numischildofsum+1) );

                     visitednodes[numvisitednodes++] = *nnodes;
                     ischildofsum[numischildofsum++] = FALSE;
#ifndef NDEBUG
                     ++currentlevel;
#endif

                     break;
                  }
               }
               /* for all other expressions, no nodes or edges have to be created */
               else
               {
                  /* do nothing here */
               }

               /* if this is the root expression, add the constraint side node (will be parent of expression node) */
               if ( SCIPexpriterGetParentDFS(it) == NULL )
               {
                  /* add the constraint side node */
                  parentnode = (*nnodes)++;
                  ++(*nnonlinearnodes);

                  /* ensure size of degrees */
                  SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes) );
                  (*degrees)[parentnode] = 0;
               }
               /* otherwise, get the parentnode stored in visitednodes */
               else
               {
                  parentnode = visitednodes[numvisitednodes - 1];
                  assert( parentnode < *nnodes );
               }

               /* in all cases apart from variable expressions, the new node is added with the corresponding color */
               if ( ! isVarExpr )
               {
                  node = (*nnodes)++;
                  ++(*nnonlinearnodes);

                  /* ensure size of degrees */
                  SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes) );
                  (*degrees)[node] = 0;
               }

               /* store the new node so that it can be used as parentnode later */
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, &visitednodes, &maxvisitednodes, numvisitednodes+1) );
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ischildofsum, &maxischildofsum, numischildofsum+1) );

               assert( node != -1 );
               visitednodes[numvisitednodes++] = node;
               ischildofsum[numischildofsum++] = FALSE;

               /* connect the current node with its parent */
               assert( parentnode != -1 );
               ++(*degrees)[node];
               ++(*degrees)[parentnode];
               ++(*nnonlinearedges);

               /* for sum expression, also add intermediate nodes for the coefficients */
               if ( SCIPisExprSum(scip, expr) )
               {
                  SCIP_Real constval;
                  int internode;

                  /* iterate over children from last to first, such that visitednodes array is in correct order */
                  for (j = SCIPexprGetNChildren(expr) - 1; j >= 0; --j)
                  {
                     /* add the intermediate node with the corresponding color */
                     internode = (*nnodes)++;
                     ++(*nnonlinearnodes);

                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &visitednodes, &maxvisitednodes, numvisitednodes+1) );
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ischildofsum, &maxischildofsum, numischildofsum+1) );

                     visitednodes[numvisitednodes++] = internode;
                     ischildofsum[numischildofsum++] = TRUE;

                     /* ensure size of degrees */
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes) );
                     (*degrees)[internode] = 0;

                     ++(*degrees)[internode];
                     ++(*degrees)[node];
                     ++(*nnonlinearedges);
                  }

                  /* add node for the constant term of the sum expression */
                  constval = SCIPgetConstantExprSum(expr);
                  if ( constval != 0.0 )
                  {
                     /* add the node with a new color */
                     internode = (*nnodes)++;
                     ++(*nnonlinearnodes);

                     /* ensure size of degrees */
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, degrees, maxdegrees, *nnodes) );
                     (*degrees)[internode] = 0;

                     ++(*degrees)[internode];
                     ++(*degrees)[node];
                     ++(*nnonlinearedges);
                  }
               }

#ifndef NDEBUG
               ++currentlevel;
#endif
               break;
            }
            /* when leaving an expression, the nodes that are not needed anymore are erased from the respective arrays */
            case SCIP_EXPRITER_LEAVEEXPR:
            {
               --numvisitednodes;
               --numischildofsum;
#ifndef NDEBUG
               currentlevel--;
#endif

               /* When leaving the child of a sum expression, we have to pop again to get rid of the intermediate nodes
                * used for the coefficients of summands
                */
               if ( numischildofsum > 0 && ischildofsum[numischildofsum - 1] )
               {
                  --numvisitednodes;
                  --numischildofsum;
               }

               break;
            }

            default:
               SCIPABORT(); /* we should never be called in this stage */
               break;
            }
         }

         assert( currentlevel == 0 );
         assert( numvisitednodes == 0 );
         assert( numischildofsum == 0 );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &vals, varssize);
      SCIPfreeBlockMemoryArrayNull(scip, &vars, varssize);

      SCIPfreeBlockMemoryArray(scip, &visitednodes, maxvisitednodes);
      SCIPfreeBlockMemoryArray(scip, &ischildofsum, maxischildofsum);
      SCIPfreeExpriter(&it);
   }

   *nedges = *nlinearedges + *nnonlinearedges;
   assert( *nnodes == matrixdata->npermvars + matrixdata->nrhscoef + *nlinearnodes + *nnonlinearnodes );

   SCIPdebugMsg(scip, "#nodes for variables: %d\n", matrixdata->npermvars);
   SCIPdebugMsg(scip, "#nodes for rhs: %d\n", matrixdata->nrhscoef);
   SCIPdebugMsg(scip, "#intermediate nodes for linear constraints: %d\n", *nlinearnodes);
   SCIPdebugMsg(scip, "#intermediate nodes for nonlinear constraints: %d\n", *nnonlinearnodes);
   SCIPdebugMsg(scip, "#edges for linear constraints: %d\n", *nlinearedges);
   SCIPdebugMsg(scip, "#edges for nonlinear constraints: %d\n", *nnonlinearedges);

   return SCIP_OKAY;
}


/** Construct linear and nonlinear part of colored graph for symmetry computations
 *
 *  Construct linear graph:
 *  - Each variable gets a different node.
 *  - Each constraint gets a different node.
 *  - Each matrix coefficient gets a different node that is connected to the two nodes
 *    corresponding to the respective constraint and variable.
 *  - Each different variable, rhs, matrix coefficient gets a different color that is attached to the corresponding entries.
 *
  *  Construct nonlinear graph:
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
 */
static
SCIP_RETCODE fillGraphByConss(
   SCIP*                 scip,               /**< SCIP instance */
   sassy::static_graph*  G,                  /**< graph to be constructed */
   SYM_MATRIXDATA*       matrixdata,         /**< data for MIP matrix */
   SYM_EXPRDATA*         exprdata,           /**< data for nonlinear constraints */
   int                   nnodes,             /**< total number of nodes in graph */
   int                   nedges,             /**< total number of edges in graph */
   int                   nlinearnodes,       /**< number of intermediate nodes for linear constraints */
   int                   nnonlinearnodes,    /**< number of intermediate nodes for nonlinear constraints */
   int                   nlinearedges,       /**< number of intermediate edges for linear constraints */
   int                   nnonlinearedges,    /**< number of intermediate edges for nonlinear constraints */
   int*                  degrees,            /**< array with the degrees of the nodes */
   int*                  nusedcolors         /**< pointer to store number of used colors in the graph so far */
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
   SCIP_EXPRITER* it;
   SCIP_Bool* ischildofsum = NULL;
   SCIP_VAR** vars = NULL;
   SCIP_Real* vals = NULL;
   SCIP_Bool groupByConstraints;
   int* internodes = NULL;
   int* visitednodes = NULL;
   int maxischildofsum;
   int maxvisitednodes;
   int numvisitednodes = 0;
   int numischildofsum = 0;
   int nconss;
   int nuniqueops = 0;
   int nuniqueconsts = 0;
   int nuniquecoefs = 0;
   int nuniquerhs = 0;
   int oparraysize;
   int constarraysize;
   int coefarraysize;
   int rhsarraysize;
   int nmaxinternodes;
   int oldcolor = -1;
   int varssize;
#ifndef NDEBUG
   SCIP_Real oldcoef = SCIP_INVALID;
   int m = 0;
#endif
   int firstcolornodenumber = -1;
   int n = 0;
   int i;
   int j;

   assert( scip != NULL );
   assert( G != NULL );
   assert( matrixdata != NULL );
   assert( exprdata != NULL );
   assert( degrees != NULL );
   assert( nusedcolors != NULL );

   SCIPdebugMsg(scip, "Filling graph with colored coefficient nodes for linear part.\n");

   /* fill in array with colors for variables */
   for (j = 0; j < matrixdata->npermvars; ++j)
   {
      assert( 0 <= matrixdata->permvarcolors[j] && matrixdata->permvarcolors[j] < matrixdata->nuniquevars );
      (void) G->add_vertex(matrixdata->permvarcolors[j], degrees[n]);
      ++n;
   }
   *nusedcolors = matrixdata->nuniquevars;

   /* fill in array with colors for rhs */
   for (i = 0; i < matrixdata->nrhscoef; ++i)
   {
      assert( 0 <= matrixdata->rhscoefcolors[i] && matrixdata->rhscoefcolors[i] < matrixdata->nuniquerhs );
      (void) G->add_vertex(*nusedcolors + matrixdata->rhscoefcolors[i], degrees[n]);
      ++n;
   }
   *nusedcolors += matrixdata->nuniquerhs;

   /* Grouping of nodes depends on the number of nodes in the bipartite graph class.  If there are more variables than
    * constraints, we group by constraints.  That is, given several variable nodes which are incident to one constraint
    * node by the same color, we join these variable nodes to the constraint node by only one intermediate node.
    */
   if ( matrixdata->nrhscoef < matrixdata->npermvars )
      groupByConstraints = TRUE;
   else
      groupByConstraints = FALSE;

   /* "colored" edges based on all matrix coefficients - loop through ordered matrix coefficients */
   if ( groupByConstraints )
      nmaxinternodes = matrixdata->nrhscoef;
   else
      nmaxinternodes = matrixdata->npermvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &internodes, nmaxinternodes) ); /*lint !e530*/
   for (j = 0; j < nmaxinternodes; ++j)
      internodes[j] = -1;

   /* We pass through the matrix coeficients, grouped by color, i.e., different coefficients. If the coeffients appear
    * in the same row or column, it suffices to only generate a single node (depending on groupByConstraints). We store
    * this node in the array internodes. In order to avoid reinitialization, we store the node number with increasing
    * numbers for each color. The smallest number for the current color is stored in firstcolornodenumber. */
   for (j = 0; j < matrixdata->nmatcoef; ++j)
   {
      int idx;
      int color;
      int rhsnode;
      int varnode;
      int varrhsidx;

      idx = matrixdata->matidx[j];
      assert( 0 <= idx && idx < matrixdata->nmatcoef );

      /* find color corresponding to matrix coefficient */
      color = matrixdata->matcoefcolors[idx];
      assert( 0 <= color && color < matrixdata->nuniquemat );

      assert( 0 <= matrixdata->matrhsidx[idx] && matrixdata->matrhsidx[idx] < matrixdata->nrhscoef );
      assert( 0 <= matrixdata->matvaridx[idx] && matrixdata->matvaridx[idx] < matrixdata->npermvars );

      rhsnode = matrixdata->npermvars + matrixdata->matrhsidx[idx];
      varnode = matrixdata->matvaridx[idx];
      assert( matrixdata->npermvars <= rhsnode && rhsnode < matrixdata->npermvars + matrixdata->nrhscoef );
      assert( rhsnode < nnodes );
      assert( varnode < nnodes );

      /* if we have only one color, we do not need intermediate nodes */
      if ( matrixdata->nuniquemat == 1 )
      {
         if ( rhsnode < varnode )
            G->add_edge((unsigned) rhsnode, (unsigned) varnode);
         else
            G->add_edge((unsigned) varnode, (unsigned) rhsnode);
#ifndef NDEBUG
         ++m;
#endif
      }
      else
      {
         SCIP_Bool newinternode = FALSE;
         int internode;

         /* if new group of coefficients has been reached */
         if ( color != oldcolor )
         {
            assert( ! SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );
            oldcolor = color;
            firstcolornodenumber = n;
#ifndef NDEBUG
            oldcoef = matrixdata->matcoef[idx];
#endif
         }
         else
            assert( SCIPisEQ(scip, oldcoef, matrixdata->matcoef[idx]) );

         if ( groupByConstraints )
            varrhsidx = matrixdata->matrhsidx[idx];
         else
            varrhsidx = matrixdata->matvaridx[idx];
         assert( 0 <= varrhsidx && varrhsidx < nmaxinternodes );

         if ( internodes[varrhsidx] < firstcolornodenumber )
         {
            (void) G->add_vertex(*nusedcolors + color, degrees[n]);
            internodes[varrhsidx] = n++;
            newinternode = TRUE;
         }
         internode = internodes[varrhsidx];
         assert( internode >= matrixdata->npermvars + matrixdata->nrhscoef );
         assert( internode >= firstcolornodenumber );
         assert( internode < nnodes );

         if ( groupByConstraints )
         {
            if ( newinternode )
            {
               G->add_edge((unsigned) rhsnode, (unsigned) internode);
#ifndef NDEBUG
               ++m;
#endif
            }
            G->add_edge((unsigned) varnode, (unsigned) internode);
#ifndef NDEBUG
            ++m;
#endif
         }
         else
         {
            if ( newinternode )
            {
               G->add_edge((unsigned) varnode, (unsigned) internode);
#ifndef NDEBUG
               ++m;
#endif
            }
            G->add_edge((unsigned) rhsnode, (unsigned) internode);
#ifndef NDEBUG
            ++m;
#endif
         }
      }
   }
   assert( n == matrixdata->npermvars + matrixdata->nrhscoef + nlinearnodes );
   assert( m == nlinearedges );

   SCIPfreeBufferArray(scip, &internodes);

   *nusedcolors += matrixdata->nuniquemat;

   /* ------------------------------------------------------------------------ */
   /* treat nonlinear constraints */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   nconss = conshdlr != NULL ? SCIPconshdlrGetNConss(conshdlr) : 0;
   if ( nconss == 0 )
   {
      return SCIP_OKAY;
   }

   conss = SCIPconshdlrGetConss(conshdlr);
   rhsarraysize = nconss;

   SCIPdebugMsg(scip, "Filling graph with colored coefficient nodes for non-linear part.\n");

   /* create maps for optypes, constants, sum coefficients and rhs to indices */
   oparraysize = exprdata->nuniqueoperators;
   constarraysize = exprdata->nuniqueconstants;
   coefarraysize = exprdata->nuniquecoefs;

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

   SCIP_CALL( SCIPcreateExpriter(scip, &it) );

   maxvisitednodes = oparraysize + constarraysize + coefarraysize;
   maxischildofsum = maxvisitednodes;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &visitednodes, maxvisitednodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ischildofsum, maxischildofsum) );

   /* get number of variables */
   varssize = SCIPgetNVars(scip);

   /* iterate over all expressions and add the corresponding nodes to the graph */
   for (i = 0; i < nconss; ++i)
   {
      SCIP_EXPR* rootexpr;
      SCIP_EXPR* expr;
      int currentlevel = 0;

      rootexpr = SCIPgetExprNonlinear(conss[i]);

      SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );
      SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_ENTEREXPR | SCIP_EXPRITER_LEAVEEXPR);

      for (expr = SCIPexpriterGetCurrent(it); ! SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it)) /*lint !e441*/ /*lint !e440*/
      {
         /* upon entering an expression, check its type and add nodes and edges if neccessary */
         switch ( SCIPexpriterGetStageDFS(it) )
         {
            case SCIP_EXPRITER_ENTEREXPR:
            {
               int node = -1;
               int parentnode = -1;
               int color = -1;

               /* for variable expressions, get the corresponding node that is already in the graph */
               if ( SCIPisExprVar(scip, expr) )
               {
                  SCIP_VAR* var;

                  var = SCIPgetVarExprVar(expr);

                  /* Check whether the variable is active; if not, then replace the inactive variable by its aggregation
                   * or its fixed value; note that this step is equivalent as representing an inactive variable as sum
                   * expression.
                   */
                  if ( SCIPvarIsActive(var) )
                  {
                     node = SCIPvarGetProbindex(var);
                     assert( node < nnodes );
                  }
                  else
                  {
                     SCIP_Real constant = 0.0;
                     int nvars;
                     int requiredsize;
                     int k;

                     if ( vars == NULL )
                     {
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, varssize) );
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vals, varssize) );
                     }
                     assert( vars != NULL && vals != NULL );

                     vars[0] = var;
                     vals[0] = 1.0;
                     nvars = 1;

                     SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, vals, &nvars, varssize, &constant, &requiredsize, TRUE) );
                     assert( requiredsize <= varssize );

                     assert( numvisitednodes > 0 );

                     parentnode = visitednodes[numvisitednodes-1];
                     assert( parentnode < nnodes );

                     /* create nodes for all aggregation variables and coefficients and connect them to the parent node */
                     for (k = 0; k < nvars; ++k)
                     {
                        SYM_CONSTTYPE* ct;
                        int internode;

                        assert( vars[k] != NULL );
                        assert( vals[k] != 0.0 );
                        assert( nuniquecoefs < coefarraysize );

                        ct = &sumcoefarray[nuniquecoefs];
                        ct->value = vals[k];

                        if ( ! SCIPhashtableExists(sumcoefmap, (void *) ct) )
                        {
                           SCIP_CALL( SCIPhashtableInsert(sumcoefmap, (void *) ct) );
                           ct->color = (*nusedcolors)++;
                           color = ct->color;
                           nuniquecoefs++;
                        }
                        else
                           color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(sumcoefmap, (void *) ct))->color;

                        /* add the intermediate node with the corresponding color */
                        (void) G->add_vertex(color, degrees[n]);
                        internode = n++;

                        assert( internode < nnodes );

                        G->add_edge((unsigned) parentnode, (unsigned) internode);
#ifndef NDEBUG
                        ++m;
#endif
                        assert( m <= nedges );

                        /* connect the intermediate node to its corresponding variable node */
                        node = SCIPvarGetProbindex(vars[k]);
                        assert( node < nnodes );

                        G->add_edge((unsigned) node, (unsigned) internode);
#ifndef NDEBUG
                        ++m;
#endif
                        assert( m <= nedges );
                     }

                     /* add the node for the constant */
                     if ( constant != 0.0 )
                     {
                        SYM_CONSTTYPE* ct;

                        /* check whether we have to resize */
                        SCIP_CALL( SCIPensureBlockMemoryArray(scip, &uniqueconstarray, &constarraysize, nuniqueconsts+1) );
                        assert( nuniqueconsts < constarraysize );

                        ct = &uniqueconstarray[nuniqueconsts];
                        ct->value = constant;

                        if ( ! SCIPhashtableExists(consttypemap, (void *) ct) )
                        {
                           SCIP_CALL( SCIPhashtableInsert(consttypemap, (void *) ct) );
                           ct->color = (*nusedcolors)++;
                           color = ct->color;
                           nuniqueconsts++;
                        }
                        else
                           color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(consttypemap, (void *) ct))->color;

                        /* add the node with a new color */
                        (void) G->add_vertex(color, degrees[n]);
                        node = n++;

                        assert( node < nnodes );

                        G->add_edge((unsigned) parentnode, (unsigned) node);
#ifndef NDEBUG
                        ++m;
#endif
                        assert( m <= nedges );
                     }

                     /* add a filler node since it will be removed in the next iteration anyway */
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &visitednodes, &maxvisitednodes, numvisitednodes+1) );
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ischildofsum, &maxischildofsum, numischildofsum+1) );

                     visitednodes[numvisitednodes++] = n;
                     ischildofsum[numischildofsum++] = FALSE;
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

                  if ( ! SCIPhashtableExists(consttypemap, (void *) ct) )
                  {
                     SCIP_CALL( SCIPhashtableInsert(consttypemap, (void *) ct) );
                     ct->color = (*nusedcolors)++;
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

                  if ( ! SCIPhashtableExists(optypemap, (void *) ot) )
                  {
                     SCIP_CALL( SCIPhashtableInsert(optypemap, (void *) ot) );
                     ot->color = (*nusedcolors)++;
                     color = ot->color;
                     nuniqueops++;
                  }
                  else
                     color = ((SYM_OPTYPE*) SCIPhashtableRetrieve(optypemap, (void *) ot))->color;
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

                  if ( ! SCIPhashtableExists(rhstypemap, (void *) rt) )
                  {
                     SCIP_CALL( SCIPhashtableInsert(rhstypemap, (void *) rt) );
                     rt->color = (*nusedcolors)++;
                     parentcolor = rt->color;
                     nuniquerhs++;
                  }
                  else
                     parentcolor = ((SYM_RHSTYPE*) SCIPhashtableRetrieve(rhstypemap, (void *) rt))->color;

                  /* add the constraint side node with the corresponding color */
                  (void) G->add_vertex(parentcolor, degrees[n]);
                  parentnode = n++;
                  assert( parentnode < nnodes );
               }
               /* otherwise, get the parentnode stored in visitednodes */
               else
               {
                  parentnode = visitednodes[numvisitednodes - 1];
                  assert( parentnode < nnodes );
               }

               /* in all cases apart from variable expressions, the new node is added with the corresponding color */
               if ( color != -1 )
               {
                  (void) G->add_vertex(color, degrees[n]);
                  node = n++;
                  assert( node < nnodes );
                  assert( n <= nnodes );
               }

               /* store the new node so that it can be used as parentnode later */
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, &visitednodes, &maxvisitednodes, numvisitednodes+1) );
               SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ischildofsum, &maxischildofsum, numischildofsum+1) );

               assert( node != -1 );
               visitednodes[numvisitednodes++] = node;
               ischildofsum[numischildofsum++] = FALSE;

               /* connect the current node with its parent */
               assert( parentnode != -1 );

               if ( parentnode < node )
                  G->add_edge((unsigned) parentnode, (unsigned) node);
               else
                  G->add_edge((unsigned) node, (unsigned) parentnode);
#ifndef NDEBUG
               ++m;
#endif
               assert( m <= nedges );

               /* for sum expression, also add intermediate nodes for the coefficients */
               if ( SCIPisExprSum(scip, expr) )
               {
                  SCIP_Real* coefs;
                  SCIP_Real constval;
                  int internode;

                  coefs = SCIPgetCoefsExprSum(expr);

                  /* iterate over children from last to first, such that visitednodes array is in correct order */
                  for (j = SCIPexprGetNChildren(expr) - 1; j >= 0; --j)
                  {
                     SYM_CONSTTYPE* ct;

                     assert( nuniquecoefs < coefarraysize );

                     ct = &sumcoefarray[nuniquecoefs];
                     ct->value = coefs[j];

                     if ( ! SCIPhashtableExists(sumcoefmap, (void *) ct) )
                     {
                        SCIP_CALL( SCIPhashtableInsert(sumcoefmap, (void *) ct) );
                        ct->color = (*nusedcolors)++;
                        color = ct->color;
                        nuniquecoefs++;
                     }
                     else
                        color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(sumcoefmap, (void *) ct))->color;

                     /* add the intermediate node with the corresponding color */
                     (void) G->add_vertex(color, degrees[n]);
                     internode = n++;

                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &visitednodes, &maxvisitednodes, numvisitednodes+1) );
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &ischildofsum, &maxischildofsum, numischildofsum+1) );

                     visitednodes[numvisitednodes++] = internode;
                     ischildofsum[numischildofsum++] = TRUE;

                     assert( internode < nnodes );

                     G->add_edge((unsigned) node, (unsigned) internode);
#ifndef NDEBUG
                     ++m;
#endif
                     assert( m <= nedges );
                  }

                  /* add node for the constant term of the sum expression */
                  constval = SCIPgetConstantExprSum(expr);
                  if ( constval != 0.0 )
                  {
                     SYM_CONSTTYPE* ct;

                     /* check whether we have to resize */
                     SCIP_CALL( SCIPensureBlockMemoryArray(scip, &uniqueconstarray, &constarraysize, nuniqueconsts + 1) );
                     assert( nuniqueconsts < constarraysize );

                     ct = &uniqueconstarray[nuniqueconsts];
                     ct->value = constval;

                     if ( ! SCIPhashtableExists(consttypemap, (void *) ct) )
                     {
                        SCIP_CALL( SCIPhashtableInsert(consttypemap, (void *) ct) );
                        ct->color = (*nusedcolors)++;
                        color = ct->color;
                        nuniqueconsts++;
                     }
                     else
                        color = ((SYM_CONSTTYPE*) SCIPhashtableRetrieve(consttypemap, (void *) ct))->color;

                     /* add the node with a new color */
                     (void) G->add_vertex(color, degrees[n]);
                     internode = n++;

                     assert( node < nnodes );

                     G->add_edge((unsigned) node, (unsigned) internode);
#ifndef NDEBUG
                     ++m;
#endif
                     assert( m <= nedges );
                  }
               }

               ++currentlevel;
               break;
            }
            /* when leaving an expression, the nodes that are not needed anymore are erased from the respective arrays */
            case SCIP_EXPRITER_LEAVEEXPR:
            {
               --numvisitednodes;
               --numischildofsum;
               currentlevel--;

               /* When leaving the child of a sum expression, we have to pop again to get rid of the intermediate nodes
                * used for the coefficients of summands
                */
               if ( numischildofsum > 0 && ischildofsum[numischildofsum - 1] )
               {
                  --numvisitednodes;
                  --numischildofsum;
               }

               break;
            }

            default:
               SCIPABORT(); /* we should never be called in this stage */
               break;
         }
      }

      assert( currentlevel == 0 );
      assert( numvisitednodes == 0 );
      assert( numischildofsum == 0 );
   }
   assert( n == nnodes );
   assert( m == nedges );
   assert( n == matrixdata->npermvars + matrixdata->nrhscoef + nlinearnodes + nnonlinearnodes );
   assert( m == nlinearedges + nnonlinearedges );

   /* free everything */
   SCIPfreeBlockMemoryArrayNull(scip, &vals, varssize);
   SCIPfreeBlockMemoryArrayNull(scip, &vars, varssize);

   SCIPfreeBlockMemoryArray(scip, &visitednodes, maxvisitednodes);
   SCIPfreeBlockMemoryArray(scip, &ischildofsum, maxischildofsum);
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

/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return TRUE;
}

/** return name of external program used to compute generators */
char*
initStaticSymmetryName( )
{
   char* blissname = new char[100];
#ifdef BLISS_PATCH_PRESENT
   (void) SCIPsnprintf(blissname, 100, "bliss %sp", bliss::version);
#else
   (void) SCIPsnprintf(blissname, 100, "bliss %s", bliss::version);
#endif
   return blissname;
}

/** return name of external program used to compute generators */
char*
initStaticSymmetryAddName( )
{
   char* sassyname = new char[100];
   (void) SCIPsnprintf(sassyname, 100, "sassy %d.%d", SASSY_VERSION_MAJOR, SASSY_VERSION_MINOR);
   return sassyname;
}

static const char* symmetryname = initStaticSymmetryName();
static const char* symmetryaddname = initStaticSymmetryAddName();

/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
   return symmetryname;
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "Computing Graph Automorphisms by T. Junttila and P. Kaski (users.aalto.fi/~tjunttil/bliss/)";
}

/** return name of additional external program used for computing symmetries */
const char* SYMsymmetryGetAddName(void)
{
   return symmetryaddname;
}

/** return description of additional external program used to compute symmetries */
const char* SYMsymmetryGetAddDesc(void)
{
   return "Symmetry preprocessor by Markus Anders (github.com/markusa4/sassy)";
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
   SCIP_Real*            log10groupsize,     /**< pointer to store size of group */
   SCIP_Real*            symcodetime         /**< pointer to store the time for symmetry code */
   )
{
   SCIP_Real oldtime;
   SCIP_Bool success = FALSE;
   int* degrees;
   int maxdegrees;
   int nnodes;
   int nedges;
   int nlinearnodes;
   int nnonlinearnodes;
   int nlinearedges;
   int nnonlinearedges;
   int nusedcolors;

   assert( scip != NULL );
   assert( matrixdata != NULL );
   assert( exprdata != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( maxgenerators >= 0 );
   assert( symcodetime != NULL );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;
   *symcodetime = 0.0;

   /* determine number of nodes and edges */
   SCIP_CALL( determineGraphSize(scip, matrixdata, exprdata,
         &nnodes, &nedges, &nlinearnodes, &nnonlinearnodes, &nlinearedges, &nnonlinearedges,
         &degrees, &maxdegrees, &success) );

   if ( ! success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, 0, "Stopped symmetry computation: Symmetry graph would become too large.\n");
      return SCIP_OKAY;
   }

   /* create sassy graph */
   sassy::static_graph sassygraph;

   /* init graph */
   sassygraph.initialize_graph((unsigned) nnodes, (unsigned) nedges);

   /* add the nodes for linear and nonlinear constraints to the graph */
   SCIP_CALL( fillGraphByConss(scip, &sassygraph, matrixdata, exprdata,
         nnodes, nedges, nlinearnodes, nnonlinearnodes, nlinearedges, nnonlinearedges,
         degrees, &nusedcolors) );

   SCIPfreeBlockMemoryArray(scip, &degrees, maxdegrees);

   SCIPdebugMsg(scip, "Symmetry detection graph has %d nodes.\n", nnodes);

   /* init data */
   struct SYMMETRY_Data data;
   data.scip = scip;
   data.npermvars = matrixdata->npermvars;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.maxgenerators = maxgenerators;
   data.perms = NULL;

   oldtime = SCIPgetSolvingTime(scip);

   /* set up sassy preprocessor */
   sassy::preprocessor sassy;

   /* turn off some preprocessing that generates redudant permutations */
   sassy::configstruct sconfig;
   sconfig.CONFIG_PREP_DEACT_PROBE = true;
   sassy.configure(&sconfig);

   /* lambda function to have access to data and pass it to sassyhook above */
   sassy::sassy_hook sassyglue = [&](int n, const int* p, int nsupp, const int* suppa) {
      sassyhook((void*)&data, n, p, nsupp, suppa);
   };

   /* call sassy to reduce graph */
   sassy.reduce(&sassygraph, &sassyglue);

   /* create bliss graph */
   bliss::Graph blissgraph(0);

   /* convert sassy to bliss graph */
   convert_sassy_to_bliss(&sassygraph, &blissgraph);

#ifdef SCIP_OUTPUT
   blissgraph.write_dot("debug.dot");
#endif

#ifdef SCIP_DISABLED_CODE
   char filename[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s.dimacs", SCIPgetProbName(scip));
   FILE* fp = fopen(filename, "w");
   if ( fp )
   {
      blissgraph.write_dimacs(fp);
      fclose(fp);
   }
#endif


   /* compute automorphisms */
   bliss::Stats stats;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   blissgraph.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   blissgraph.set_component_recursion(false);

   /* do not use a node limit, but set generator limit */
#ifdef BLISS_PATCH_PRESENT
   blissgraph.set_search_limits(0, (unsigned) maxgenerators);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   /* lambda function to have access to stats and terminate the search if maxgenerators are reached */
   auto term = [&]() {
      return (stats.get_nof_generators() >= (unsigned) maxgenerators);
   };

   auto hook = [&](unsigned int n, const unsigned int* aut) {
      sassy.bliss_hook(n, aut);
   };

   /* start search */
   blissgraph.find_automorphisms(stats, hook, term);
#else
   /* start search */
   blissgraph.find_automorphisms(stats, sassy::preprocessor::bliss_hook, (void*) &sassy);
#endif
   *symcodetime = SCIPgetSolvingTime(scip) - oldtime;

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
