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
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   ConshdlrSubtour.cpp
 * @brief  Subtour elimination constraint handler for TSP problems, written in C++
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <string>
#include <iostream>
#include "ConshdlrSubtour.h"
#include "GomoryHuTree.h"

#include "objscip/objscip.h"

#include "scip/cons_linear.h"

using namespace tsp;
using namespace scip;
using namespace std;

struct SCIP_ConsData
{
   GRAPH* graph;
};

/* checks whether proposed solution contains a subtour */
static
SCIP_Bool findSubtour( 
   SCIP*              scip,               /**< SCIP data structure */
   GRAPH*             graph,              /**< underlying graph */
   SCIP_SOL*          sol                 /**< proposed solution */
   )
{  
   GRAPHNODE* node;
   GRAPHNODE* startnode;
   GRAPHEDGE* lastedge;
   GRAPHEDGE* edge;
   GRAPHEDGE* nextedge;
   int tourlength;
   SCIP_Bool foundnextedge;

   if(graph->nnodes <= 1)
      return FALSE;

   startnode = &graph->nodes[0];

   tourlength = 0;
   lastedge = NULL;
   node = startnode;

   // follow the (sub?)tour until you come back to the startnode
   do
   {
      edge = node->first_edge;
      foundnextedge = FALSE;
      nextedge = NULL;

      // look for an outgoing edge to proceed
      while( edge != NULL )
      {
         // if a new edge with value numerical equal to one is found, we proceed
         if( edge->back != lastedge && SCIPgetSolVal(scip, sol, edge->var) > 0.5 )
         {          
            tourlength++;
            
            if( foundnextedge || tourlength > graph->nnodes )
            {
               /* we found a subtour without the starting node, e.g. 0 - 1 - 2 - 3 - 1 - 2 - ...;
                * this can only happen, if the degree constraints are violated;
                * start again with the last visited node as starting node, because this must be member of the subtour;
                * thus, in the second run we will find the subtour!
                */
               return TRUE;
            }

            foundnextedge= TRUE;
            nextedge = edge;            
            
            if( node == startnode )
               break;
         }    
      
         edge = edge->next;        
      }
   
      /* we didn't find an outgoing edge in the solution: the degree constraints must be violated; abort! */
      if( nextedge == NULL )
         return TRUE;

      node = nextedge->adjac;
      lastedge = nextedge;
   }
   while( node != startnode );

   assert(tourlength <= graph->nnodes);

   return ( graph->nnodes != tourlength );
}

/* separates subtour elemination cuts */
static
SCIP_RETCODE sepaSubtour(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints to process */
   int                nconss,             /**< number of constraints to process */
   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*          sol,                /**< primal solution that should be separated */
   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
   )
{
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   for( int c = 0; c < nusefulconss; ++c )
   {
      // get all required structures
      SCIP_CONSDATA* consdata;
      GRAPH* graph;
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);
   
      double cap;
       
      // store the suggested, but infeasible solution into the capacity of the edges
      for( int i = 0; i < graph->nedges; i++)
      {
         cap = SCIPgetSolVal(scip, sol, graph->edges[i].var);
         graph->edges[i].rcap = cap;
         graph->edges[i].cap = cap;
         graph->edges[i].back->rcap = cap;
         graph->edges[i].back->cap = cap;   
      }
           
      SCIP_Bool** cuts;
      int ncuts;

      SCIP_CALL( SCIPallocBufferArray(scip, &cuts, graph->nnodes) );
      for(int i = 0; i < graph->nnodes; i++)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &cuts[i], graph->nnodes) );
      }

      // try to find cuts
      if( ghc_tree( graph, cuts, &ncuts, SCIPfeastol(scip) ) )
      { 
         int i = 0;

         // create a new cutting plane for every suitable arc (representing a cut with value < 2) of the Gomory Hu Tree
         while( i < ncuts )
         {
            SCIP_ROW* row; 
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "sepa_con", 2.0, SCIPinfinity(scip), FALSE, FALSE, TRUE) ); 

            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            for( int j = 0; j < graph->nnodes; j++)
            { 
               // in gmincut the graph has been partitioned into two parts, represented by bools
               if( cuts[i][j] )
               {
                  GRAPHEDGE* edge = graph->nodes[j].first_edge;
                                        
                  // take every edge with nodes in different parts into account
                  while( edge != NULL )
                  {
                     if( !cuts[i][edge->adjac->id] )
                     {
                        SCIP_CALL( SCIPaddVarToRow(scip, row, edge->var, 1.0) ); 
                     }
                     edge = edge->next;
                  }
               }
            } 

            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            // add cut
            if( SCIPisCutEfficacious(scip, sol, row) )
            {
               SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
               *result = SCIP_SEPARATED;    
            }
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
               
            i++;
         }            
      }
      for( int i = graph->nnodes - 1; i >= 0; i-- )
         SCIPfreeBufferArray( scip, &cuts[i] );
      SCIPfreeBufferArray( scip, &cuts );
                
   }

   return SCIP_OKAY;
}


/** frees specific constraint data
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 */
SCIP_DECL_CONSDELETE(ConshdlrSubtour::scip_delete)
{
   assert(consdata != NULL);

   release_graph(&(*consdata)->graph);
   SCIPfreeMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrSubtour::scip_trans)
{
   SCIP_CONSDATA* sourcedata = NULL;
   SCIP_CONSDATA* targetdata = NULL;

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocMemory(scip, &targetdata) );
   targetdata->graph = sourcedata->graph;
   capture_graph(targetdata->graph);

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solution
 *
 *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
 *  which means that a valid LP solution exists.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
SCIP_DECL_CONSSEPALP(ConshdlrSubtour::scip_sepalp)
{
   SCIP_CALL( sepaSubtour(scip, conshdlr, conss, nconss, nusefulconss, NULL, result) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solution
 *
 *  Separates all constraints of the constraint handler. The method is called outside the LP solution loop (e.g., by
 *  a relaxator or a primal heuristic), which means that there is no valid LP solution.
 *  Instead, the method should produce cuts that separate the given solution.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 *  - SCIP_DELAYED    : the separator was skipped, but should be called again
 */
SCIP_DECL_CONSSEPASOL(ConshdlrSubtour::scip_sepasol)
{
   SCIP_CALL( sepaSubtour(scip, conshdlr, conss, nconss, nusefulconss, sol, result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was solved.
 *  The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or separating the solution with a valid
 *  cutting plane.
 *
 *  The enforcing methods of the active constraint handlers are called in decreasing order of their enforcing
 *  priorities until the first constraint handler returned with the value SCIP_CUTOFF, SCIP_SEPARATED,
 *  SCIP_REDUCEDDOM, SCIP_CONSADDED, or SCIP_BRANCHED.
 *  The integrality constraint handler has an enforcing priority of zero. A constraint handler which can
 *  (or wants) to enforce its constraints only for integral solutions should have a negative enforcing priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to incorporate its own branching strategy even on non-integral
 *  solutions must have an enforcing priority greater than zero (e.g. the SOS-constraint incorporates
 *  SOS-branching on non-integral solutions).
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
 *  be enforced, if no violation was found in the useful constraints.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
SCIP_DECL_CONSENFOLP(ConshdlrSubtour::scip_enfolp)
{
   *result = SCIP_FEASIBLE;

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_CONSDATA* consdata;
      GRAPH* graph;
      SCIP_Bool found;
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);

      found = findSubtour(scip, graph, NULL);
      
      // if a subtour was found, we generate a cut constraint saying that there must be at least two outgoing edges
      if( found )
         *result = SCIP_INFEASIBLE;
   }
   
   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was not solved.
 *  The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or adding an additional constraint.
 *  Separation is not possible, since the LP is not processed at the current node. All LP informations like
 *  LP solution, slack values, or reduced costs are invalid and must not be accessed.
 *
 *  Like in the enforcing method for LP solutions, the enforcing methods of the active constraint handlers are
 *  called in decreasing order of their enforcing priorities until the first constraint handler returned with
 *  the value SCIP_CUTOFF, SCIP_REDUCEDDOM, SCIP_CONSADDED, SCIP_BRANCHED, or SCIP_SOLVELP.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
 *  be enforced, if no violation was found in the useful constraints.
 *
 *  If the pseudo solution's objective value is lower than the lower bound of the node, it cannot be feasible
 *  and the enforcing method may skip it's check and set *result to SCIP_DIDNOTRUN. However, it can also process
 *  its constraints and return any other possible result code.
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the SCIP_LP
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is true)
 */
SCIP_DECL_CONSENFOPS(ConshdlrSubtour::scip_enfops)
{
   *result = SCIP_FEASIBLE;
 
   for( int i = 0; i < nconss; ++i )
   {
      SCIP_CONSDATA* consdata;
      GRAPH* graph;
      SCIP_Bool found;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);
  
      // if a subtour is found, the solution must be infeasible
      found = findSubtour(scip, graph, NULL);      
      if( found )
         *result = SCIP_INFEASIBLE;
   }

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions
 *
 *  The given solution has to be checked for feasibility.
 *  
 *  The check methods of the active constraint handlers are called in decreasing order of their check
 *  priorities until the first constraint handler returned with the result SCIP_INFEASIBLE.
 *  The integrality constraint handler has a check priority of zero. A constraint handler which can
 *  (or wants) to check its constraints only for integral solutions should have a negative check priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to check feasibility even on non-integral solutions must have a
 *  check priority greater than zero (e.g. if the check is much faster than testing all variables for
 *  integrality).
 *
 *  In some cases, integrality conditions or rows of the current LP don't have to be checked, because their
 *  feasibility is already checked or implicitly given. In these cases, 'checkintegrality' or
 *  'checklprows' is FALSE.
 *
 *  possible return values for *result:
 *  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
SCIP_DECL_CONSCHECK(ConshdlrSubtour::scip_check)
{
   *result = SCIP_FEASIBLE;

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_CONSDATA* consdata;
      GRAPH* graph;
      SCIP_Bool found;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);
     
      // if a subtour is found, the solution must be infeasible
      found = findSubtour(scip, graph, sol);      
      if( found )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
            SCIPinfoMessage(scip, NULL, "violation: graph has a subtour\n");
         }
      }
   }   


   return SCIP_OKAY;
}

/** domain propagation method of constraint handler
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched, but did not find any domain reductions
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
 */
SCIP_DECL_CONSPROP(ConshdlrSubtour::scip_prop)
{
   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler
 *
 *  This method is called, after a constraint is added or removed from the transformed problem.
 *  It should update the rounding locks of all associated variables with calls to SCIPaddVarLocks(),
 *  depending on the way, the variable is involved in the constraint:
 *  - If the constraint may get violated by decreasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the
 *    (positive) constraint infeasible and rounding up is potentially rendering the negation of the constraint
 *    infeasible.
 *  - If the constraint may get violated by increasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlocksneg, nlockspos), saying that rounding up is potentially rendering the
 *    constraint's negation infeasible and rounding up is potentially rendering the constraint itself
 *    infeasible.
 *  - If the constraint may get violated by changing the variable in any direction, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg).
 *
 *  Consider the linear constraint "3x -5y +2z <= 7" as an example. The variable rounding lock method of the
 *  linear constraint handler should call SCIPaddVarLocks(scip, x, nlocksneg, nlockspos), 
 *  SCIPaddVarLocks(scip, y, nlockspos, nlocksneg) and SCIPaddVarLocks(scip, z, nlocksneg, nlockspos) to tell SCIP,
 *  that rounding up of x and z and rounding down of y can destroy the feasibility of the constraint, while rounding
 *  down of x and z and rounding up of y can destroy the feasibility of the constraint's negation "3x -5y +2z > 7".
 *  A linear constraint "2 <= 3x -5y +2z <= 7" should call
 *  SCIPaddVarLocks(scip, ..., nlockspos + nlocksneg, nlockspos + nlocksneg) on all variables, since rounding in both
 *  directions of each variable can destroy both the feasibility of the constraint and it's negation
 *  "3x -5y +2z < 2  or  3x -5y +2z > 7".
 *
 *  If the constraint itself contains other constraints as sub constraints (e.g. the "or" constraint concatenation
 *  "c(x) or d(x)"), the rounding lock methods of these constraints should be called in a proper way.
 *  - If the constraint may get violated by the violation of the sub constraint c, it should call
 *    SCIPaddConsLocks(scip, c, nlockspos, nlocksneg), saying that infeasibility of c may lead to infeasibility of
 *    the (positive) constraint, and infeasibility of c's negation (i.e. feasibility of c) may lead to infeasibility
 *    of the constraint's negation (i.e. feasibility of the constraint).
 *  - If the constraint may get violated by the feasibility of the sub constraint c, it should call
 *    SCIPaddConsLocks(scip, c, nlocksneg, nlockspos), saying that infeasibility of c may lead to infeasibility of
 *    the constraint's negation (i.e. feasibility of the constraint), and infeasibility of c's negation (i.e. feasibility
 *    of c) may lead to infeasibility of the (positive) constraint.
 *  - If the constraint may get violated by any change in the feasibility of the sub constraint c, it should call
 *    SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg).
 *
 *  Consider the or concatenation "c(x) or d(x)". The variable rounding lock method of the or constraint handler
 *  should call SCIPaddConsLocks(scip, c, nlockspos, nlocksneg) and SCIPaddConsLocks(scip, d, nlockspos, nlocksneg)
 *  to tell SCIP, that infeasibility of c and d can lead to infeasibility of "c(x) or d(x)".
 *
 *  As a second example, consider the equivalence constraint "y <-> c(x)" with variable y and constraint c. The
 *  constraint demands, that y == 1 if and only if c(x) is satisfied. The variable lock method of the corresponding
 *  constraint handler should call SCIPaddVarLocks(scip, y, nlockspos + nlocksneg, nlockspos + nlocksneg) and
 *  SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg), because any modification to the
 *  value of y or to the feasibility of c can alter the feasibility of the equivalence constraint.
 */
SCIP_DECL_CONSLOCK(ConshdlrSubtour::scip_lock)
{
   SCIP_CONSDATA* consdata;
   GRAPH* g;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   g = consdata->graph;
   assert(g != NULL);

   for( int i = 0; i < g->nedges; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, g->edges[i].var, nlocksneg, nlockspos) );
   }

   return SCIP_OKAY;
}

/** variable deletion method of constraint handler
 *
 *  This method should iterate over all constraints of the constraint handler and delete all variables
 *  that were marked for deletion by SCIPdelVar().
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints in transformed problem
 *  - nconss          : number of constraints in transformed problem
 */
SCIP_DECL_CONSDELVARS(ConshdlrSubtour::scip_delvars)
{
   return SCIP_OKAY;
}


/** constraint display method of constraint handler
 *
 *  The constraint handler should store a representation of the constraint into the given text file.
 */
SCIP_DECL_CONSPRINT(ConshdlrSubtour::scip_print)
{
   SCIP_CONSDATA* consdata;
   GRAPH* g;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
      
   g = consdata->graph;
   assert(g != NULL);

   SCIPinfoMessage(scip, file, "subtour of Graph G with %d nodes and %d edges\n", g->nnodes, g->nedges);
   
   return SCIP_OKAY;
}

/** clone method which will be used to copy a objective plugin */
SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable* ConshdlrSubtour::clone)
{
   *valid = true;
   return new ConshdlrSubtour(scip);
}

/** constraint copying method of constraint handler
 *
 *  The constraint handler can provide a copy method, which copies a constraint from one SCIP data structure into a other
 *  SCIP data structure.
 */
SCIP_DECL_CONSCOPY(ConshdlrSubtour::scip_copy)
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSDATA* consdata = NULL;

   /* find the subtour constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "subtour");
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("subtour constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   ProbDataTSP * probdatatsp = NULL;
   probdatatsp = dynamic_cast<ProbDataTSP *>(SCIPgetObjProbData(scip));
   assert( probdatatsp != NULL );
   GRAPH * graph = probdatatsp->getGraph();
   consdata->graph = graph;
   capture_graph( consdata->graph );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, (name == NULL) ? SCIPconsGetName(sourcecons) : name, 
         conshdlr, consdata, initial, separate, enforce, check, 
         propagate, local, modifiable, dynamic, removable, FALSE) );

   *valid = true;
   return SCIP_OKAY;
}

/** creates and captures a TSP subtour constraint */
SCIP_RETCODE tsp::SCIPcreateConsSubtour(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   GRAPH*                graph,              /**< the underlying graph */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSDATA* consdata = NULL;

   /* find the subtour constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "subtour");
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("subtour constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   consdata->graph = graph;
   capture_graph( consdata->graph );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   return SCIP_OKAY;
}
