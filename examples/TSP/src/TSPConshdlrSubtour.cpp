/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: TSPConshdlrSubtour.cpp,v 1.1 2005/03/03 16:43:34 bzfberth Exp $"

/**@file   TSPReader.cpp
 * @brief  C++ file reader for TSP data files
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <string>
#include <sstream>

#include "TSPConshdlrSubtour.h"
#include "gminucut.h"

extern "C" {
#include "scip/cons_linear.h"
}
using namespace tsp;
using namespace scip;
using namespace std;


struct ConsData
{
   GRAPH* graph;
};
 


bool findSubtour( 
   SCIP*         scip,               /**< SCIP data structure */
   GRAPH*        graph,              /**< underlying graph */
   SOL*          sol,                /**< proposed solution */
   bool*         subtour             /**< if a subtour elimination cut is to be created, the nodes in a subtour are 
                                      * labeled true */
   )
{  
   if(graph->nnodes <= 1)
      return false;
   int tourlength;
   GRAPHNODE* node;
   GRAPHNODE* startnode;
   GRAPHEDGE* startedge;
   GRAPHEDGE* lastedge;
   GRAPHEDGE* edge;

   startnode = &graph->nodes[0];

   //has to be restarted at maximum once
 SUBTOURLOOP:

   // if a subtour exits, every node is part of some subtour 
   if(subtour != NULL)
   {
      assert(0 <= startnode->id && startnode->id < graph->nnodes);
      clearMemoryArray(subtour, graph->nnodes);
      subtour[startnode->id] = true;
   }
   tourlength = 0;
   lastedge = NULL;
   node = startnode;

   // walk through the tour until you come back to the startnode
   do
   {
      startedge = node->first_edge;
      assert(startedge != NULL);
      edge = startedge;

      // look for an outgoing edge to proceed
      do
      {
         // if a new edge with value numerical equal to one is found, we proceed
         if( edge->back != lastedge && SCIPgetSolVal(scip, sol, edge->var) > 0.5 )
         {
            tourlength++;
            node = edge->adjac;
            lastedge = edge;
            
            // the new node is stored to be part of a possible subtour
            if( subtour != NULL )
            {
               assert(0 <= node->id && node->id < graph->nnodes);
               assert(node == startnode || !subtour[node->id]);
               subtour[node->id] = true;
            }
            edge = NULL;
            if( tourlength > graph->nnodes )
            {
               /* we found a subtour without the starting node, e.g. 0 - 1 - 2 - 3 - 1 - 2 - ...;
                * this can only happen, if the degree constraints are violated;
                * start again with the last visited node as starting node, because this must be member of the subtour;
                * thus, in the second run we will find the subtour!
                */
               startnode = node;

               goto SUBTOURLOOP;
            }
            break;
         }    
         edge = edge->next;        
      }
      while( edge != startedge );

      /* we didn't find an outgoing edge in the solution: the degree constraints must be violated; abort! */
      if( edge == startedge )
         return false;
      
   }
   while( node != startnode );
  
   assert(tourlength <= graph->nnodes);

   if( graph->nnodes != tourlength )
      return true;
   else 
      return false; 
}

/** frees specific constraint data
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 */
RETCODE TSPConshdlrSubtour::scip_delete(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONSDATA**    consdata            /**< pointer to the constraint data to free */
   )
{
   assert(consdata != NULL);

   release_graph(&(*consdata)->graph);
   SCIPfreeMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
RETCODE TSPConshdlrSubtour::scip_trans(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS*         sourcecons,         /**< source constraint to transform */
   CONS**        targetcons          /**< pointer to store created target constraint */
   )
{
   CONSDATA* sourcedata;
   CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);

   CHECK_OKAY( SCIPallocMemory(scip,&targetdata) );
   targetdata->graph = sourcedata->graph;
   capture_graph(targetdata->graph);

   /* create target constraint */
   CHECK_OKAY( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), SCIPconsIsRemoveable(sourcecons)) );

   return SCIP_OKAY;
}

/** separation method of constraint handler
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
RETCODE TSPConshdlrSubtour::scip_sepa(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS**        conss,              /**< array of constraints to process */
   int           nconss,             /**< number of constraints to process */
   int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   RESULT*       result              /**< pointer to store the result of the separation call */
   )
{
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   for( int c = 0; c < nusefulconss; ++c )
   {
      // get all required structures
      CONSDATA* consdata;
      GRAPH* graph;
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);
  
      double cutvalue;
      long cutsize;

      // store the suggested, but infeasible solution into the capacity of the edges
      for( int i = 0; i < graph->nedges; i++)
      {
         graph->edges[i].rcap = SCIPgetVarSol(scip, graph->edges[i].var);
         graph->edges[i].cap = graph->edges[i].rcap;
      }
     
      // try to cut of this solution
      if( gmincut(graph, &cutvalue, &cutsize) && cutvalue < 2.0 )
      {
         
         ROW* row;
         stringstream name;

         // a new seperation constraint is created out of a partition of the graph with a cut value less than 2  
         name << "sepa_con";
         CHECK_OKAY( SCIPcreateEmptyRow(scip, &row, name.str().c_str(), 2.0, SCIPinfinity(scip), 
               FALSE, FALSE, TRUE) ); 
         for( int i = 0; i < graph->nnodes;i++)
         { 
            // in gmincut the graph has been partitioned into two parts
            if( graph->nodes[i].partition )
            {
               GRAPHEDGE* edge;
               GRAPHEDGE* startedge;

               edge = graph->nodes[i].first_edge;
               startedge = edge;

               // take every edge with nodes in different parts into account
               do
               {
                  if( !edge->adjac->partition )
                  {
                     CHECK_OKAY( SCIPaddVarToRow(scip, row, edge->var, 1.0) );
                  }
                  edge = edge->next;
               }
               while( edge != startedge );
            }
         }
        
         // add cut
         if( SCIPisCutEfficacious(scip, row) )
         {
            CHECK_OKAY( SCIPaddCut(scip, row, FALSE) );
            *result = SCIP_SEPARATED;    
         }
         CHECK_OKAY( SCIPreleaseRow(scip, &row) );
      }
          
   }

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
RETCODE TSPConshdlrSubtour::scip_enfolp(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS**        conss,              /**< array of constraints to process */
   int           nconss,             /**< number of constraints to process */
   int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   RESULT*       result              /**< pointer to store the result of the enforcing call */
   )
{
   bool* subtour;

   *result = SCIP_FEASIBLE;
   subtour = NULL;
   for( int i = 0; i < nconss; ++i )
   {
      CONSDATA* consdata;
      GRAPH* graph;
      Bool found;
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);

      CHECK_OKAY( SCIPreallocBufferArray(scip, &subtour, graph->nnodes) );
      found = findSubtour(scip, graph, NULL, subtour);
      
      // if a subtour was found, we generate a cut constraint saying that there must be at least two outgoing edges
      if( found )
      {
         ROW* row;
         stringstream name;

         // a new cut constraint is created 
         name << "loop_con";
         CHECK_OKAY( SCIPcreateEmptyRow(scip, &row, name.str().c_str(), 2.0, SCIPinfinity(scip), 
               FALSE, FALSE, TRUE) ); 
         for( int j = i+1; j < graph->nnodes; j++ )
         {
            if( subtour[j] )
            {
               GRAPHEDGE* edge;
               GRAPHEDGE* startedge;
               
               edge = graph->nodes[j].first_edge;
               startedge = edge;

               // find edges going out of the subtour
               do
               {
                  if( !subtour[edge->adjac->id] )
                  {
                     CHECK_OKAY( SCIPaddVarToRow(scip, row, edge->var, 1.0) );
                  }
                  edge = edge->next;
               }
               while( edge != startedge );
            }
         }
         
         // add the constraint to SCIP
         CHECK_OKAY( SCIPaddCut(scip, row, FALSE) );
         CHECK_OKAY( SCIPreleaseRow(scip, &row) );
         
         *result = SCIP_SEPARATED;
      }
   }
   SCIPfreeBufferArrayNull(scip, &subtour);


   /* weil die enfoprio des conshdlrs < 0 ist, sind hier die Variablen in jedem Falle ganzzahlig;
    * von knoten 0 aus kreis laufen und gucken wie lang der ist;
    * falls kleiner n: schnitt gefunden, einfuegen, *result = SCIP_SEPARATED
    */

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
 *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the LP
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is true)
 */
RETCODE TSPConshdlrSubtour::scip_enfops(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS**        conss,              /**< array of constraints to process */
   int           nconss,             /**< number of constraints to process */
   int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   RESULT*       result              /**< pointer to store the result of the enforcing call */
   )
{
   *result = SCIP_FEASIBLE;
 
   for( int i = 0; i < nconss; ++i )
   {
      CONSDATA* consdata;
      GRAPH* graph;
      Bool found;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);
  
      //if a subtour is found, the solution must be infeasible
      found = findSubtour(scip, graph, NULL, NULL);      
      if( found )
      {
         *result = SCIP_INFEASIBLE;
      }
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
RETCODE TSPConshdlrSubtour::scip_check(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS**        conss,              /**< array of constraints to process */
   int           nconss,             /**< number of constraints to process */
   SOL*          sol,                /**< the solution to check feasibility for */
   Bool          checkintegrality,   /**< has integrality to be checked? */
   Bool          checklprows,        /**< have current LP rows to be checked? */
   RESULT*       result              /**< pointer to store the result of the feasibility checking call */
   )
{
   *result = SCIP_FEASIBLE;

   for( int i = 0; i < nconss; ++i )
   {
      CONSDATA* consdata;
      GRAPH* graph;
      Bool found;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);
      graph = consdata->graph;
      assert(graph != NULL);
     
      // if a subtour is found, the solution must be infeasible
      found = findSubtour(scip, graph, sol, NULL);      
      if( found )
      {
         *result = SCIP_INFEASIBLE;
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
RETCODE TSPConshdlrSubtour::scip_prop(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS**        conss,              /**< array of constraints to process */
   int           nconss,             /**< number of constraints to process */
   int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   RESULT*       result              /**< pointer to store the result of the propagation call */
   )
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
RETCODE TSPConshdlrSubtour::scip_lock(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
                                      *   constraint handler does not need constraints */
   int           nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
   int           nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
   )
{
   CONSDATA* consdata;
   GRAPH* g;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   
   g = consdata->graph;
   assert(g != NULL);

   for( int i = 0; i < g->nedges; ++i )
   {
      CHECK_OKAY( SCIPaddVarLocks(scip, g->edges[i].var, nlocksneg, nlockspos) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler
 *
 *  The constraint handler should store a representation of the constraint into the given text file.
 */
RETCODE TSPConshdlrSubtour::scip_print(
   SCIP*         scip,               /**< SCIP data structure */
   CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   CONS*         cons,               /**< the constraint that should be displayed */
   FILE*         file                /**< the text file to store the information into */
   )
{
   CONSDATA* consdata;
   GRAPH* g;
      
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
      
   g = consdata->graph;
   assert(g != NULL);

   fprintf(file, "subtour of Graph G with %d nodes and %d edges", g->nnodes, g->nedges);

   return SCIP_OKAY;
}



/** creates and captures a TSP subtour constraint */
RETCODE tsp::SCIPcreateConsSubtour(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   GRAPH*           graph,              /**< the underlying graph */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   /* find the knapsack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "subtour");
   if( conshdlr == NULL )
   {
      errorMessage("subtour constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   CHECK_OKAY( SCIPallocMemory( scip, &consdata) );
   consdata->graph = graph;
   capture_graph( consdata->graph );

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, removeable) );

   return SCIP_OKAY;
}
