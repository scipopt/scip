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

/**@file   cons_storeGraph.c
 * @brief  constraint handler for storing the graph at each node of the tree
 * @author Gerald Gamrath
 *
 * This file implements the constraints that are used for the branching in the coloring algorithm.
 *
 * For each node in the branch-and-bound tree, a constraint of this type is created, which stores
 * all restrictions related to that branch-and-bound node.
 *
 * First of all, it stores the type of the constraint ("same" or "differ", the root has type root)
 * and the two nodes in the graph on which this restriction is applied.  When the branch-and-bound
 * node corresponding to the constraint is examined for the first time, the constraint creates a
 * graph that takes into account all the restrictions, which are active at this node.
 * At the root, this is the original (preprocessed) graph.  At any other branch-and-bound node, it
 * takes the graph of the constraint related to the branch-and-bound parent node of the current node and
 * modifies it so that all restrictions up to this node are respected.  Since the graph in the
 * branch-and-bound parent respects all restrictions on the path to that node, only the last
 * requirement, the one saved at the current branch-and-bound node, must be added.
 * This is done as follows: Adding a DIFFER(v,w) constraint is easy, since it suffices to add
 * an edge between v and w. For a SAME(v,w) constraint, the original idea is to collapse the nodes v
 * and w into one single vertex. Since this is not possible in the tclique-graph data structure, we
 * introduce new edges in the graph, so that v and w have the same neighborhood.  Hence, in the
 * pricing routine, each new stable set will either contain both nodes or none of them, since we
 * create (inclusion-) maximal sets.
 *
 * This does of course not hold for sets created in a higher level of the branch-and-bound tree or
 * in another subtree. In order to forbid all of these sets, which do not fulfill the current
 * restrictions, a propagation is started when the node is entered the first time and repeated
 * later, if the node is reentered after the creation of new variables in another subtree. The
 * propagation simply fixes  all variables to 0 which represent a stable set that does not
 * fulfill the restriction at the current node.
 *
 * The information about all fusions of nodes (caused by the SAME() operation) is stored, so that the nodes
 * constituting a union can be accessed easily. Each union has a representative and a set of nodes, whereas
 * each node knows the representative of the union it belongs to. At the beginning, each node forms its own
 * union and therefore each node also represents this union, consisting of only this node.  Later on, some
 * nodes represent unions of several nodes, while other nodes are part of a union which they do not represent,
 * so they have another node as representative. The representatives of the nodes are returned by the methods
 * COLORconsGetRepresentative() / COLORconsGetRepresentatives(), the union represented by a node is returned
 * by COLORconsGetUnion(), the array of unions, indexed by the representing node, is returned by
 * COLORconsGetUnions().
 */

#include <assert.h>
#include <string.h>

#include "scip/type_cons.h"
#include "cons_storeGraph.h"
#include "probdata_coloring.h"
#include "tclique/tclique.h"
#include "reader_col.h"
#include "scip/cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "storeGraph"
#define CONSHDLR_DESC          "storing graph at nodes of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              * propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP


/** constraint data for storing graph constraints */
struct SCIP_ConsData
{
   TCLIQUE_GRAPH*     graph;                 /* the current graph in the B&B-node belonging to this constraint */
   TCLIQUE_GRAPH*     cgraph;                /* the complementary graph of the current graph */
   SCIP_CONS*         fathercons;            /* the constraint sticking at the B&B-node's father */
   int*               representativeofnode;  /* r...[i] = j if node j is representative of the union containing node i */
   int**              unionofnode;           /* for all represantatives of a union an array with all the union's members */
   int*               nnodesinunion;         /* value at position i = #elements in unionofnode[i] */
   int                node1;                 /* first node for DIFFER / SAME */
   int                node2;                 /* second node for DIFFER / SAME */
   int                type;                  /* type of the branching operation: COLOR_CONSTYPE_DIFFER oder COLOR_CONSTYPE_SAME */
   int                propagatedvars;        /* number of Vars that existed, the last time, the related node was propagated,
                                                used to determine whether the constraint should be repropagated*/
   SCIP_Bool          created;               /* flag for saving the creation status of the graph saved in the cons,
                                                at the beginning false, after the first activation set to true */
   SCIP_NODE*         stickingatnode;        /* the node in the B&B-tree at which the cons is sticking */
};



/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONS**        stack;                 /**< stack for storing active constraints */
   int                nstack;                /**< number of elements on the stack */
   int                maxstacksize;          /**< maximum size of the stack */
};


/*
 * Local methods
 */

/** creates and captures the storeGraph constraint for the root node*/
static
SCIP_RETCODE createConsStoreGraphAtRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   TCLIQUE_GRAPH*        graph               /**< the original graph */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;
   int nnodes;

   assert(scip != NULL);
   assert(graph != NULL);
   nnodes = tcliqueGetNNodes(graph);
   /* find the storeGraph constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   SCIPdebugMessage("Creating graph storage constraint at root node.\n");

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->graph = graph;
   consdata->node1 = -1;
   consdata->node2 = -1;
   consdata->type = COLOR_CONSTYPE_ROOT;
   consdata->fathercons = NULL;
   consdata->propagatedvars = 0;
   consdata->stickingatnode = NULL;
   consdata->created = TRUE;

   /* allocate memory for the arrays and fill them */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->representativeofnode), nnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->nnodesinunion), nnodes) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->unionofnode), nnodes) );
   for ( i = 0; i < nnodes; i++ )
   {
      consdata->representativeofnode[i] = i;
      consdata->nnodesinunion[i] = 1;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->unionofnode[i]), 1) );
      consdata->unionofnode[i][0] = i;
   }

   /* create the complementary graph */
   tcliqueCreate(&(consdata->cgraph));
   assert(consdata->cgraph != NULL);

   COLORprobGetComplementaryGraph(scip, graph, consdata->cgraph);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
			     TRUE, FALSE, TRUE, FALSE, FALSE));

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
/** We do not want to copy store graph constraints into subSCIPs since they just store information about
 *  branching decisions and are used to enforce those.
 *  However, in subSCIPs, we only want to solve the current MIP with a branch-and-cut approach.
 */
#define conshdlrCopyStoreGraph NULL

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeStoreGraph)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   SCIPdebugMessage("freeing store graph constraint handler\n");

   /* free constraint handler storage */
   assert(conshdlrData->stack == NULL);
   SCIPfreeMemory(scip, &conshdlrData);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolStoreGraph)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);

   /* prepare stack */
   SCIP_CALL( SCIPallocMemoryArray(scip, &conshdlrData->stack, conshdlrData->maxstacksize) );
   SCIP_CALL( createConsStoreGraphAtRoot(scip, &cons, "root", COLORprobGetGraph(scip)) );

   /* release constraints */
   conshdlrData->stack[0] = cons;
   conshdlrData->nstack = 1;

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolStoreGraph)
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->nstack == 1); /* at this point the stack should only have the root-constraint on it */
   SCIP_CALL( SCIPreleaseCons(scip, &(conshdlrData->stack[0])) );
   conshdlrData->stack[0] = NULL;
   SCIPdebugMessage("exiting store graph constraint handler\n");

   /* free stack */
   SCIPfreeMemoryArray(scip, &conshdlrData->stack);

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteStoreGraph)
{
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(*consdata != NULL);

   SCIPdebugMessage("Deleting store graph constraint: <%s(%d,%d)>.\n", SCIPconsGetName(cons), (*consdata)->node1+1, (*consdata)->node2+1);

   /* free constraint data */
   if ( (*consdata)->type == COLOR_CONSTYPE_ROOT )
   {
      for ( i = tcliqueGetNNodes((*consdata)->graph)-1; i >= 0; i-- )
      {
         SCIPfreeBlockMemoryArray(scip, &((*consdata)->unionofnode[i]), (*consdata)->nnodesinunion[i]);
         assert((*consdata)->nnodesinunion[i] == 1);
      }
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->unionofnode), tcliqueGetNNodes((*consdata)->graph));
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->nnodesinunion), tcliqueGetNNodes((*consdata)->graph));
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->representativeofnode), tcliqueGetNNodes((*consdata)->graph));
      tcliqueFree(&((*consdata)->cgraph));
   }
   else
   {
      if ((*consdata)->created)
      {
         for ( i = tcliqueGetNNodes((*consdata)->graph)-1; i >= 0; i-- )
         {
            if ( (*consdata)->nnodesinunion[i] > 0 )
            {
               SCIPfreeBlockMemoryArray(scip, &((*consdata)->unionofnode[i]), (*consdata)->nnodesinunion[i]);
               (*consdata)->unionofnode[i] = NULL;
            }
         }
         SCIPfreeBlockMemoryArray(scip, &((*consdata)->unionofnode), tcliqueGetNNodes((*consdata)->graph));
         SCIPfreeBlockMemoryArray(scip, &((*consdata)->nnodesinunion), tcliqueGetNNodes((*consdata)->graph));
         SCIPfreeBlockMemoryArray(scip, &((*consdata)->representativeofnode), tcliqueGetNNodes((*consdata)->graph));

         (*consdata)->unionofnode = NULL;
         (*consdata)->representativeofnode = NULL;
         (*consdata)->nnodesinunion = NULL;

         if ((*consdata)->graph != NULL)
         {
            tcliqueFree(&((*consdata)->graph));
         }
         if ((*consdata)->cgraph != NULL)
         {
            tcliqueFree(&((*consdata)->cgraph));
         }
      }
   }
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpStoreGraph)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsStoreGraph)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckStoreGraph)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockStoreGraph)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   SCIPdebugMessage("Locking method for store graph constraint: <%s>.\n", SCIPconsGetName(cons));

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveStoreGraph)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSDATA*     consdata;
   SCIP_CONSDATA*     olddata;
   TCLIQUE_GRAPH*     fathergraph;
   int   i;
   int   j;
   int*  firstedge;
   int*  lastedge;
   int   inserted;
   int   nnodes;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert((consdata->type == COLOR_CONSTYPE_ROOT) || (consdata->fathercons != NULL));

   SCIPdebugMessage("Activating store graph constraint: <%s(%d,%d)> [stack size: %d].\n", SCIPconsGetName(cons),
                         (consdata->node1+1), (consdata->node2+1), conshdlrData->nstack+1);

   /* put constraint on the stack */
   if ( conshdlrData->nstack >= conshdlrData->maxstacksize )
   {
      SCIPreallocMemoryArray(scip, &(conshdlrData->stack), 2*(conshdlrData->maxstacksize));
      conshdlrData->maxstacksize = 2*(conshdlrData->maxstacksize);
      SCIPdebugMessage("reallocating Memory for stack! %d --> %d\n", conshdlrData->maxstacksize/2, conshdlrData->maxstacksize);
   }
   conshdlrData->stack[conshdlrData->nstack] = cons;
   ++(conshdlrData->nstack);

   /* if the current graph was not yet created, create it now */
   if ( consdata->created == FALSE )
   {
      consdata->created = TRUE;
      olddata = SCIPconsGetData(consdata->fathercons);
      assert((consdata->type == COLOR_CONSTYPE_ROOT)
         || (consdata->node1 == olddata->representativeofnode[consdata->node1]
            && consdata->node2 == olddata->representativeofnode[consdata->node2]));
      nnodes = tcliqueGetNNodes(olddata->graph);
      fathergraph = olddata->graph;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->representativeofnode), nnodes) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->nnodesinunion), nnodes) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->unionofnode), nnodes) );

      for ( i = 0; i < nnodes; i++ )
      {
         consdata->representativeofnode[i] = olddata->representativeofnode[i];
         consdata->nnodesinunion[i] = olddata->nnodesinunion[i];
         if ( consdata->nnodesinunion[i] > 0 )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->unionofnode[i]), consdata->nnodesinunion[i]) );
            for ( j = 0; j < consdata->nnodesinunion[i]; j++ )
            {
               consdata->unionofnode[i][j] = olddata->unionofnode[i][j];
            }
         }
      }

      /* copy the graph */
      tcliqueCreate(&(consdata->graph));
      tcliqueAddNode((consdata)->graph, nnodes-1, 0);
      for ( i = 0; i < nnodes; i++ )
      {
         /* get adjacent nodes for node i and add them to new graph*/
         firstedge = tcliqueGetFirstAdjedge(fathergraph, i);
         lastedge = tcliqueGetLastAdjedge(fathergraph, i);
         while ( firstedge <= lastedge )
         {
            if ( *firstedge > i )
            {
               tcliqueAddEdge(consdata->graph, i, *firstedge);
            }
            firstedge++;
         }
      }
      tcliqueFlush(consdata->graph);
      assert(consdata->representativeofnode[consdata->node2] == consdata->node2);
      assert(consdata->representativeofnode[consdata->node1] == consdata->node1);

      /* type == COLOR_CONSTYPE_DIFFER --> insert edge between node1 and node2 */
      if (consdata->type == COLOR_CONSTYPE_DIFFER)
      {
         for ( i = 0; i < consdata->nnodesinunion[consdata->representativeofnode[consdata->node2]]; i++ )
         {
            for ( j = 0; j < consdata->nnodesinunion[consdata->representativeofnode[consdata->node1]]; j++ )
            {
               tcliqueAddEdge(consdata->graph, consdata->unionofnode[consdata->representativeofnode[consdata->node1]][j],
                  consdata->unionofnode[consdata->representativeofnode[consdata->node2]][i]);
            }
         }
         tcliqueFlush(consdata->graph);
      }
      /* type == COLOR_CONSTYPE_SAME --> insert edge (node2, i) - if not yet existing - if there exists an edge (node1, i) and vice versa */
      else
      {
         assert(consdata->type == COLOR_CONSTYPE_SAME);
         inserted = 0;

         /* add edges from all nodes of union2 to all nodes adjacent to union1 */
         for ( i = 0; i < consdata->nnodesinunion[consdata->node2]; i++ )
         {
            /* set representative of nodes in the union of node2 */
            consdata->representativeofnode[consdata->unionofnode[consdata->node2][i]] = consdata->node1;

            /* insert edges to all nodes adjacent to node1 */
            firstedge = tcliqueGetFirstAdjedge(fathergraph, consdata->node1);
            lastedge = tcliqueGetLastAdjedge(fathergraph, consdata->node1);
            while ( firstedge <= lastedge )
            {
               if ( !tcliqueIsEdge(fathergraph, *firstedge, consdata->node2) )
               {
                  inserted++;
                  tcliqueAddEdge(consdata->graph, consdata->unionofnode[consdata->node2][i], *firstedge);
               }
               firstedge++;
            }
         }
         /* add edges from all nodes of union1 to all nodes adjacent to union2 */
         for ( i = 0; i < consdata->nnodesinunion[consdata->node1]; i++ )
         {
            /* insert edges to all nodes adjacent to node2 */
            firstedge = tcliqueGetFirstAdjedge(fathergraph, consdata->node2);
            lastedge = tcliqueGetLastAdjedge(fathergraph, consdata->node2);
            while ( firstedge <= lastedge )
            {
               if ( !tcliqueIsEdge(fathergraph, *firstedge, consdata->node1) )
               {
                  inserted++;
                  tcliqueAddEdge(consdata->graph, consdata->unionofnode[consdata->node1][i], *firstedge);
               }
               firstedge++;
            }
         }
         if ( inserted > 0 )
         {
            tcliqueFlush(consdata->graph);
         }

         /* update union represented by node1 */
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->unionofnode[consdata->node1]),
               consdata->nnodesinunion[consdata->node1],
               (consdata->nnodesinunion[consdata->node1]) + (consdata->nnodesinunion[consdata->node2])) );
         for ( i = 0; i < consdata->nnodesinunion[consdata->node2]; i ++ )
         {
            consdata->unionofnode[consdata->node1][consdata->nnodesinunion[consdata->node1]+i]
               = consdata->unionofnode[consdata->node2][i];
         }
         SCIPfreeBlockMemoryArray(scip, &(consdata->unionofnode[consdata->node2]),
            consdata->nnodesinunion[consdata->node2]);
         consdata->nnodesinunion[consdata->node1] =
            (consdata->nnodesinunion[consdata->node1]) + (consdata->nnodesinunion[consdata->node2]);
         consdata->nnodesinunion[consdata->node2] = 0;
         consdata->unionofnode[consdata->node2] = NULL;
      }

      /* create the complementary graph */
      tcliqueCreate(&(consdata->cgraph));
      assert(consdata->cgraph != NULL);
      COLORprobGetComplementaryGraph(scip, consdata->graph, consdata->cgraph);
   }
   /* if new variables where created after the last propagation of this cons, repropagate it */
   else
   {
      if ( (consdata->type != COLOR_CONSTYPE_ROOT) && (consdata->propagatedvars < SCIPgetNTotalVars(scip)) )
      {
         SCIPrepropagateNode(scip, consdata->stickingatnode);
      }
   }

   return SCIP_OKAY;
}



/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveStoreGraph)
{
   SCIP_CONSHDLRDATA* conshdlrData;
#ifdef SCIP_DEBUG
   SCIP_CONSDATA*     consdata;
#endif

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);
   assert(cons == conshdlrData->stack[conshdlrData->nstack-1]);

#ifdef SCIP_DEBUG
   consdata = SCIPconsGetData(cons);

   SCIPdebugMessage("Deactivating store graph constraint: <%s(%d,%d)> [stack size: %d].\n", SCIPconsGetName(cons), (consdata->node1+1), (consdata->node2+1), conshdlrData->nstack-1);
#endif

   /* remove constraint from the stack */
   --conshdlrData->nstack;

   return SCIP_OKAY;
}



/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropStoreGraph)
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;
   SCIP_VAR*          var;
   int**              sets;
   int*               nsetelements;
   int                nsets;
   int                i;
   int                propcount;

   assert(conshdlr != NULL);
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   /* get all stable sets */
   COLORprobGetStableSets(scip, &sets, &nsetelements, &nsets);
   *result = SCIP_DIDNOTFIND;
   propcount = 0;

   /* the constraint data of the cons related to the current node */
   cons = conshdlrData->stack[conshdlrData->nstack-1];
   consdata = SCIPconsGetData(cons);

   SCIPdebugMessage( "Starting propagation of store graph constraint <%s(%d,%d)> .\n", SCIPconsGetName(cons), (consdata->node1+1), (consdata->node2+1));

   /* propagation for differ: set upper bound to 0 for all stable sets, which contain both nodes */
   if (consdata->type == COLOR_CONSTYPE_DIFFER)
   {
      for ( i = 0; i < nsets; i++ )
      {
         if ( !SCIPisFeasZero(scip, SCIPvarGetUbLocal(COLORprobGetVarForStableSet(scip, i))) )
         {
            if ( COLORprobIsNodeInStableSet(scip, i, consdata->node1) && COLORprobIsNodeInStableSet(scip, i, consdata->node2) )
            {
               var = COLORprobGetVarForStableSet(scip, i);
               SCIP_CALL( SCIPchgVarUb(scip, var, 0) );
               propcount++;
            }
         }
      }
   }

   /* propagation for same: set upper bound to 0 for all stable sets, which do not contain both nodes */
   if ( consdata->type == COLOR_CONSTYPE_SAME )
   {
      for ( i = 0; i < nsets; i++ )
      {
         if ( !SCIPisFeasZero(scip, SCIPvarGetUbLocal(COLORprobGetVarForStableSet(scip, i))) )
         {
            if ( (COLORprobIsNodeInStableSet(scip, i, consdata->node1) || COLORprobIsNodeInStableSet(scip, i, consdata->node2))
               && !(COLORprobIsNodeInStableSet(scip, i, consdata->node1) && COLORprobIsNodeInStableSet(scip, i, consdata->node2)) )
            {
               var = COLORprobGetVarForStableSet(scip, i);
               SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
               propcount++;
            }
         }
      }
   }

   SCIPdebugMessage( "Finished propagation of store graph constraint <%s(%d,%d)>, %d vars fixed.\n", SCIPconsGetName(cons), (consdata->node1+1), (consdata->node2+1), propcount);

   consdata = SCIPconsGetData(COLORconsGetActiveStoreGraphCons(scip));
   consdata->propagatedvars = SCIPgetNTotalVars(scip);

   return SCIP_OKAY;
}

/*
 * interface methods
 */


/** creates the handler for storeGraph constraints and includes it in SCIP */
SCIP_RETCODE COLORincludeConshdlrStoreGraph(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONSHDLR* conshdlr;

   SCIPdebugMessage("Including graph storage constraint handler.\n");

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrData) );
   conshdlrData->stack = NULL;
   conshdlrData->nstack = 0;
   conshdlrData->maxstacksize = 25;

   conshdlr = NULL;
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpStoreGraph, consEnfopsStoreGraph, consCheckStoreGraph, consLockStoreGraph,
         conshdlrData) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteStoreGraph) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeStoreGraph) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolStoreGraph) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolStoreGraph) );
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveStoreGraph) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveStoreGraph) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropStoreGraph, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );

   return SCIP_OKAY;
}

/** creates and captures a storeGraph constraint, uses knowledge of the B&B-father*/
SCIP_RETCODE COLORcreateConsStoreGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONS*            fatherconstraint,   /**< constraint in B&B-father */
   int                   type,               /**< type of the constraint: COLOR_CONSTYPE_SAME or COLOR_CONSTYPE_DIFFER */
   int                   node1,              /**< the first node of the constraint */
   int                   node2,              /**< the second node of the constraint */
   SCIP_NODE*            stickingnode        /**< the B&B-tree node at which the constraint will be sticking */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int   temp;

   assert(scip != NULL);
   assert(fatherconstraint != NULL);
   assert(type == COLOR_CONSTYPE_SAME || type == COLOR_CONSTYPE_DIFFER);
   assert(stickingnode != NULL);

   /* find the storeGraph constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   if ( node1 > node2 )
   {
      temp = node1;
      node1 = node2;
      node2 = temp;
   }
   SCIPdebugMessage("Creating store graph constraint: <%s(%d,%d)>. \n", name, (node1+1), (node2+1));

   consdata->node1 = node1;
   consdata->node2 = node2;
   consdata->type = type;
   consdata->fathercons = fatherconstraint;
   consdata->propagatedvars = 0;
   consdata->stickingatnode = stickingnode;
   consdata->created = FALSE;


   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         TRUE, FALSE, TRUE, FALSE, TRUE) );

   return SCIP_OKAY;
}




/* ----------------------------------- external methods -------------------------- */

/** returns the store graph constraint of the current node, needs the pointer to the constraint handler */
SCIP_CONS* COLORconsGetActiveStoreGraphConsFromHandler(
   SCIP_CONSHDLR*        conshdlr            /**< constaint handler for store-graph constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(conshdlr != NULL);
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the store graph constraint of the current node, only needs the pointer to scip */
SCIP_CONS* COLORconsGetActiveStoreGraphCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return NULL;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   assert(conshdlrData->nstack > 0);

   return conshdlrData->stack[conshdlrData->nstack-1];
}


/** returns the current graph */
TCLIQUE_GRAPH* COLORconsGetCurrentGraph(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return NULL;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   cons = conshdlrData->stack[conshdlrData->nstack-1];
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   return consdata->graph;
}


/** returns the complementary graph */
TCLIQUE_GRAPH* COLORconsGetComplementaryGraph(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return NULL;
   }

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   cons = conshdlrData->stack[conshdlrData->nstack-1];
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   return consdata->cgraph;
}


/** returns array of representatives of all nodes */
int* COLORconsGetRepresentatives(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return NULL;
   }

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   cons = conshdlrData->stack[conshdlrData->nstack-1];
   consdata = SCIPconsGetData(cons);
   return consdata->representativeofnode;
}

/** returns the representative of the union which contains a given node */
int COLORconsGetRepresentative(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< the node, for wich the representative is searched */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return -1;
   }

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   cons = conshdlrData->stack[conshdlrData->nstack-1];
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(node >= 0 && node < tcliqueGetNNodes(consdata->graph));

   return consdata->representativeofnode[node];
}

/** returns the array of all unions, a union is saved in the array at the position of its representative */
void COLORconsGetUnions(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                unions,             /**< output: array containing array which contains nodes in the union */
   int**                 lengths             /**< output: lengths of the unions */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return;
   }

   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   cons = conshdlrData->stack[conshdlrData->nstack-1];
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *unions = consdata->unionofnode;
   *lengths = consdata->nnodesinunion;
}

/** returns the union which has a given node as representative */
void COLORconsGetUnion(
   SCIP*                 scip,               /**< SCIP data structure */
   int**                 nodesinunion,       /**< output: array containig nodes in the union */
   int*                  nnodesinunion,      /**< output: length of the union */
   int                   node                /**< the node, whose union we want to get */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;
   SCIP_CONS*         cons;
   SCIP_CONSDATA*     consdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);
   cons = conshdlrData->stack[conshdlrData->nstack-1];
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *nodesinunion = consdata->unionofnode[node];
   *nnodesinunion = consdata->nnodesinunion[node];
}

/** returns the stack and the number of elements on it */
void COLORconsGetStack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          stack,              /**< return value: pointer to the stack */
   int*                  nstackelements      /**< return value: pointer to int, for number of elements on the stack */
   )
{
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrData;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, "storeGraph");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("storeGraph constraint handler not found\n");
      return;
   }
   conshdlrData = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrData != NULL);
   assert(conshdlrData != NULL);
   assert(conshdlrData->stack != NULL);

   *stack = conshdlrData->stack;
   *nstackelements = conshdlrData->nstack;
}


