/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_pseudo.c
 * @brief  pseudo branching rule
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>
#include "scip/branch_pseudo.h"
#include "scip/branch_nodereopt.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/scip.h"
#include "scip/tree.h"

#define BRANCHRULE_NAME            "pseudo"
#define BRANCHRULE_DESC            "branching rule for pseudo-branched nodes in reoptimization"
#define BRANCHRULE_PRIORITY        -536870911
#define BRANCHRULE_MAXDEPTH        0
#define BRANCHRULE_MAXBOUNDDIST    1.0

/*
 * Data structures
 */

struct LogicOrData
{
   SCIP_VAR**            vars;
   SCIP_Real*            vals;
   REOPT_CONSTYPE        constype;
   int                   allocmem;
   int                   nvars;
};

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_Bool             init;                    /** is data initialized? */

   /** data of pseudo-branched variables */
   LOGICORDATA**         consdata;
   int*                  npseudobranchvars;       /** list of numbers of pseudo-branched variables. */
   int                   nrpseudobranchednodes;   /** number of already saved pseudo-branched nodes in this pricing round. */

   /** data to handle free IDs */
   SCIP_QUEUE*           openIDs;                 /** queue of empty IDs in node-data array */
   int*                  exIDtoinID;              /** array which links an ID from branch_nodereopt to an ID from branch_pseudo. */
   int*                  nodetoid;                /** array of all IDs corresponding to saved nodes; we need this to find the IDs of the
                                                      parent node. key: node number */
   int                   nodeID;                  /** the current nodeID, -1 if no current node exist */
   int                   allocmemsizenodeID;      /** size of allocated memory for pseudo-branched nodes */
   int                   allocmemsizeIDID;        /** size of allocated memory for IDs and linked node numbers */
   int                   allocmemsizeconsdata;    /** size of allocated memory for consdata */

   /** statistic stuff */
   SCIP_Longint          lastseennode;            /** number of the last seen node */
   SCIP_Bool             newnode;

   SCIP_Bool             reopt;                   /** is reoptimization enabled? */
};

/*
 * allocate memory for consdata and initialize all node parameters
 */
static
SCIP_RETCODE initCons(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   nodeID
)
{
   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(branchruledata->reopt);
   assert(nodeID >= 0);
   assert(branchruledata->consdata[nodeID] == NULL);

   SCIP_CALL( SCIPallocMemory(scip, &branchruledata->consdata[nodeID]) );
   branchruledata->consdata[nodeID]->allocmem = 0;
   branchruledata->consdata[nodeID]->nvars = 0;
   branchruledata->consdata[nodeID]->vars = NULL;
   branchruledata->consdata[nodeID]->vals = NULL;

   return SCIP_OKAY;
}

/**
 * delete the data for node nodeID in consdata
 */
static
SCIP_RETCODE deleteConsData(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   internID,
   SCIP_Bool             exitsolving
)
{
   assert(scip != NULL);
   assert(branchruledata != NULL);
   assert(branchruledata->consdata[internID] != NULL);

   if( exitsolving )
   {
      if( branchruledata->consdata[internID]->vars != NULL )
      {
        SCIPfreeMemoryArray(scip, &branchruledata->consdata[internID]->vars);
      }

      if( branchruledata->consdata[internID]->vals != NULL )
      {
        SCIPfreeMemoryArray(scip, &branchruledata->consdata[internID]->vals);
      }

      SCIPfreeMemory(scip, &branchruledata->consdata[internID] );
      branchruledata->consdata[internID] = NULL;
   }
   else
      branchruledata->consdata[internID]->nvars = 0;

   return SCIP_OKAY;
}

static
SCIP_RETCODE reallocNodedata(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata
)
{
   int oldsize;
   int pos;

   assert( scip != NULL );
   assert( branchruledata != NULL );

   oldsize = branchruledata->allocmemsizeconsdata;
   branchruledata->allocmemsizeconsdata *= 2;

   SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->consdata), branchruledata->allocmemsizeconsdata) );

   /* write 0 empty slots and fill the queue*/
   for(pos = branchruledata->allocmemsizeconsdata-1; pos >= oldsize; pos--)
   {
      SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) pos) );
      branchruledata->consdata[pos] = NULL;
   }
   assert( SCIPqueueNElems(branchruledata->openIDs) == branchruledata->allocmemsizeconsdata/2 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE reallocIDID(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   int                   externID
)
{
   int oldsize;
   int pos;

   assert( scip != NULL );
   assert( branchruledata != NULL );

   oldsize = branchruledata->allocmemsizeIDID;
   branchruledata->allocmemsizeIDID = 2*externID;

   SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->exIDtoinID), branchruledata->allocmemsizeIDID) );

   /* write 0 empty slots and fill the queue*/
   for(pos = branchruledata->allocmemsizeIDID-1; pos >= oldsize; pos--)
      branchruledata->exIDtoinID[pos] = 0;

   return SCIP_OKAY;
}

static
SCIP_RETCODE reallocNodeID(
   SCIP*                 scip,
   SCIP_BRANCHRULEDATA*  branchruledata,
   SCIP_Longint          nodenumber
)
{
   int oldsize;
   int pos;

   assert( scip != NULL );
   assert( branchruledata != NULL );

   oldsize = branchruledata->allocmemsizenodeID;
   branchruledata->allocmemsizenodeID = 2*nodenumber;

   SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->nodetoid), branchruledata->allocmemsizenodeID) );

   /* write 0 empty slots and fill the queue*/
   for(pos = branchruledata->allocmemsizenodeID-1; pos >= oldsize; pos--)
      branchruledata->nodetoid[pos] = 0;

   return SCIP_OKAY;
}

/*
 * Local methods
 */
SCIP_RETCODE SCIPbranchrulePseudoGenerateCons(
   SCIP*                 scip,
   LOGICORDATA*          consdata,
   int*                  nvars,
   int                   nallocvars,
   int                   externID,
   SCIP_Bool             local,
   SCIP_Bool             cleardata
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int nodeID;
   int varnr;

   assert( scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL);
   assert(branchruledata->reopt);
   assert(externID >= 0);
   assert(consdata != NULL);

   if( externID == 0 )
      nodeID = externID;
   else
   {
      assert( branchruledata->exIDtoinID[externID] != 0 );
      nodeID = branchruledata->exIDtoinID[externID];
   }
   assert(branchruledata->consdata[nodeID] != NULL);
   assert(branchruledata->consdata[nodeID]->vars != NULL);
   assert(branchruledata->consdata[nodeID]->nvars > 0);
   assert(nodeID >= 0);

   (*nvars) = branchruledata->consdata[nodeID]->nvars;
   consdata->constype = branchruledata->consdata[nodeID]->constype;

   if( nallocvars < branchruledata->consdata[nodeID]->nvars )
      return SCIP_OKAY;

   /** copy logic-or data */
   for(varnr = 0; varnr < (*nvars); varnr++)
   {
      consdata->vars[varnr] = branchruledata->consdata[nodeID]->vars[varnr];
      consdata->vals[varnr] = branchruledata->consdata[nodeID]->vals[varnr];
   }

   if( cleardata )
   {
      SCIP_CALL( deleteConsData(scip, branchruledata, nodeID, FALSE) );

      if( nodeID >= 1 )
      {
         branchruledata->exIDtoinID[externID] = 0;
         SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) nodeID) );
      }
   }

   return SCIP_OKAY;
}

/*
 * link an extern ID to an intern ID
 */
SCIP_RETCODE SCIPbranchrulePseudoLinkIDs(
   SCIP*                 scip,
   int                   externID
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert( scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert( branchruledata != NULL );
   assert( branchruledata->nodeID > 0 ); /* we don't link the root IDs, because the root gets always ID 0 */

   if( externID >= branchruledata->allocmemsizeIDID )
   {
      SCIP_CALL( reallocIDID(scip, branchruledata, externID) );
      assert(branchruledata->exIDtoinID[externID] == 0);
   }

   branchruledata->exIDtoinID[externID] = branchruledata->nodeID;

   return SCIP_OKAY;
}

/*
 * return if a given node is pseudo-branched or not
 */
SCIP_Bool SCIPbranchrulePseudoIsPseudoBranched(
   SCIP*             scip,
   SCIP_NODE*        node
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert( scip != NULL );
   assert( node != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert( branchruledata != NULL );

   if( node == SCIPgetRootNode(scip)
    && branchruledata->consdata[0] != NULL
    && branchruledata->consdata[0]->vars != NULL
    && branchruledata->consdata[0]->nvars > 0
    && branchruledata->nodeID == 0 )
      return TRUE;

   else if( node != SCIPgetRootNode(scip) && SCIPnodeGetNumber(node) < branchruledata->allocmemsizenodeID )
   {
      if( branchruledata->nodetoid[SCIPnodeGetNumber(node)-1] != 0
       && branchruledata->consdata[branchruledata->nodetoid[SCIPnodeGetNumber(node)-1]]->nvars > 0 )
      return TRUE;
   }
   return FALSE;
}

/*
 * add a variable to consdata
 */
SCIP_RETCODE SCIPbranchrulePseudoAddPseudoVar(
      SCIP*             scip,
      SCIP_NODE*        node,
      SCIP_VAR*         var,
      SCIP_BOUNDTYPE    boundtype,
      int               newbound
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Real scalar;
   SCIP_Real constant;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert( branchruledata != NULL );

   if( branchruledata->newnode )
   {
      assert( branchruledata->lastseennode != SCIPnodeGetNumber(node) );
      assert( branchruledata->nodeID == -1 );

      /** reallocate memory */
      if( SCIPqueueIsEmpty(branchruledata->openIDs) )
      {
         SCIP_CALL( reallocNodedata(scip, branchruledata) );
      }
      assert( branchruledata->nrpseudobranchednodes < branchruledata->allocmemsizeconsdata );

      /** get an empty slot (to ensure that ID != NULL, the queue should contains only IDs >= 1); the 0 is always reserved for the root */
      if( SCIPgetFocusDepth(scip) == 0 )
         branchruledata->nodeID = 0;
      else
      {
         branchruledata->nodeID = (int) (size_t) SCIPqueueRemove(branchruledata->openIDs);
         assert( branchruledata->nodeID >= 1 && branchruledata->nodeID <= branchruledata->allocmemsizeconsdata );
      }

      /**
       * if var == NULL, we save all information in the next step in NodeFinished().
       */
      if( var != NULL )
      {
         assert(newbound > -1);
         assert(SCIPgetFocusDepth(scip) <= SCIPgetEffectiveRootDepth(scip) || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE);

         if( branchruledata->consdata[branchruledata->nodeID] == NULL )
         {
            /** number of pseudo-branched variables could be much bigger than the depth,
             * thats way use #variable-depth as default value */
            SCIP_CALL( initCons(scip, branchruledata, branchruledata->nodeID) );
            branchruledata->consdata[branchruledata->nodeID]->allocmem = SCIPgetNBinVars(scip);
            SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->consdata[branchruledata->nodeID]->vars), branchruledata->consdata[branchruledata->nodeID]->allocmem) );
            SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->consdata[branchruledata->nodeID]->vals), branchruledata->consdata[branchruledata->nodeID]->allocmem) );
         }
         assert(branchruledata->consdata[branchruledata->nodeID]->nvars == 0);

         scalar = 1;
         constant = 0;

         /*
          * Retransform the bound and ensure that the variable is original.
          */
         if( !SCIPvarIsOriginal(var) )
         {
            SCIP_CALL(SCIPvarGetOrigvarSum(&var, &scalar, &constant));
            newbound = (newbound - constant) / scalar;
         }
         assert(SCIPvarIsOriginal(var));

         branchruledata->consdata[branchruledata->nodeID]->vars[0] = var;
         branchruledata->consdata[branchruledata->nodeID]->vals[0] = newbound;
         branchruledata->consdata[branchruledata->nodeID]->nvars = 1;

         assert( SCIPisFeasEQ(scip, branchruledata->consdata[branchruledata->nodeID]->vals[0], 0 )
              || SCIPisFeasEQ(scip, branchruledata->consdata[branchruledata->nodeID]->vals[0], 1) );
      }

      branchruledata->newnode = FALSE;
      branchruledata->lastseennode = SCIPnodeGetNumber(node);

      if( branchruledata->nodeID > 0 ) /* ID 0 would be the root */
      {
         if( SCIPnodeGetNumber(node) >= branchruledata->allocmemsizenodeID )
         {
            SCIP_CALL( reallocNodeID(scip, branchruledata, SCIPnodeGetNumber(node)) );
         }
         assert( branchruledata->nodetoid[SCIPnodeGetNumber(node)-1] == 0
              || branchruledata->consdata[branchruledata->nodetoid[SCIPnodeGetNumber(node)-1]]->nvars == 0);
         branchruledata->nodetoid[SCIPnodeGetNumber(node)-1] = branchruledata->nodeID;
      }
   }
   else
   {
      int pbvars;

      assert(var != NULL);
      assert(newbound != -1);
      assert( SCIPisFeasEQ(scip, newbound, 0) || SCIPisFeasEQ(scip, newbound, 1) );
      assert(branchruledata->consdata[branchruledata->nodeID] != NULL);
      assert(branchruledata->consdata[branchruledata->nodeID]->allocmem > 0);
      assert(branchruledata->consdata[branchruledata->nodeID]->nvars >= 1);


      pbvars = branchruledata->consdata[branchruledata->nodeID]->nvars;

      /**
       * reallocate if the list is completely filled
       */
      if(branchruledata->consdata[branchruledata->nodeID]->allocmem == pbvars)
      {
         assert(FALSE); /* this should never be the case */
         branchruledata->consdata[branchruledata->nodeID]->allocmem *= 2;

         SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->consdata[branchruledata->nodeID]->vars), branchruledata->consdata[branchruledata->nodeID]->allocmem));
         SCIP_CALL( SCIPreallocMemoryArray(scip, &(branchruledata->consdata[branchruledata->nodeID]->vals), branchruledata->consdata[branchruledata->nodeID]->allocmem));
      }

      assert( branchruledata->consdata[branchruledata->nodeID]->vars != NULL );
      assert( branchruledata->consdata[branchruledata->nodeID]->vals != NULL );

      scalar = 1;
      constant = 0;

      /*
       * Retransform the bound and ensure that the variable is original.
       */
      if( !SCIPvarIsOriginal(var) )
      {
         SCIP_CALL(SCIPvarGetOrigvarSum(&var, &scalar, &constant));
         newbound = (newbound - constant) / scalar;
      }
      assert(SCIPvarIsOriginal(var));

      branchruledata->consdata[branchruledata->nodeID]->vars[pbvars] = var;
      branchruledata->consdata[branchruledata->nodeID]->vals[pbvars] = newbound;
      branchruledata->consdata[branchruledata->nodeID]->nvars++;

      assert( SCIPisFeasEQ(scip, branchruledata->consdata[branchruledata->nodeID]->vals[pbvars], 0)
           || SCIPisFeasEQ(scip, branchruledata->consdata[branchruledata->nodeID]->vals[pbvars], 1) );
   }

   return SCIP_OKAY;
}

/*
 * finish the collecting of information
 */
SCIP_RETCODE SCIPbranchrulePseudoNodeFinished(
   SCIP*                 scip,
   SCIP_NODE*            node,
   REOPT_CONSTYPE        constype
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_Real constant;
   SCIP_Real scalar;
   int nbranchvars;
   int npseudobranchvars;
   int var;

   assert( scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert( branchruledata != NULL );

   /**
    * we know the exact number of pseudo-branchings
    * all boundchanges are local an we can look at the domchg data-structure
    *
    * if the node depth is equal to the effective root depth all boundchanges are already known.
    */
   if( SCIPgetFocusDepth(scip) > SCIPgetEffectiveRootDepth(scip) )
   {
      npseudobranchvars = SCIPgetNBinVars(scip);
      assert( npseudobranchvars > 0 );

      if( branchruledata->consdata[branchruledata->nodeID] == NULL )
      {
         /** allocate memory */
         SCIP_CALL( initCons(scip, branchruledata, branchruledata->nodeID) );
         branchruledata->consdata[branchruledata->nodeID]->allocmem = npseudobranchvars;
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->consdata[branchruledata->nodeID]->vars), branchruledata->consdata[branchruledata->nodeID]->allocmem) );
         SCIP_CALL( SCIPallocClearMemoryArray(scip, &(branchruledata->consdata[branchruledata->nodeID]->vals), branchruledata->consdata[branchruledata->nodeID]->allocmem) );
      }
      assert(branchruledata->consdata[branchruledata->nodeID]->nvars == 0);

      /** save all boundchanges which don't correspond to branching decisions */
      SCIPnodeGetPseudoBranchings(node,
                                 branchruledata->consdata[branchruledata->nodeID]->vars,
                                 branchruledata->consdata[branchruledata->nodeID]->vals,
                                 &nbranchvars,
                                 npseudobranchvars);

      assert(nbranchvars <= branchruledata->consdata[branchruledata->nodeID]->allocmem);
      branchruledata->consdata[branchruledata->nodeID]->nvars = nbranchvars;

      /*
       * Retransform the bound and ensure that the variable is original.
       */
      for(var = 0; var < branchruledata->consdata[branchruledata->nodeID]->nvars; var++)
      {
         constant = 0;
         scalar = 1;

         if( !SCIPvarIsOriginal(branchruledata->consdata[branchruledata->nodeID]->vars[var]) )
         {
            SCIP_CALL(SCIPvarGetOrigvarSum(&(branchruledata->consdata[branchruledata->nodeID]->vars[var]), &scalar, &constant));
            branchruledata->consdata[branchruledata->nodeID]->vals[var] = (branchruledata->consdata[branchruledata->nodeID]->vals[var] - constant) / scalar;
         }
         assert( SCIPvarIsOriginal(branchruledata->consdata[branchruledata->nodeID]->vars[var]) );
         assert( 0 <= branchruledata->consdata[branchruledata->nodeID]->vals[var] && branchruledata->consdata[branchruledata->nodeID]->vals[var] <= 1 );
      }
   }

   /* set REOPT_CONSTYPE */
   branchruledata->consdata[branchruledata->nodeID]->constype = constype;

   assert( branchruledata->consdata[branchruledata->nodeID]->nvars > 0 );
   assert( branchruledata->consdata[branchruledata->nodeID]->allocmem > 0 );
   assert( branchruledata->consdata[branchruledata->nodeID]->vars != NULL );
   assert( branchruledata->consdata[branchruledata->nodeID]->vals != NULL );

#ifdef SCIP_DEBUG
   {
      int varnr;
      SCIPdebugMessage("finish collecting strong branching information about node %lld.\n", SCIPnodeGetNumber(node));
      SCIPdebugMessage(" -> nvars: %d\n", branchruledata->consdata[branchruledata->nodeID]->nvars);
      for(varnr = 0; varnr < branchruledata->consdata[branchruledata->nodeID]->nvars; ++varnr)
      {
         SCIPdebugMessage("  %s %s %f\n", SCIPvarGetName( branchruledata->consdata[branchruledata->nodeID]->vars[varnr] ), branchruledata->consdata[branchruledata->nodeID]->vals[varnr] == 1 ? "=>" : "<=", branchruledata->consdata[branchruledata->nodeID]->vals[varnr] );
      }
   }
#endif

   branchruledata->nodetoid[SCIPnodeGetNumber(node)-1] = 0;
   branchruledata->nrpseudobranchednodes++;
   branchruledata->newnode = TRUE;
   branchruledata->nodeID = -1;

   return SCIP_OKAY;
}

/*
 * delete branching information for a given nodeID
 */
SCIP_RETCODE SCIPbranchrulePseudoDelInformation(
   SCIP*                 scip,
   int                   externID
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int internID;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL);
   assert(externID >= 0);

   /* get the liked ID*/
   if( externID == 0 )
      internID = 0;
   else
   {
      assert( 1 <= branchruledata->exIDtoinID[externID] && branchruledata->exIDtoinID[externID] <= branchruledata->allocmemsizeconsdata );
      internID = branchruledata->exIDtoinID[externID];
   }

   SCIP_CALL( deleteConsData(scip, branchruledata, internID, FALSE) );

   if( externID >= 1 )
   {
      branchruledata->exIDtoinID[externID] = 0;
      SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) internID) );
   }

   return SCIP_OKAY;
}

/*
 * reset the links between node numbers and IDs
 */
SCIP_RETCODE SCIPbranchrulePseudoReset(
   SCIP*                 scip,
   SCIP_Bool             restart,
   SCIP_Bool             exceptroot
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int i;

   assert( scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL);
   assert(branchruledata->init);
   assert(branchruledata->reopt);

   BMSclearMemoryArray(branchruledata->nodetoid, branchruledata->allocmemsizenodeID);

   if( restart )
   {
      /* clear externID to internID links */
      BMSclearMemoryArray(branchruledata->exIDtoinID, branchruledata->allocmemsizeIDID);

      /** clear queue with open IDs */
      SCIPqueueClear(branchruledata->openIDs);

      for(i = (int) exceptroot; i < branchruledata->allocmemsizeconsdata; i++)
      {
         if( branchruledata->consdata[i] != NULL )
         {
            SCIP_CALL( deleteConsData(scip, branchruledata, i, FALSE) );
         }

         if( i > 0 )
         {
            SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) i) );
         }
      }
   }

   assert(!restart || SCIPqueueNElems(branchruledata->openIDs) == branchruledata->allocmemsizeconsdata-1);

   branchruledata->lastseennode = -1;
   branchruledata->newnode = TRUE;
   if( exceptroot )
      branchruledata->nrpseudobranchednodes = ( branchruledata->consdata[0] == NULL ? 0 : branchruledata->consdata[0]->nvars > 0 ? 1 : 0 );
   else
      branchruledata->nrpseudobranchednodes = 0;

   return SCIP_OKAY;
}

/*
 * delete information of the last seen node
 */
SCIP_RETCODE SCIPbranchrulePseudoDeleteLastNodeInfo(
      SCIP*          scip,
      SCIP_NODE*     node
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int internID;

   assert( scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert( branchruledata != NULL );
   assert( branchruledata->nodeID >= 0 );

   if( SCIPnodeGetNumber(node) != branchruledata->lastseennode )
      return SCIP_OKAY;

   if( SCIPgetRootNode(scip) != node )
   {
      assert( branchruledata->nodetoid[SCIPnodeGetNumber(node)-1] != 0 );
      assert( branchruledata->nodetoid[SCIPnodeGetNumber(node)-1] == branchruledata->nodeID );
      internID = branchruledata->nodeID;
      branchruledata->nodetoid[SCIPnodeGetNumber(node)-1] = 0;
   }
   else
   {
      assert( branchruledata->nodetoid[0] == 0 );
      assert( branchruledata->nodeID == 0);
      internID = 0;
   }

   SCIPdebugMessage("delete %d dual variable information about node %lld\n", branchruledata->consdata[internID]->nvars, SCIPnodeGetNumber(node));

   SCIP_CALL( deleteConsData(scip, branchruledata, internID, FALSE) );

   branchruledata->newnode = TRUE;
   branchruledata->lastseennode = -1;
   branchruledata->nodeID = -1;


   return SCIP_OKAY;
}

/*
 * return the number of saved pseudo-branched nodes
 */
int SCIPbranchrulePseudoGetNPseudoNodes(
      SCIP*             scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert( scip != NULL );

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert( branchrule != NULL );

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert( branchruledata != NULL );

   return branchruledata->nrpseudobranchednodes;

}

int SCIPbranchrulePseudoGetNPseudoVars(
   SCIP*                 scip,
   int                   externID
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;
   int internID;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);

   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);

   assert(branchruledata != NULL);
   assert(externID < branchruledata->allocmemsizeIDID);

   if( externID == 0 )
      internID = 0;
   else if( branchruledata->exIDtoinID[externID] != 0 )
      internID = branchruledata->exIDtoinID[externID];
   else
      return 0;

   if( branchruledata->consdata[internID] != NULL )
      return branchruledata->consdata[internID]->nvars;
   else
      return 0;
}

void SCIPbranchrulePseudoDisable(
   SCIP*                 scip
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   assert(scip != NULL);

   branchrule = SCIPfindBranchrule(scip, BRANCHRULE_NAME);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   branchruledata->reopt = FALSE;

   return;
}

/** copy method for branchrule plugins (called when SCIP copies plugins) */
#define branchCopyPseudo NULL;

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreePseudo)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   int i;

   assert(scip != NULL);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   if( branchruledata->init )
   {
      for(i = 0; i < branchruledata->allocmemsizeconsdata; ++i)
      {
         if( branchruledata->consdata[i] != NULL )
         {
            SCIP_CALL( deleteConsData(scip, branchruledata, i, TRUE) );
         }
      }

      SCIPfreeMemoryArray(scip, &branchruledata->consdata);
      SCIPqueueFree(&branchruledata->openIDs);
      SCIPfreeMemoryArray(scip, &branchruledata->exIDtoinID);
      SCIPfreeMemoryArray(scip, &branchruledata->nodetoid);

      branchruledata->init = FALSE;
      branchruledata->reopt = FALSE;
   }
   assert(!branchruledata->init);

   SCIPfreeMemory(scip, &branchruledata);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitPseudo)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /** HACK */
   /** check if all variable are binary, if not, disable reoptimization */
   if( !branchruledata->init )
   {
      int maxsavednodes;

      SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/enable", &branchruledata->reopt) );
      if( branchruledata->reopt && SCIPgetNOrigImplVars(scip) + SCIPgetNOrigIntVars(scip) > 0 )
         branchruledata->reopt = FALSE;

      SCIP_CALL( SCIPgetIntParam(scip, "reoptimization/maxsavednodes", &maxsavednodes) );

      if( maxsavednodes == 0 )
         branchruledata->reopt = FALSE;
   }

   /** initialize the data and change parameters */
   if( !branchruledata->init && branchruledata->reopt )
   {
      int slot;

      branchruledata->allocmemsizenodeID = 1000;
      branchruledata->allocmemsizeIDID = 1000;
      branchruledata->allocmemsizeconsdata = 1000;
      branchruledata->nodeID = -1;
      branchruledata->lastseennode = -1;
      branchruledata->newnode = TRUE;
      branchruledata->nrpseudobranchednodes = 0;

      SCIP_CALL(SCIPallocClearMemoryArray(scip, &(branchruledata->consdata), branchruledata->allocmemsizeconsdata));
      SCIP_CALL(SCIPallocClearMemoryArray(scip, &(branchruledata->nodetoid), branchruledata->allocmemsizenodeID));
      SCIP_CALL(SCIPallocClearMemoryArray(scip, &(branchruledata->exIDtoinID), branchruledata->allocmemsizeIDID));

      /** data structure to handle IDs */
      SCIP_CALL( SCIPqueueCreate(&(branchruledata->openIDs), branchruledata->allocmemsizenodeID-1, 2) );

      /** fill the queue with free IDs 1,...,allocnodeIDs-1; the ID 0 is always reserved for the root node */
      slot = 1;
      while( slot < branchruledata->allocmemsizenodeID )
      {
         SCIP_CALL( SCIPqueueInsert(branchruledata->openIDs, (void*) (size_t) slot) );
         slot++;
      }
      assert( SCIPqueueNElems(branchruledata->openIDs) == branchruledata->allocmemsizenodeID-1 );

      branchruledata->init = TRUE;
   }
   return SCIP_OKAY;
}

/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitPseudo NULL;

/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolPseudo NULL;

/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolPseudo)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   branchruledata->lastseennode = -1;
   branchruledata->newnode = TRUE;
   branchruledata->nodeID = -1;

   return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpPseudo)
{
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}
/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextPseudo)
{
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsPseudo)
{
   *result = SCIP_DIDNOTRUN;
   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the pseudo branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchrulePseudo(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULE* branchrule;
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create pseudo branching rule data */
   branchruledata = NULL;

   SCIP_CALL(SCIPallocMemory(scip, &branchruledata));

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);
   branchruledata->init = FALSE;

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetBranchruleFree(scip, branchrule, branchFreePseudo));
   SCIP_CALL(SCIPsetBranchruleInit(scip, branchrule, branchInitPseudo));
   SCIP_CALL(SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpPseudo));
   SCIP_CALL(SCIPsetBranchruleExecExt(scip, branchrule, branchExecextPseudo));
   SCIP_CALL(SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsPseudo));
   SCIP_CALL(SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolPseudo));

   return SCIP_OKAY;
}
