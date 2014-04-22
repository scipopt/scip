/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   visual.c
 * @brief  methods for creating output for visualization tools (VBC, BAK)
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 *
 * Output can be generated for the following visualization tools:
 *
 * - VBCTOOL - a graphical interface for Visualization of Branch Cut algorithms @n
 *   See <a href="http://www.informatik.uni-koeln.de/ls_juenger/research/vbctool">VBCTOOL</a>.
 * - BAK: Branch-and-bound Analysis Kit @n
 *   BAK is available through COIN-OR, see <a href="https://projects.coin-or.org/CoinBazaar/wiki/Projects/BAK">BAK</a>.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "blockmemshell/memory.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/var.h"
#include "scip/tree.h"
#include "scip/visual.h"
#include "scip/struct_visual.h"


/** returns the branching variable of the node, or NULL */
static
void getBranchInfo(
   SCIP_NODE*            node,               /**< node */
   SCIP_VAR**            var,                /**< pointer to store the branching variable */
   SCIP_BOUNDTYPE*       boundtype,          /**< pointer to store the branching type: lower or upper bound */
   SCIP_Real*            bound               /**< pointer to store the new bound of the branching variable */
   )
{
   SCIP_DOMCHGBOUND* domchgbound;

   (*var) = NULL;
   (*bound) = 0.0;
   (*boundtype) = SCIP_BOUNDTYPE_LOWER;

   assert(node != NULL);
   if( node->domchg == NULL )
      return;

   domchgbound = &node->domchg->domchgbound;
   if( domchgbound->nboundchgs == 0 )
      return;

   (*var) = domchgbound->boundchgs[0].var;
   (*bound) = domchgbound->boundchgs[0].newbound;
   (*boundtype) = (SCIP_BOUNDTYPE) domchgbound->boundchgs[0].boundtype;
}

/** creates visualization data structure */
SCIP_RETCODE SCIPvisualCreate(
   SCIP_VISUAL**         visual,             /**< pointer to store visualization information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_ALLOC( BMSallocMemory(visual) );

   (*visual)->vbcfile = NULL;
   (*visual)->bakfile = NULL;
   (*visual)->messagehdlr = messagehdlr;
   (*visual)->nodenum = NULL;
   (*visual)->timestep = 0;
   (*visual)->lastnode = NULL;
   (*visual)->lastcolor = SCIP_VBCCOLOR_NONE;
   (*visual)->userealtime = FALSE;

   return SCIP_OKAY;
}

/** frees visualization data structure */
void SCIPvisualFree(
   SCIP_VISUAL**         visual              /**< pointer to store visualization information */
   )
{
   assert( visual != NULL );
   assert( *visual != NULL );
   assert( (*visual)->vbcfile == NULL );
   assert( (*visual)->bakfile == NULL );
   assert( (*visual)->nodenum == NULL );

   BMSfreeMemory(visual);
}

/** initializes visualization information and creates a file for visualization output */
SCIP_RETCODE SCIPvisualInit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert( visual != NULL );
   assert( set != NULL );
   assert( set->visual_vbcfilename != NULL );
   assert( set->visual_bakfilename != NULL );
   assert( visual->nodenum == NULL );

   /* check whether we should initialize VBC output */
   if ( set->visual_vbcfilename[0] != '-' || set->visual_vbcfilename[1] != '\0' )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
         "storing VBC information in file <%s>\n", set->visual_vbcfilename);
      visual->vbcfile = fopen(set->visual_vbcfilename, "w");
      visual->timestep = 0;
      visual->lastnode = NULL;
      visual->lastcolor = SCIP_VBCCOLOR_NONE;
      visual->userealtime = set->visual_realtime;

      if( visual->vbcfile == NULL )
      {
         SCIPerrorMessage("error creating file <%s>\n", set->visual_vbcfilename);
         SCIPprintSysError(set->visual_vbcfilename);
         return SCIP_FILECREATEERROR;
      }

      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#TYPE: COMPLETE TREE\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#TIME: SET\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#BOUNDS: SET\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#INFORMATION: STANDARD\n");
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "#NODE_NUMBER: NONE\n");
   }

   /* check whether we should initialize BAK output */
   if ( set->visual_bakfilename[0] != '-' || set->visual_bakfilename[1] != '\0' )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
         "storing BAK information in file <%s>\n", set->visual_bakfilename);
      visual->bakfile = fopen(set->visual_bakfilename, "w");
      visual->timestep = 0;
      visual->lastnode = NULL;
      visual->lastcolor = SCIP_VBCCOLOR_NONE;
      visual->userealtime = set->visual_realtime;

      if ( visual->bakfile == NULL )
      {
         SCIPerrorMessage("error creating file <%s>\n", set->visual_bakfilename);
         SCIPprintSysError(set->visual_bakfilename);
         return SCIP_FILECREATEERROR;
      }
   }

   /* possibly init hashmap for nodes */
   if ( visual->vbcfile != NULL || visual->bakfile != NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&visual->nodenum, blkmem, SCIP_HASHSIZE_VBC) );
   }

   return SCIP_OKAY;
}

/** closes the visualization output file */
void SCIPvisualExit(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
  )
{
   assert( visual != NULL );
   assert( set != NULL );

   if ( visual->vbcfile != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL, "closing VBC information file\n");

      fclose(visual->vbcfile);
      visual->vbcfile = NULL;
   }

   if ( visual->bakfile != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL, "closing BAK information file\n");

      fclose(visual->bakfile);
      visual->bakfile = NULL;
   }

   if ( visual->nodenum )
      SCIPhashmapFree(&visual->nodenum);
}

/** prints current solution time to visualization output file */
static
void printTime(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Longint step;
   int hours;
   int mins;
   int secs;
   int hunds;

   assert( visual != NULL );
   assert( stat != NULL );

   if( visual->userealtime )
   {
      double time;
      time = SCIPclockGetTime(stat->solvingtime);
      step = (SCIP_Longint)(time * 100.0);
   }
   else
   {
      step = visual->timestep;
      visual->timestep++;
   }

   if ( visual->vbcfile != NULL )
   {
      hours = (int)(step / (60*60*100));
      step %= 60*60*100;
      mins = (int)(step / (60*100));
      step %= 60*100;
      secs = (int)(step / 100);
      step %= 100;
      hunds = (int)step;

      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "%02d:%02d:%02d.%02d ", hours, mins, secs, hunds);
   }

   if ( visual->bakfile != NULL )
   {
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "%f ", (SCIP_Real) step/100.0);
   }
}

/** creates a new node entry in the visualization output file */
SCIP_RETCODE SCIPvisualNewChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   size_t parentnodenum;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return SCIP_OKAY;

   /* insert mapping node -> nodenum into hash map */
   if( stat->ncreatednodesrun >= (SCIP_Longint)INT_MAX )
   {
      SCIPerrorMessage("too many nodes to store in the visualization file\n");
      return SCIP_INVALIDDATA;
   }

   nodenum = (size_t)stat->ncreatednodesrun;
   assert(nodenum > 0);
   SCIP_CALL( SCIPhashmapInsert(visual->nodenum, node, (void*)nodenum) );

   /* get nodenum of parent node from hash map */
   parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, node->parent) : 0);
   assert(node->parent == NULL || parentnodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat);
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "N %d %d %d\n", (int)parentnodenum, (int)nodenum, SCIP_VBCCOLOR_UNSOLVED);
      printTime(visual, stat);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar), SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, SCIPnodeGetLowerbound(node));
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node));
      }
   }

   if ( visual->bakfile != NULL )
   {
      char t = 'M';

      if ( branchvar != NULL )
         t = branchtype == SCIP_BOUNDTYPE_LOWER ? 'R' : 'L';

      /* todo: get information about fractionalities */
      if ( branchvar != NULL || parentnodenum == 0 )
      {
         printTime(visual, stat);
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "branched %d %d %c %f %f %d\n", (int)nodenum, (int)parentnodenum, t,
            SCIPnodeGetLowerbound(node), 0.0, 0);
      }
   }

   return SCIP_OKAY;
}

/** updates a node entry in the visualization output file */
SCIP_RETCODE SCIPvisualUpdateChild(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return SCIP_OKAY;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar), SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, SCIPnodeGetLowerbound(node));
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node));
      }
   }

   if ( visual->bakfile != NULL )
   {
      size_t parentnodenum;
      char t;

      /* determine branching type */
      if ( branchvar != NULL )
         t = branchtype == SCIP_BOUNDTYPE_LOWER ? 'R' : 'L';
      else
         t = 'M';

      /* get nodenum of parent node from hash map */
      parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, node->parent) : 0);
      assert(node->parent == NULL || parentnodenum > 0);

      /* todo: get information about fractionalities */
      printTime(visual, stat);
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "branched %d %d %c %f %f %d\n", (int)nodenum, (int)parentnodenum, t,
         SCIPnodeGetLowerbound(node), 0.0, 0);
   }

   return SCIP_OKAY;
}

/** changes the color of the node to the given color */
static
void vbcSetColor(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node to change color for */
   SCIP_VBCCOLOR         color               /**< new color of node, or SCIP_VBCCOLOR_NONE */
   )
{
   assert( visual != NULL );
   assert( node != NULL );

   if( visual->vbcfile != NULL && color != SCIP_VBCCOLOR_NONE && (node != visual->lastnode || color != visual->lastcolor) )
   {
      size_t nodenum;

      nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
      assert(nodenum > 0);
      printTime(visual, stat);
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "P %d %d\n", (int)nodenum, color);
      visual->lastnode = node;
      visual->lastcolor = color;
   }
}

/** marks node as solved in visualization output file */
void SCIPvisualSolvedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was solved */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\\nnr:\\t%"SCIP_LONGINT_FORMAT"\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar),  SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, SCIPnodeGetLowerbound(node), stat->nnodes);
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\\nnr:\\t%"SCIP_LONGINT_FORMAT"\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node), stat->nnodes);
      }
      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_SOLVED);
   }

   /* do nothing for BAK */
}

/** changes the color of the node to the color of cutoff nodes */
void SCIPvisualCutoffNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was cut off */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   size_t nodenum;

   assert( visual != NULL );
   assert( stat != NULL );
   assert( node != NULL );

   /* check whether output should be created */
   if ( visual->vbcfile == NULL && visual->bakfile == NULL )
      return;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   if ( visual->vbcfile != NULL )
   {
      printTime(visual, stat);
      if( branchvar != NULL )
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\\nnr:\\t%"SCIP_LONGINT_FORMAT"\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
            SCIPvarGetName(branchvar),  SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
            branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, SCIPnodeGetLowerbound(node), stat->nnodes);
      }
      else
      {
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\\nnr:\\t%"SCIP_LONGINT_FORMAT"\n",
            (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node), stat->nnodes);
      }
      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_CUTOFF);
   }

   if ( visual->bakfile != NULL )
   {
      size_t parentnodenum;
      char t;

      /* determine branching type */
      if ( branchvar != NULL )
         t = branchtype == SCIP_BOUNDTYPE_LOWER ? 'R' : 'L';
      else
         t = 'M';

      /* get nodenum of parent node from hash map */
      parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, node->parent) : 0);
      assert(node->parent == NULL || parentnodenum > 0);

      /* todo: distinguish between infeasible and fathomed nodes */
      printTime(visual, stat);
      SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "infeasible %d %d %c\n", (int)nodenum, (int)parentnodenum, t);
   }
}

/** changes the color of the node to the color of nodes where a conflict constraint was found */
void SCIPvisualFoundConflict(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, where the conflict was found */
   )
{
   assert(node != NULL);

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_CONFLICT);

   /* do nothing for BAK */
}

/** changes the color of the node to the color of nodes that were marked to be repropagated */
void SCIPvisualMarkedRepropagateNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was marked to be repropagated */
   )
{
   assert(node != NULL);

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* if the node number is zero, then SCIP is currently in probing and wants to mark a probing node; however this node
    * is not part of the search tree */
   if( SCIPnodeGetNumber(node) > 0 )
      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_MARKREPROP);

   /* do nothing for BAK */
}

/** changes the color of the node to the color of repropagated nodes */
void SCIPvisualRepropagatedNode(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was repropagated */
   )
{
   assert(node != NULL);

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_REPROP);

   /* do nothing for BAK */
}

/** changes the color of the node to the color of nodes with a primal solution */
void SCIPvisualFoundSolution(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node where the solution was found, or NULL */
   )
{
   if ( visual->vbcfile != NULL )
   {
      if ( node == NULL || ! set->visual_dispsols )
         return;

      if ( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
         return;

      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_SOLUTION);
   }

   if ( visual->bakfile != NULL )
   {
      SCIP_SOL* bestsol = NULL;

      bestsol = SCIPgetBestSol(set->scip);
      assert( bestsol != NULL );

      /* todo: also handle case that we find primal solution within probing - this seems to happen often. */

      /* if we cannot determine node information or a heuristic found the solution */
      if ( node == NULL || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE || SCIPsolGetHeur(bestsol) != NULL )
      {
         printTime(visual, stat);
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "heuristic %f\n", SCIPgetSolTransObj(set->scip, bestsol));
      }
      else
      {
         /* if LP solution was feasible ... */
         SCIP_VAR* branchvar;
         SCIP_BOUNDTYPE branchtype;
         SCIP_Real branchbound;
         size_t parentnodenum;
         size_t nodenum;
         char t = 'M';

         /* get node num from hash map */
         nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);

         /* get nodenum of parent node from hash map */
         parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(visual->nodenum, node->parent) : 0);
         assert(node->parent == NULL || parentnodenum > 0);

         /* get branching information */
         getBranchInfo(node, &branchvar, &branchtype, &branchbound);

         /* determine branching type */
         if ( branchvar != NULL )
            t = branchtype == SCIP_BOUNDTYPE_LOWER ? 'R' : 'L';

         printTime(visual, stat);
         SCIPmessageFPrintInfo(visual->messagehdlr, visual->bakfile, "integer %d %d %c\n", (int)nodenum, (int)parentnodenum, t, SCIPgetUpperbound(set->scip));
      }
   }
}

/** outputs a new global lower bound to the visualization output file */
void SCIPvisualLowerbound(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             lowerbound          /**< new lower bound */
   )
{
   assert( visual != NULL );

   /* check, if VBC output should be created */
   if( visual->vbcfile == NULL )
      return;

   printTime(visual, stat);
   SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "L %f\n", lowerbound);

   /* do nothing for BAK */
}

/** outputs a new global upper bound to the visualization output file */
void SCIPvisualUpperbound(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             upperbound          /**< new upper bound */
   )
{
   assert( visual != NULL );

   /* check, if VBC output should be created */
   if( visual->vbcfile == NULL )
      return;

   printTime(visual, stat);
   SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "U %f\n", upperbound);

   /* do nothing for BAK */
}
