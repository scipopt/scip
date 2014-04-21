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
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "blockmemshell/memory.h"
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

   if( set->visual_vbcfilename[0] == '-' && set->visual_vbcfilename[1] == '\0' )
      return SCIP_OKAY;

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

   SCIP_CALL( SCIPhashmapCreate(&visual->nodenum, blkmem, SCIP_HASHSIZE_VBC) );

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

   if( visual->vbcfile != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL, "closing VBC information file\n");

      fclose(visual->vbcfile);
      visual->vbcfile = NULL;

      SCIPhashmapFree(&visual->nodenum);
   }
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
   hours = (int)(step / (60*60*100));
   step %= 60*60*100;
   mins = (int)(step / (60*100));
   step %= 60*100;
   secs = (int)(step / 100);
   step %= 100;
   hunds = (int)step;

   SCIPmessageFPrintInfo(visual->messagehdlr, visual->vbcfile, "%02d:%02d:%02d.%02d ", hours, mins, secs, hunds);
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

   /* check whether VBC output should be created */
   if( visual->vbcfile == NULL )
      return SCIP_OKAY;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
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

   /* check, if VBC output should be created */
   if( visual->vbcfile == NULL )
      return SCIP_OKAY;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

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
   assert(node != NULL);

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

/** changes the color of the node to the color of solved nodes */
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
   assert(stat != NULL);
   assert(node != NULL);

   /* check, if VBC output should be created */
   if( visual->vbcfile == NULL )
      return;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

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
   assert(stat != NULL);
   assert(node != NULL);

   /* check, if VBC output should be created */
   if( visual->vbcfile == NULL )
      return;

   /* visualization is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(visual->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

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
}

/** changes the color of the node to the color of nodes with a primal solution */
void SCIPvisualFoundSolution(
   SCIP_VISUAL*          visual,             /**< visualization information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node where the solution was found, or NULL */
   )
{
   if( node != NULL && set->visual_dispsols )
   {
      /* visualization is disabled on probing nodes */
      if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
	 return;

      vbcSetColor(visual, stat, node, SCIP_VBCCOLOR_SOLUTION);
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
}
