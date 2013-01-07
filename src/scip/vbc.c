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

/**@file   vbc.c
 * @brief  methods for VBC Tool output
 * @author Tobias Achterberg
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
#include "scip/vbc.h"
#include "scip/struct_vbc.h"


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

/** creates VBC Tool data structure */
SCIP_RETCODE SCIPvbcCreate(
   SCIP_VBC**            vbc,                /**< pointer to store the VBC information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_ALLOC( BMSallocMemory(vbc) );

   (*vbc)->file = NULL;
   (*vbc)->messagehdlr = messagehdlr;
   (*vbc)->nodenum = NULL;
   (*vbc)->timestep = 0;
   (*vbc)->lastnode = NULL;
   (*vbc)->lastcolor = SCIP_VBCCOLOR_NONE;
   (*vbc)->userealtime = FALSE;

   return SCIP_OKAY;
}

/** frees VBC Tool data structure */
void SCIPvbcFree(
   SCIP_VBC**            vbc                 /**< pointer to store the VBC information */
   )
{
   assert(vbc != NULL);
   assert(*vbc != NULL);
   assert((*vbc)->file == NULL);
   assert((*vbc)->nodenum == NULL);

   BMSfreeMemory(vbc);
}

/** initializes VBC information and creates a file for VBC output */
SCIP_RETCODE SCIPvbcInit(
   SCIP_VBC*             vbc,                /**< VBC information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert(vbc != NULL);
   assert(set != NULL);
   assert(set->vbc_filename != NULL);

   if( set->vbc_filename[0] == '-' && set->vbc_filename[1] == '\0' )
      return SCIP_OKAY;

   SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "storing VBC information in file <%s>\n", set->vbc_filename);
   vbc->file = fopen(set->vbc_filename, "w");
   vbc->timestep = 0;
   vbc->lastnode = NULL;
   vbc->lastcolor = SCIP_VBCCOLOR_NONE;
   vbc->userealtime = set->vbc_realtime;

   if( vbc->file == NULL )
   {
      SCIPerrorMessage("error creating file <%s>\n", set->vbc_filename);
      SCIPprintSysError(set->vbc_filename);
      return SCIP_FILECREATEERROR;
   }

   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "#TYPE: COMPLETE TREE\n");
   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "#TIME: SET\n");
   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "#BOUNDS: SET\n");
   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "#INFORMATION: STANDARD\n");
   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "#NODE_NUMBER: NONE\n");

   SCIP_CALL( SCIPhashmapCreate(&vbc->nodenum, blkmem, SCIP_HASHSIZE_VBC) );

   return SCIP_OKAY;
}

/** closes the VBC output file */
void SCIPvbcExit(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
  )
{
   assert(vbc != NULL);
   assert(set != NULL);

   if( vbc->file != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_FULL, "closing VBC information file\n");

      fclose(vbc->file);
      vbc->file = NULL;

      SCIPhashmapFree(&vbc->nodenum);
   }
}

/** prints current solution time to VBC output file */
static
void printTime(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   SCIP_Longint step;
   int hours;
   int mins;
   int secs;
   int hunds;

   assert(vbc != NULL);
   assert(stat != NULL);

   if( vbc->userealtime )
   {
      double time;
      time = SCIPclockGetTime(stat->solvingtime);
      step = (SCIP_Longint)(time * 100.0);
   }
   else
   {
      step = vbc->timestep;
      vbc->timestep++;
   }
   hours = (int)(step / (60*60*100));
   step %= 60*60*100;
   mins = (int)(step / (60*100));
   step %= 60*100;
   secs = (int)(step / 100);
   step %= 100;
   hunds = (int)step;

   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "%02d:%02d:%02d.%02d ", hours, mins, secs, hunds);
}

/** creates a new node entry in the VBC output file */
SCIP_RETCODE SCIPvbcNewChild(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   size_t parentnodenum;
   size_t nodenum;

   assert(vbc != NULL);
   assert(stat != NULL);
   assert(node != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return SCIP_OKAY;

   /* vbc is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* insert mapping node -> nodenum into hash map */
   if( stat->ncreatednodesrun >= (SCIP_Longint)INT_MAX )
   {
      SCIPerrorMessage("too many nodes to store in the VBC file\n");
      return SCIP_INVALIDDATA;
   }

   nodenum = (size_t)stat->ncreatednodesrun;
   assert(nodenum > 0);
   SCIP_CALL( SCIPhashmapInsert(vbc->nodenum, node, (void*)nodenum) );

   /* get nodenum of parent node from hash map */
   parentnodenum = (node->parent != NULL ? (size_t)SCIPhashmapGetImage(vbc->nodenum, node->parent) : 0);
   assert(node->parent == NULL || parentnodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   printTime(vbc, stat);
   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "N %d %d %d\n", (int)parentnodenum, (int)nodenum, SCIP_VBCCOLOR_UNSOLVED);
   printTime(vbc, stat);
   if( branchvar != NULL )
   {
      SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\n",
         (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
         SCIPvarGetName(branchvar), SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
         branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, SCIPnodeGetLowerbound(node));
   }
   else
   {
      SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\n",
         (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node));
   }

   return SCIP_OKAY;
}

/** changes the color of the node to the given color */
static
void vbcSetColor(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node,               /**< node to change color for */
   SCIP_VBCCOLOR         color               /**< new color of node, or SCIP_VBCCOLOR_NONE */
   )
{
   assert(vbc != NULL);
   assert(node != NULL);

   if( vbc->file != NULL && color != SCIP_VBCCOLOR_NONE && (node != vbc->lastnode || color != vbc->lastcolor) )
   {
      size_t nodenum;

      nodenum = (size_t)SCIPhashmapGetImage(vbc->nodenum, node);
      assert(nodenum > 0);
      printTime(vbc, stat);
      SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "P %d %d\n", (int)nodenum, color);
      vbc->lastnode = node;
      vbc->lastcolor = color;
   }
}

/** changes the color of the node to the color of solved nodes */
void SCIPvbcSolvedNode(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was solved */
   )
{
   SCIP_VAR* branchvar;
   SCIP_BOUNDTYPE branchtype;
   SCIP_Real branchbound;
   size_t nodenum;

   assert(vbc != NULL);
   assert(stat != NULL);
   assert(node != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return;

   /* vbc is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* get node num from hash map */
   nodenum = (size_t)SCIPhashmapGetImage(vbc->nodenum, node);
   assert(nodenum > 0);

   /* get branching information */
   getBranchInfo(node, &branchvar, &branchtype, &branchbound);

   printTime(vbc, stat);

   if( branchvar != NULL )
   {
      SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s [%g,%g] %s %f\\nbound:\\t%f\\nnr:\\t%"SCIP_LONGINT_FORMAT"\n",
         (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node),
         SCIPvarGetName(branchvar),  SCIPvarGetLbLocal(branchvar), SCIPvarGetUbLocal(branchvar),
         branchtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",  branchbound, SCIPnodeGetLowerbound(node), stat->nnodes);
   }
   else
   {
      SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t-\\nbound:\\t%f\\nnr:\\t%"SCIP_LONGINT_FORMAT"\n",
         (int)nodenum, (int)nodenum, node, SCIPnodeGetDepth(node), SCIPnodeGetLowerbound(node), stat->nnodes);
   }

   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_SOLVED);
}

/** changes the color of the node to the color of cutoff nodes */
void SCIPvbcCutoffNode(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was cut off */
   )
{
   assert(node != NULL);

   /* vbc is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_CUTOFF);
}

/** changes the color of the node to the color of nodes where a conflict constraint was found */
void SCIPvbcFoundConflict(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, where the conflict was found */
   )
{
   assert(node != NULL);

   /* vbc is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_CONFLICT);
}

/** changes the color of the node to the color of nodes that were marked to be repropagated */
void SCIPvbcMarkedRepropagateNode(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was marked to be repropagated */
   )
{
   assert(node != NULL);

   /* vbc is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   /* if the node number is zero, then SCIP is currently in probing and wants to mark a probing node; however this node
    * is not part of the search tree */
   if( SCIPnodeGetNumber(node) > 0 )
      vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_MARKREPROP);
}

/** changes the color of the node to the color of repropagated nodes */
void SCIPvbcRepropagatedNode(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node, that was repropagated */
   )
{
   assert(node != NULL);

   /* vbc is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return;

   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_REPROP);
}

/** changes the color of the node to the color of nodes with a primal solution */
void SCIPvbcFoundSolution(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< node where the solution was found, or NULL */
   )
{
   if( node != NULL && set->vbc_dispsols )
   {
      /* vbc is disabled on probing nodes */
      if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
	 return;

      vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_SOLUTION);
   }
}

/** outputs a new global lower bound to the VBC output file */
void SCIPvbcLowerbound(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             lowerbound          /**< new lower bound */
   )
{
   assert(vbc != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return;

   printTime(vbc, stat);
   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "L %f\n", lowerbound);
}

/** outputs a new global upper bound to the VBC output file */
void SCIPvbcUpperbound(
   SCIP_VBC*             vbc,                /**< VBC information */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_Real             upperbound          /**< new upper bound */
   )
{
   assert(vbc != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return;

   printTime(vbc, stat);
   SCIPmessageFPrintInfo(vbc->messagehdlr, vbc->file, "U %f\n", upperbound);
}

