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
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: vbc.c,v 1.8 2005/01/18 09:26:59 bzfpfend Exp $"

/**@file   vbc.c
 * @brief  methods for VBC Tool output
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>

#include "memory.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
#include "misc.h"
#include "var.h"
#include "tree.h"
#include "vbc.h"
#include "struct_vbc.h"




/** node colors in VBC output:
 *   1: indian red
 *   2: green
 *   3: light gray
 *   4: red
 *   5: blue
 *   6: black
 *   7: light pink
 *   8: cyan
 *   9: dark green
 *  10: brown
 *  11: orange
 *  12: yellow
 *  13: pink
 *  14: purple
 *  15: light blue
 *  16: muddy green
 *  17: white
 *  18: light grey
 *  19: light grey
 *  20: light grey
 */
enum VBCColor
{
   SCIP_VBCCOLOR_UNSOLVED   =  3,       /**< color for newly created, unsolved nodes */
   SCIP_VBCCOLOR_SOLVED     =  2,       /**< color for solved nodes */
   SCIP_VBCCOLOR_CUTOFF     =  4,       /**< color for nodes that were cut off */
   SCIP_VBCCOLOR_CONFLICT   = 15,       /**< color for nodes where a conflict clause was found */
   SCIP_VBCCOLOR_MARKREPROP = 11,       /**< color for nodes that were marked to be repropagated */
   SCIP_VBCCOLOR_REPROP     = 12,       /**< color for repropagated nodes */
   SCIP_VBCCOLOR_SOLUTION   = -1        /**< color for solved nodes, where a solution has been found */
};
typedef enum VBCColor VBCCOLOR;


/** returns the branching variable of the node, or NULL */
static
VAR* getBranchVar(
   NODE*            node                /**< new node, that was created */
   )
{
   DOMCHGBOUND* domchgbound;

   assert(node != NULL);
   if( node->domchg == NULL )
      return NULL;
   
   domchgbound = &node->domchg->domchgbound;
   if( domchgbound->nboundchgs == 0 )
      return NULL;

   return domchgbound->boundchgs[0].var;
}

/** creates VBC Tool data structure */
RETCODE SCIPvbcCreate(
   VBC**            vbc                 /**< pointer to store the VBC information */
   )
{
   ALLOC_OKAY( allocMemory(vbc) );

   (*vbc)->file = NULL;
   (*vbc)->nodenum = NULL;
   (*vbc)->timestep = 0;
   (*vbc)->userealtime = FALSE;

   return SCIP_OKAY;
}

/** frees VBC Tool data structure */
void SCIPvbcFree(
   VBC**            vbc                 /**< pointer to store the VBC information */
   )
{
   assert(vbc != NULL);
   assert(*vbc != NULL);
   assert((*vbc)->file == NULL);
   assert((*vbc)->nodenum == NULL);

   freeMemory(vbc);
}

/** initializes VBC information and creates a file for VBC output */
RETCODE SCIPvbcInit(
   VBC*             vbc,                /**< VBC information */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(vbc != NULL);
   assert(set != NULL);
   assert(set->vbc_filename != NULL);

   if( set->vbc_filename[0] == '-' && set->vbc_filename[1] == '\0' )
      return SCIP_OKAY;

   infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "storing VBC information in file <%s>\n", set->vbc_filename);
   vbc->file = fopen(set->vbc_filename, "w");
   vbc->timestep = 0;
   vbc->userealtime = set->vbc_realtime;

   if( vbc->file == NULL )
   {
      errorMessage("error creating file <%s>\n", set->vbc_filename);
      return SCIP_FILECREATEERROR;
   }

   fprintf(vbc->file, "#TYPE: COMPLETE TREE\n");
   fprintf(vbc->file, "#TIME: SET\n");
   fprintf(vbc->file, "#BOUNDS: SET\n");
   fprintf(vbc->file, "#INFORMATION: STANDARD\n");
   fprintf(vbc->file, "#NODE_NUMBER: NONE\n");

   CHECK_OKAY( SCIPhashmapCreate(&vbc->nodenum, memhdr, SCIP_HASHSIZE_VBC) );

   return SCIP_OKAY;
}

/** closes the VBC output file */
void SCIPvbcExit(
   VBC*             vbc,                /**< VBC information */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(vbc != NULL);
   assert(set != NULL);

   if( vbc->file != NULL )
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL, "closing VBC information file\n");
      
      fclose(vbc->file);
      vbc->file = NULL;

      SCIPhashmapFree(&vbc->nodenum);
   }
}

/** prints current solution time to VBC output file */
static
void printTime(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat                /**< problem statistics */
   )
{
   Longint step;
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
      step = (Longint)(time * 100.0);
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

   fprintf(vbc->file, "%02d:%02d:%02d.%02d ", hours, mins, secs, hunds);
}

/** creates a new node entry in the VBC output file */
RETCODE SCIPvbcNewChild(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   VAR* branchvar;
   int parentnodenum;
   int nodenum;

   assert(vbc != NULL);
   assert(stat != NULL);
   assert(node != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return SCIP_OKAY;

   /* insert mapping node -> nodenum into hash map */
   nodenum = stat->ncreatednodesrun;
   assert(nodenum > 0);
   CHECK_OKAY( SCIPhashmapInsert(vbc->nodenum, node, (void*)nodenum) );

   /* get nodenum of parent node from hash map */
   parentnodenum = node->parent != NULL ? (int)SCIPhashmapGetImage(vbc->nodenum, node->parent) : 0;
   assert(node->parent == NULL || parentnodenum > 0);

   /* get branching variable */
   branchvar = getBranchVar(node);

   printTime(vbc, stat);
   fprintf(vbc->file, "N %d %d %d\n", parentnodenum, nodenum, SCIP_VBCCOLOR_UNSOLVED);
   printTime(vbc, stat);
   fprintf(vbc->file, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s\\nbound:\\t%f\n",
      nodenum, nodenum, node, SCIPnodeGetDepth(node),
      branchvar == NULL ? "-" : SCIPvarGetName(branchvar), SCIPnodeGetLowerbound(node));

   return SCIP_OKAY;
}

/** changes the color of the node to the given color */
static
void vbcSetColor(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node,               /**< node to change color for */
   VBCCOLOR         color               /**< new color of node, or -1 */
   )
{
   assert(vbc != NULL);
   assert(node != NULL);

   if( vbc->file != NULL && color != -1 )
   {
      int nodenum;

      nodenum = (int)SCIPhashmapGetImage(vbc->nodenum, node);
      assert(nodenum > 0);
      printTime(vbc, stat);
      fprintf(vbc->file, "P %d %d\n", nodenum, color);
   }
}

/** changes the color of the node to the color of solved nodes */
void SCIPvbcSolvedNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   VAR* branchvar;
   int nodenum;

   assert(vbc != NULL);
   assert(stat != NULL);
   assert(node != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return;

   /* get node num from hash map */
   nodenum = (int)SCIPhashmapGetImage(vbc->nodenum, node);
   assert(nodenum > 0);

   /* get branching variable */
   branchvar = getBranchVar(node);

   printTime(vbc, stat);
   fprintf(vbc->file, "I %d \\inode:\\t%d (%p)\\idepth:\\t%d\\nvar:\\t%s\\nbound:\\t%f\\nnr:\\t%lld\n", 
      nodenum, nodenum, node, SCIPnodeGetDepth(node),
      branchvar == NULL ? "-" : SCIPvarGetName(branchvar), SCIPnodeGetLowerbound(node), stat->nnodes);

   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_SOLVED);
}

/** changes the color of the node to the color of cutoff nodes */
void SCIPvbcCutoffNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_CUTOFF);
}

/** changes the color of the node to the color of nodes where a conflict clause was found */
void SCIPvbcFoundConflict(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_CONFLICT);
}

/** changes the color of the node to the color of nodes that were marked to be repropagated */
void SCIPvbcMarkedRepropagateNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_MARKREPROP);
}

/** changes the color of the node to the color of repropagated nodes */
void SCIPvbcRepropagatedNode(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_REPROP);
}

/** changes the color of the node to the color of nodes with a primal solution */
void SCIPvbcFoundSolution(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   vbcSetColor(vbc, stat, node, SCIP_VBCCOLOR_SOLUTION);
}

/** outputs a new global lower bound to the VBC output file */
void SCIPvbcLowerbound(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   Real             lowerbound          /**< new lower bound */
   )
{
   assert(vbc != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return;

   printTime(vbc, stat);
   fprintf(vbc->file, "L %f\n", lowerbound);
}

/** outputs a new global upper bound to the VBC output file */
void SCIPvbcUpperbound(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   Real             upperbound          /**< new upper bound */
   )
{
   assert(vbc != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return;

   printTime(vbc, stat);
   fprintf(vbc->file, "U %f\n", upperbound);
}

