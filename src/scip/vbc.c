/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: vbc.c,v 1.4 2004/05/03 11:26:57 bzfpfend Exp $"

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
   assert(set->vbcfilename != NULL);

   if( set->vbcfilename[0] == '-' && set->vbcfilename[1] == '\0' )
      return SCIP_OKAY;

   infoMessage(set->verblevel, SCIP_VERBLEVEL_NORMAL, "storing VBC information in file <%s>\n", set->vbcfilename);
   vbc->file = fopen(set->vbcfilename, "w");

   if( vbc->file == NULL )
   {
      errorMessage("error creating file <%s>\n", set->vbcfilename);
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
      infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL, "closing VBC information file\n");
      
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
   double time;
   int hours;
   int mins;
   int secs;
   int hunds;
   
   assert(vbc != NULL);
   assert(stat != NULL);

   time = SCIPclockGetTime(stat->solvingtime);
   hunds = (int)(time * 100.0);
   hours = hunds / (60*60*100);
   hunds %= 60*60*100;
   mins = hunds / (60*100);
   hunds %= 60*100;
   secs = hunds / 100;
   hunds %= 100;

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
   nodenum = stat->ncreatednodes;
   CHECK_OKAY( SCIPhashmapInsert(vbc->nodenum, node, (void*)nodenum) );

   /* get nodenum of parent node from hash map */
   parentnodenum = node->parent != NULL ? (int)SCIPhashmapGetImage(vbc->nodenum, node->parent) : 0;

   /* get branching variable */
   branchvar = getBranchVar(node);

   printTime(vbc, stat);
   fprintf(vbc->file, "N %d %d %d\n", parentnodenum, nodenum, SCIP_VBCCOLOR_UNSOLVED);
   printTime(vbc, stat);
   fprintf(vbc->file, "I %d \\inode:\\t%d\\ivar:\\t%s\\nbound:\\t%f\n",
      nodenum, nodenum, branchvar == NULL ? "-" : SCIPvarGetName(branchvar), SCIPnodeGetLowerbound(node));

   return SCIP_OKAY;
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

   /* get branching variable */
   branchvar = getBranchVar(node);

   printTime(vbc, stat);
   fprintf(vbc->file, "I %d \\inode:\\t%d\\ivar:\\t%s\\nbound:\\t%f\\nnr:\\t%lld\n", 
      nodenum, nodenum, branchvar == NULL ? "-" : SCIPvarGetName(branchvar), SCIPnodeGetLowerbound(node), stat->nnodes);

   printTime(vbc, stat);
   fprintf(vbc->file, "P %d %d\n", nodenum, SCIP_VBCCOLOR_SOLVED);
}

/** changes the color of the node to the color of nodes with a primal solution */
void SCIPvbcFoundSolution(
   VBC*             vbc,                /**< VBC information */
   STAT*            stat,               /**< problem statistics */
   NODE*            node                /**< new node, that was created */
   )
{
   int nodenum;

   assert(vbc != NULL);
   assert(stat != NULL);
   assert(node != NULL);

   /* check, if VBC output should be created */
   if( vbc->file == NULL )
      return;

   /* get node num from hash map */
   nodenum = (int)SCIPhashmapGetImage(vbc->nodenum, node);

   printTime(vbc, stat);
   fprintf(vbc->file, "P %d %d\n", nodenum, SCIP_VBCCOLOR_SOLUTION);
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

