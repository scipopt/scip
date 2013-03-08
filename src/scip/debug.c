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

/**@file   debug.c
 * @brief  methods for debugging
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/debug.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/struct_scip.h"

#ifdef SCIP_DEBUG_SOLUTION

#define SCIP_HASHSIZE_DEBUG        131101    /**< minimum size of hash map for storing whether a solution is valid for the node */

static char** solnames = NULL;
static SCIP_Real* solvals = NULL;
static int nsolvals = 0;
static int solsize = 0;
static SCIP_SET* mainscipset = NULL;
static SCIP_HASHMAP* solinnode = NULL;       /**< maps nodes to bools, storing whether the solution is valid for the node */
static SCIP_Bool falseptr = FALSE;
static SCIP_Bool trueptr = TRUE;
static SCIP_Bool solisachieved = FALSE;      /**< means if current best solution is better than the given debug solution */
static SCIP_Real debugsolval = 0.0;          /**< objective value for debug solution */
static SCIP_Bool debugsoldisabled = FALSE;   /**< flag indicating if debugging of solution was disabled or not */

/** reads solution from given file into given arrays */
static
SCIP_RETCODE readSolfile(
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           solfilename,        /**< solution filename to read */
   char***               names,              /**< pointer to store the array of variable names */
   SCIP_Real**           vals,               /**< pointer to store the array of solution values */
   int*                  nvals,              /**< pointer to store the number of non-zero elements */
   int*                  valssize            /**< pointer to store the length of the variable names and solution values arrays */
   )
{
   FILE* file;
   int nonvalues;
   int i;

   assert(set != NULL);
   assert(solfilename != NULL);
   assert(names != NULL);
   assert(*names == NULL);
   assert(vals != NULL);
   assert(*vals == NULL);
   assert(nvals != NULL);
   assert(valssize != NULL);

   printf("***** debug: reading solution file <%s>\n", solfilename);

   /* open solution file */
   file = fopen(solfilename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open solution file <%s> specified in scip/debug.h\n", solfilename);
      SCIPprintSysError(solfilename);
      return SCIP_NOFILE;
   }

   /* read data */
   nonvalues = 0;
   *valssize = 0;

   while( !feof(file) )
   {
      char buf[SCIP_MAXSTRLEN];
      char name[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      SCIP_Real val;
      int nread;

      if( fgets(buf, SCIP_MAXSTRLEN, file) == NULL )
      {
         if( feof(file) )
            break;
         else
            return SCIP_READERROR;
      }

      /* the lines "solution status: ..." and "objective value: ..." may preceed the solution information */
      if( strncmp(buf, "solution", 8) == 0 || strncmp(buf, "objective", 9) == 0 )
      {
         nonvalues++;
         continue;
      }

      /* skip empty lines */
      if( strlen(buf) == 1 )
      {
         nonvalues++;
         continue;
      }


      nread = sscanf(buf, "%s %lf %s\n", name, &val, objstring);
      if( nread < 2 )
      {
         printf("invalid input line %d in solution file <%s>: <%s>\n", *nvals + nonvalues, SCIP_DEBUG_SOLUTION, name);
         fclose(file);
         return SCIP_READERROR;
      }

      /* allocate memory */
      if( *nvals >= *valssize )
      {
         *valssize = MAX(2 * *valssize, (*nvals)+1);
         SCIP_ALLOC( BMSreallocMemoryArray(names, *valssize) );
         SCIP_ALLOC( BMSreallocMemoryArray(vals, *valssize) );
      }
      assert(*nvals < *valssize);

      /* store solution value in sorted list */
      for( i = *nvals; i > 0 && strcmp(name, (*names)[i-1]) < 0; --i )
      {
         (*names)[i] = (*names)[i-1];
         (*vals)[i] = (*vals)[i-1];
      }
      SCIP_ALLOC( BMSduplicateMemoryArray(&(*names)[i], name, strlen(name)+1) );
      SCIPdebugMessage("found variable <%s>: value <%g>\n", (*names)[i], val);
      (*vals)[i] = val;
      (*nvals)++;
   }

   debugsolval = 0.0;

   /* get solution value */
   for( i = *nvals - 1; i >= 0; --i)
   {
      SCIP_VAR* var;
      var = SCIPfindVar(set->scip, (*names)[i]);
      if( var != NULL )
         debugsolval += (*vals)[i] * SCIPvarGetObj(var);
   }
   SCIPdebugMessage("Debug Solution value is %g.\n", debugsolval);

   /* close file */
   fclose(file);

   /* remember the set pointer to identify sub-MIP calls */
   mainscipset = set;

   printf("***** debug: read %d non-zero entries\n", *nvals);

   return SCIP_OKAY;
}

/** reads feasible solution to check from file */
static
SCIP_RETCODE readSolution(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(set != NULL);

   if( nsolvals > 0 )
      return SCIP_OKAY;

   SCIP_CALL( readSolfile(set, SCIP_DEBUG_SOLUTION, &solnames, &solvals, &nsolvals, &solsize) );

   return SCIP_OKAY;
}

/** gets value of given variable in debugging solution */
static
SCIP_RETCODE getSolutionValue(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to get solution value for */
   SCIP_Real*            val                 /**< pointer to store solution value */
   )
{
   SCIP_VAR* solvar;
   SCIP_Real scalar;
   SCIP_Real constant;
   const char* name;
   int left;
   int right;
   int middle;
   int cmp;

   assert(set != NULL);
   assert(var != NULL);
   assert(val != NULL);

   SCIP_CALL( readSolution(set) );
   SCIPdebugMessage("Now handling variable <%s>, which has status %d, is of type %d, and was deleted: %d, negated: %d, transformed: %d\n",
      SCIPvarGetName(var), SCIPvarGetStatus(var), SCIPvarGetType(var), SCIPvarIsDeleted(var), SCIPvarIsNegated(var),SCIPvarIsTransformedOrigvar(var));
   /* ignore deleted variables */
   if( SCIPvarIsDeleted(var) )
   {
      SCIPdebugMessage("**** unknown solution value for deleted variable <%s>\n", SCIPvarGetName(var));
      *val = SCIP_UNKNOWN;
      return SCIP_OKAY;
   }
   /* retransform variable onto original variable space */
   solvar = var;
   scalar = 1.0;
   constant = 0.0;
   if( SCIPvarIsNegated(solvar) )
   {
      scalar = -1.0;
      constant = SCIPvarGetNegationConstant(solvar);
      solvar = SCIPvarGetNegationVar(solvar);
   }
   if( SCIPvarIsTransformed(solvar) )
   {
      SCIP_CALL( SCIPvarGetOrigvarSum(&solvar, &scalar, &constant) );
      if( solvar == NULL )
      {
         /* if no original counterpart, then maybe someone added a value for the transformed variable, so search for var (or its negation) */
         SCIPdebugMessage("variable <%s> has no original counterpart\n", SCIPvarGetName(var));
         solvar = var;
         scalar = 1.0;
         constant = 0.0;
         if( SCIPvarIsNegated(solvar) )
         {
            scalar = -1.0;
            constant = SCIPvarGetNegationConstant(solvar);
            solvar = SCIPvarGetNegationVar(solvar);
         }
      }
   }
   /* perform a binary search for the variable */
   name = SCIPvarGetName(solvar);
   left = 0;
   right = nsolvals-1;
   while( left <= right )
   {
      middle = (left+right)/2;
      cmp = strcmp(name, solnames[middle]);
      if( cmp < 0 )
         right = middle-1;
      else if( cmp > 0 )
         left = middle+1;
      else
      {
         *val = scalar * solvals[middle] + constant;
         return SCIP_OKAY;
      }
   }
   *val = constant;

   if( *val < SCIPvarGetLbGlobal(var) - 1e-06 || *val > SCIPvarGetUbGlobal(var) + 1e-06 )
   {
      SCIPmessagePrintWarning(SCIPgetMessagehdlr(set->scip), "invalid solution value %.15g for variable <%s>[%.15g,%.15g]\n",
         *val, SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   return SCIP_OKAY;
}

/** gets value for a variable in the debug solution
 * if no value is stored for the variable, gives 0.0
 */
SCIP_RETCODE SCIPdebugGetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which to get the value */
   SCIP_Real*            val                 /**< buffer to store solution value */
   )
{
   SCIP_CALL( getSolutionValue(scip->set, var, val) );

   return SCIP_OKAY;
}

/** returns whether the debug solution is worse as the best known solution or if the debug solution was found */
static
SCIP_Bool debugSolIsAchieved(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_SOL* bestsol;
   SCIP* scip;

   if( solisachieved )
      return TRUE;

   assert(set != NULL);

   scip = set->scip;
   assert(scip != NULL);

   bestsol = SCIPgetBestSol(scip);

   if( bestsol != NULL )
   {
      SCIP_Real solvalue;

      /* don't check solution while in problem creation stage */
      if( SCIPsetGetStage(set) == SCIP_STAGE_PROBLEM )
         return TRUE;

      solvalue = SCIPgetSolOrigObj(scip, bestsol);

      /* make sure a debug solution has been read, so we do not compare against the initial debugsolval == 0 */
      SCIP_CALL( readSolution(set) );

      if( (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE && SCIPsetIsLE(set, solvalue, debugsolval)) || (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE && SCIPsetIsGE(set, solvalue, debugsolval)) )
         solisachieved = TRUE;
   }

   return solisachieved;
}

/** returns whether the solution belongs to the current SCIP instance */
static
SCIP_Bool isSolutionInMip(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   return (mainscipset == NULL || mainscipset == set);
}

/** returns whether the solution is contained in node's subproblem */
static
SCIP_RETCODE isSolutionInNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< local node where this bound change was applied */
   SCIP_Bool*            solcontained        /**< pointer to store whether the solution is contained in node's subproblem */
   )
{
   SCIP_Bool* boolptr;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(solcontained != NULL);

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
   {
      *solcontained = FALSE;
      return SCIP_OKAY;
   }

   /* generate the hashmap */
   if( solinnode == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&solinnode, blkmem, SCIPcalcHashtableSize(SCIP_HASHSIZE_DEBUG)) );
   }

   /* check, whether we know already whether the solution is contained in the given node */
   boolptr = (SCIP_Bool*)SCIPhashmapGetImage(solinnode, (void*)node);
   if( boolptr != NULL )
   {
      if( boolptr != &falseptr && boolptr != &trueptr )
      {
         SCIPerrorMessage("wrong value in node hashmap\n");
         SCIPABORT();
      }
      *solcontained = *boolptr;
      return SCIP_OKAY;
   }

   /* if the solution is not contained in the parent of the node, it cannot be contained in the current node */
   *solcontained = TRUE;
   if( node->parent != NULL )
   {
      SCIP_CALL( isSolutionInNode(blkmem, set, node->parent, solcontained) );
   }

   if( *solcontained )
   {
      /* check whether the bound changes at the current node remove the debugging solution from the subproblem */
      if( node->domchg != NULL )
      {
         SCIP_DOMCHGBOUND* domchgbound;
         SCIP_BOUNDCHG* boundchgs;
         int i;

         domchgbound = &node->domchg->domchgbound;
         boundchgs = domchgbound->boundchgs;
         for( i = 0; i < (int)domchgbound->nboundchgs && *solcontained; ++i )
         {
            SCIP_Real varsol;

            /* get solution value of variable */
            SCIP_CALL( getSolutionValue(set, boundchgs[i].var, &varsol) );

            if( varsol != SCIP_UNKNOWN ) /*lint !e777*/
            {
               /* compare the bound change with the solution value */
               if( SCIPboundchgGetBoundtype(&boundchgs[i]) == SCIP_BOUNDTYPE_LOWER )
                  *solcontained = SCIPsetIsFeasGE(set, varsol, boundchgs[i].newbound);
               else
                  *solcontained = SCIPsetIsFeasLE(set, varsol, boundchgs[i].newbound);

               if( !(*solcontained) && SCIPboundchgGetBoundchgtype(&boundchgs[i]) != SCIP_BOUNDCHGTYPE_BRANCHING )
               {
                  SCIPerrorMessage("debugging solution was cut off in local node %p at depth %d by inference <%s>[%.15g] %s %.15g\n",
                     node, SCIPnodeGetDepth(node), SCIPvarGetName(boundchgs[i].var), varsol,
                     SCIPboundchgGetBoundtype(&boundchgs[i]) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", boundchgs[i].newbound);
                  SCIPABORT();
               }
            }
            else if( SCIPboundchgGetBoundchgtype(&boundchgs[i]) == SCIP_BOUNDCHGTYPE_BRANCHING )
            {
               /* we branched on a variable were we don't know the solution: no debugging can be applied in this subtree */
               *solcontained = FALSE;
            }
         }
      }
   }

   /* remember the status of the current node */
   SCIP_CALL( SCIPhashmapSetImage(solinnode, (void*)node, *solcontained ? (void*)(&trueptr) : (void*)(&falseptr)) );

   return SCIP_OKAY;
}

/** frees all debugging solution data*/
SCIP_RETCODE SCIPdebugFreeDebugData(
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int s;

   assert(set != NULL);

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   for( s = nsolvals - 1; s >= 0; --s )
      BMSfreeMemoryArrayNull(&solnames[s]);

   BMSfreeMemoryArrayNull(&solnames);
   BMSfreeMemoryArrayNull(&solvals);

   nsolvals = 0;
   debugsolval = 0.0;
   mainscipset = NULL;
   solisachieved = FALSE;

   if( solinnode != NULL)
      SCIPhashmapFree(&solinnode);

   return SCIP_OKAY;
}

/** checks whether given row is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row                 /**< row to check for validity */
   )
{
   SCIP_COL** cols;
   SCIP_Real* vals;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nnonz;
   int i;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real solval;

   assert(set != NULL);
   assert(row != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* if the row is only locally valid, check whether the debugging solution is contained in the local subproblem */
   if( SCIProwIsLocal(row) )
   {
      SCIP_Bool solcontained;

      SCIP_CALL( isSolutionInNode(SCIPblkmem(set->scip), set, SCIPgetCurrentNode(set->scip), &solcontained) );
      if( !solcontained )
         return SCIP_OKAY;
   }

   cols = SCIProwGetCols(row);
   vals = SCIProwGetVals(row);
   nnonz = SCIProwGetNNonz(row);
   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);

   /* calculate row's activity on debugging solution */
   minactivity = SCIProwGetConstant(row);
   maxactivity = minactivity;
   for( i = 0; i < nnonz; ++i )
   {
      SCIP_VAR* var;

      /* get solution value of variable in debugging solution */
      var = SCIPcolGetVar(cols[i]);
      SCIP_CALL( getSolutionValue(set, var, &solval) );

      if( solval != SCIP_UNKNOWN ) /*lint !e777*/
      {
         minactivity += vals[i] * solval;
         maxactivity += vals[i] * solval;
      }
      else if( vals[i] > 0.0 )
      {
         minactivity += vals[i] * SCIPvarGetLbGlobal(var);
         maxactivity += vals[i] * SCIPvarGetUbGlobal(var);
      }
      else if( vals[i] < 0.0 )
      {
         minactivity += vals[i] * SCIPvarGetUbGlobal(var);
         maxactivity += vals[i] * SCIPvarGetLbGlobal(var);
      }
   }
   SCIPdebugMessage("debugging solution on row <%s>: %g <= [%g,%g] <= %g\n",
      SCIProwGetName(row), lhs, minactivity, maxactivity, rhs);

   /* check row for violation */
   if( SCIPsetIsFeasLT(set, maxactivity, lhs) || SCIPsetIsFeasGT(set, minactivity, rhs) )
   {
      printf("***** debug: row <%s> violates debugging solution (lhs=%.15g, rhs=%.15g, activity=[%.15g,%.15g], local=%d)\n",
         SCIProwGetName(row), lhs, rhs, minactivity, maxactivity, SCIProwIsLocal(row));
      SCIProwPrint(row, SCIPgetMessagehdlr(set->scip), NULL);

      /* output row with solution values */
      printf("\n\n");
      printf("***** debug: violated row <%s>:\n", SCIProwGetName(row));
      printf(" %.15g <= %.15g", lhs, SCIProwGetConstant(row));
      for( i = 0; i < nnonz; ++i )
      {
         /* get solution value of variable in debugging solution */
         SCIP_CALL( getSolutionValue(set, SCIPcolGetVar(cols[i]), &solval) );
         printf(" %+.15g<%s>[%.15g]", vals[i], SCIPvarGetName(SCIPcolGetVar(cols[i])), solval);
      }
      printf(" <= %.15g\n", rhs);

      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given global lower bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckLbGlobal(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lb                  /**< lower bound */
   )
{
   SCIP_Real varsol;

   assert(set != NULL);
   assert(var != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(set, var, &varsol) );
   SCIPdebugMessage("debugging solution on lower bound of <%s>[%g] >= %g\n", SCIPvarGetName(var), varsol, lb);

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN && SCIPsetIsFeasLT(set, varsol, lb) ) /*lint !e777*/
   {
      SCIPerrorMessage("invalid global lower bound: <%s>[%.15g] >= %.15g\n", SCIPvarGetName(var), varsol, lb);
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given global upper bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckUbGlobal(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   SCIP_Real varsol;

   assert(set != NULL);
   assert(var != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(set, var, &varsol) );
   SCIPdebugMessage("debugging solution on upper bound of <%s>[%g] <= %g\n", SCIPvarGetName(var), varsol, ub);

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN && SCIPsetIsFeasGT(set, varsol, ub) ) /*lint !e777*/
   {
      SCIPerrorMessage("invalid global upper bound: <%s>[%.15g] <= %.15g\n", SCIPvarGetName(var), varsol, ub);
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given local bound implication is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckInference(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< local node where this bound change was applied */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   SCIP_Real varsol;
   SCIP_Bool solcontained;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(var != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* in case we are in probing or diving we have to avoid checking the solution */
   if( SCIPlpDiving(set->scip->lp) || SCIPtreeProbing(set->scip->tree) )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(set, var, &varsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN ) /*lint !e777*/
   {
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasLT(set, varsol, newbound) )
      {
         SCIPerrorMessage("invalid local lower bound implication: <%s>[%.15g] >= %.15g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasGT(set, varsol, newbound) )
      {
         SCIPerrorMessage("invalid local upper bound implication: <%s>[%.15g] <= %.15g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** informs solution debugger, that the given node will be freed */
SCIP_RETCODE SCIPdebugRemoveNode(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node                /**< node that will be freed */
   )
{
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* check if a solution will be cutoff in tree */
   if( SCIPgetStage(set->scip) != SCIP_STAGE_EXITSOLVE && SCIPgetStage(set->scip) != SCIP_STAGE_EXITPRESOLVE && SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE )
   {
      SCIP_Bool solisinnode;

      solisinnode = FALSE;

      SCIP_CALL( isSolutionInNode(blkmem, set, node, &solisinnode) );
      /* wrong node will be cutoff */
      if( solisinnode )
      {
         SCIPerrorMessage("debugging solution was cut off in local node #%"SCIP_LONGINT_FORMAT" (%p) at depth %d\n",
            node->number, node, SCIPnodeGetDepth(node));
         SCIPABORT();
      }
   }

   /* remove node from the hash map */
   if( solinnode != NULL )
   {
      SCIP_CALL( SCIPhashmapRemove(solinnode, (void*)node) );
   }

   return SCIP_OKAY;
}

/** checks whether given variable bound is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckVbound(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable x in x <= b*z + d  or  x >= b*z + d */
   SCIP_BOUNDTYPE        vbtype,             /**< type of variable bound (LOWER or UPPER) */
   SCIP_VAR*             vbvar,              /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbcoef,             /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   SCIP_Real             vbconstant          /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   )
{
   SCIP_Real varsol;
   SCIP_Real vbvarsol;
   SCIP_Real vb;

   assert(set != NULL);
   assert(var != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* get solution value of variables */
   SCIP_CALL( getSolutionValue(set, var, &varsol) );
   SCIP_CALL( getSolutionValue(set, vbvar, &vbvarsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_UNKNOWN && vbvarsol != SCIP_UNKNOWN ) /*lint !e777*/
   {
      vb = vbcoef * vbvarsol + vbconstant;
      if( (vbtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasLT(set, varsol, vb))
         || (vbtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasGT(set, varsol, vb)) )
      {
         SCIPerrorMessage("invalid variable bound: <%s>[%.15g] %s %.15g<%s>[%.15g] %+.15g\n",
            SCIPvarGetName(var), varsol, vbtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", vbcoef,
            SCIPvarGetName(vbvar), vbvarsol, vbconstant);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** checks whether given implication is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckImplic(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound           /**< bound b    in implication y <= b or y >= b */
   )
{
   SCIP_Real solval;

   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(set, var, &solval) );
   if( solval == SCIP_UNKNOWN ) /*lint !e777*/
      return SCIP_OKAY;
   assert(SCIPsetIsFeasZero(set, solval) || SCIPsetIsFeasEQ(set, solval, 1.0));

   /* check, whether the implication applies for the debugging solution */
   if( (solval > 0.5) != varfixing )
      return SCIP_OKAY;

   /* get solution value of implied variable */
   SCIP_CALL( getSolutionValue(set, implvar, &solval) );
   if( solval == SCIP_UNKNOWN ) /*lint !e777*/
      return SCIP_OKAY;

   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      if( SCIPsetIsFeasLT(set, solval, implbound) )
      {
         SCIPerrorMessage("invalid implication <%s> == %d -> <%s> >= %.15g (variable has value %.15g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }
   else
   {
      if( SCIPsetIsFeasGT(set, solval, implbound) )
      {
         SCIPerrorMessage("invalid implication <%s> == %d -> <%s> <= %.15g (variable has value %.15g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** check whether given clique is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckClique(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< binary variables in the clique: at most one can be set to the given value */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars               /**< number of variables in the clique */
   )
{
   SCIP_Real solval;
   int pos1;
   int pos2;
   int v;

   assert(set != NULL);
   assert(vars != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   pos1 = -1;
   pos2 = -1;

   for( v = 0; v < nvars; ++v )
   {
      assert(vars[v] != NULL);
      assert(SCIPvarIsBinary(vars[v]));

      /* get solution value of variable */
      SCIP_CALL( getSolutionValue(set, vars[v], &solval) );

      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         continue;

      assert(SCIPsetIsFeasZero(set, solval) || SCIPsetIsFeasEQ(set, solval, 1.0));

      /* negated solution value if negated variable is in clique */
      if( values != NULL && values[v] == 0 )
         solval = 1.0 - solval;

      if( SCIPsetIsFeasEQ(set, solval, 1.0) )
      {
         if( pos1 == -1 )
            pos1 = v;
         else
         {
            assert(pos2 == -1);
            pos2 = v;
            break;
         }
      }
   }

   /* print debug message if the clique violates the debugging solution */
   if( pos2 != -1 )
   {
      assert(pos1 != -1);
      SCIPerrorMessage("clique violates debugging solution, (at least) variable <%s%s> and variable <%s%s> are both one in the debugging solution\n",
         (values == NULL || values[pos1]) ? "" : "~", SCIPvarGetName(vars[pos1]), (values == NULL || values[pos2]) ? "" : "~", SCIPvarGetName(vars[pos2]));
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given conflict is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckConflict(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node where the conflict clause is added */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change informations of the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos         /**< number of bound changes in the conflict set */
   )
{
   SCIP_Real solval;
   SCIP_Bool solcontained;
   int i;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(node != NULL);
   assert(nbdchginfos == 0 || bdchginfos != NULL);

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(set) )
      return SCIP_OKAY;

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* check, whether at least one literals is TRUE in the debugging solution */
   for( i = 0; i < nbdchginfos; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real newbound;

      var = SCIPbdchginfoGetVar(bdchginfos[i]);
      newbound = relaxedbds[i];

      SCIP_CALL( getSolutionValue(set, var, &solval) );
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
         return SCIP_OKAY;
      if( SCIPbdchginfoGetBoundtype(bdchginfos[i]) == SCIP_BOUNDTYPE_LOWER )
      {
         assert(SCIPsetIsLE(set, newbound, SCIPbdchginfoGetNewbound(bdchginfos[i])));

         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            if( SCIPsetIsLE(set, solval, newbound) )
               return SCIP_OKAY;
         }
         else
         {
            if( SCIPsetIsLT(set, solval, newbound) )
               return SCIP_OKAY;
         }
      }
      else
      {
         assert(SCIPsetIsGE(set, newbound, SCIPbdchginfoGetNewbound(bdchginfos[i])));

         if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
         {
            if( SCIPsetIsGE(set, solval, newbound) )
               return SCIP_OKAY;
         }
         else
         {
            if( SCIPsetIsGT(set, solval, newbound) )
               return SCIP_OKAY;
         }
      }
   }

   SCIPerrorMessage("invalid conflict set:");
   for( i = 0; i < nbdchginfos; ++i )
   {
      SCIP_CALL( getSolutionValue(set, SCIPbdchginfoGetVar(bdchginfos[i]), &solval) );
      printf(" <%s>[%.15g] %s %g(%g)", SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfos[i])), solval,
         SCIPbdchginfoGetBoundtype(bdchginfos[i]) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(bdchginfos[i]), relaxedbds[i]);
   }
   printf("\n");
   SCIPABORT();

   return SCIP_OKAY; /*lint !e527*/
}


/** check whether the debugging solution is valid in the current node */
SCIP_RETCODE SCIPdebugSolIsValidInSubtree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            isvalidinsubtree    /**< pointer to store whether the solution is valid in the current
                                              *   subtree
                                              */
   )
{
   SCIP_Bool solcontained;

   *isvalidinsubtree = FALSE;

   /* when debugging was disabled the solution is not defined to be not valid in the current subtree */
   if( debugsoldisabled )
      return SCIP_OKAY;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(scip->set) )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(scip->set) )
      return SCIP_OKAY;

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(SCIPblkmem(scip), scip->set, SCIPgetCurrentNode(scip), &solcontained) );

   if( solcontained )
      *isvalidinsubtree = TRUE;

   return SCIP_OKAY;
}


/** enabling solution debugging mechanism */
void SCIPdebugSolEnable(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   debugsoldisabled = FALSE;
}

/** disabling solution debugging mechanism */
void SCIPdebugSolDisable(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   debugsoldisabled = TRUE;
}

/** check if solution debugging mechanism is enabled */
SCIP_Bool SCIPdebugSolIsEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return (!debugsoldisabled);
}

/** propagator to force finding the debugging solution */
static
SCIP_DECL_PROPEXEC(propExecDebug)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   /* check if we are in the original problem and not in a sub MIP */
   if( !isSolutionInMip(scip->set) )
      return SCIP_OKAY;

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* check if the incumbent solution is at least as good as the debug solution, so we can stop to check the debug solution */
   if( debugSolIsAchieved(scip->set) )
      return SCIP_OKAY;

#if 1
   /* solve at least one LP */
   if( SCIPgetNLPIterations(scip) == 0 )
      return SCIP_OKAY;
#endif

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      SCIP_CALL( getSolutionValue(scip->set, vars[i], &solval) );
      if( solval == SCIP_UNKNOWN ) /*lint !e777*/
      {
         SCIPerrorMessage("original variable without debugging solution value\n");
         SCIPABORT();
      }

      lb = SCIPvarGetLbGlobal(vars[i]);
      ub = SCIPvarGetUbGlobal(vars[i]);
      if( SCIPisLT(scip, solval, lb) || SCIPisGT(scip, solval, ub) )
      {
         SCIPerrorMessage("solution value %.15g of <%s> outside bounds loc=[%.15g,%.15g], glb=[%.15g,%.15g]\n",
            solval, SCIPvarGetName(vars[i]), lb, ub, SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]));
         SCIPABORT();
      }

      SCIP_CALL( SCIPfixVar(scip, vars[i], solval, &infeasible, &fixed) );
      if( infeasible )
         *result = SCIP_CUTOFF;
      else if( fixed )
         *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

/** creates the debugging propagator and includes it in SCIP */
SCIP_RETCODE SCIPdebugIncludeProp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, "debug", "debugging propagator", 99999999, -1, FALSE,
         SCIP_PROPTIMING_ALWAYS, 99999999, 0, FALSE, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, propExecDebug, NULL, NULL) );

   return SCIP_OKAY;
}

/** adds a solution value for a new variable in the transformed problem that has no original counterpart
 * a value can only be set if no value has been set for this variable before
 */
SCIP_RETCODE SCIPdebugAddSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable for which to add a value */
   SCIP_Real             val                 /**< solution value for variable */
   )
{
   const char* varname;
   int i;

   assert(var != NULL);

   /* check if we are in the SCIP instance that we are debugging and not some different (subSCIP, auxiliary CIP, ...) */
   if( !isSolutionInMip(scip->set) )
      return SCIP_OKAY;

   if( SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("adding solution values for original variables is forbidden\n");
      return SCIP_ERROR;
   }

   if( SCIPvarIsTransformedOrigvar(var) )
   {
      SCIPerrorMessage("adding solution values for variable that are direct counterparts of original variables is forbidden\n");
      return SCIP_ERROR;
   }

   /* allocate memory */
   if( nsolvals >= solsize )
   {
      solsize = MAX(2*solsize, nsolvals+1);
      SCIP_ALLOC( BMSreallocMemoryArray(&solnames, solsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&solvals,  solsize) );
   }
   assert(nsolvals < solsize);

   /* store solution value in sorted list */
   varname = SCIPvarGetName(var);
   for( i = nsolvals; i > 0 && strcmp(varname, solnames[i-1]) < 0; --i )
   {
      solnames[i] = solnames[i-1];
      solvals[i]  = solvals[i-1];
   }
   if( i > 0 && strcmp(varname, solnames[i-1]) == 0 )
   {
      if( REALABS(solvals[i-1] - val) > 1e-9 )
      {
         SCIPerrorMessage("already have stored different debugging solution value (%g) for variable <%s>, cannot store %g\n", solvals[i-1], varname, val);
         return SCIP_ERROR;
      }
      else
      {
         SCIPdebugMessage("already have stored debugging solution value %g for variable <%s>, do not store same value again\n", val, varname);
         for( ; i < nsolvals; ++i )
         {
            solnames[i] = solnames[i+1];
            solvals[i]  = solvals[i+1];
         }
         return SCIP_OKAY;
      }
   }

   /* insert new solution value */
   SCIP_ALLOC( BMSduplicateMemoryArray(&solnames[i], varname, strlen(varname)+1) );
   SCIPdebugMessage("add variable <%s>: value <%g>\n", solnames[i], val);
   solvals[i] = val;
   nsolvals++;

   /* update objective function value of debug solution */
   debugsolval += solvals[i] * SCIPvarGetObj(var);
   SCIPdebugMessage("Debug Solution value is now %g.\n", debugsolval);

   return SCIP_OKAY;
}

#else

/** this is a dummy method to make the SunOS gcc linker happy */
extern void SCIPdummyDebugMethodForSun(void);
void SCIPdummyDebugMethodForSun(void)
{
      return;
}

#endif


/*
 * debug method for LP interface, to check if the LP interface works correct
 */
#ifdef SCIP_DEBUG_LP_INTERFACE

/* check whether coef is the r-th row of the inverse basis matrix B^-1; this is
 * the case if( coef * B ) is the r-th unit vector */
SCIP_RETCODE SCIPdebugCheckBInvRow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< r-th row of the inverse basis matrix */
   )
{
   SCIP_Real vecval;
   SCIP_Real matrixval;
   int* basisind;
   int nrows;
   int idx;
   int i;
   int k;

   assert(scip != NULL);

   nrows = SCIPgetNLPRows(scip);

   /* get basic indices for the basic matrix B */
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );


   /* loop over the columns of B */
   for( k = 0; k < nrows; ++k )
   {
      vecval = 0.0;

      /* indices of basic columns and rows:
       * - index i >= 0 corresponds to column i,
       * - index i < 0 to row -i-1
       */
      idx = basisind[k];

      /* check if we have a slack variable; this is the case if idx < 0 */
      if( idx >= 0 )
      {
         /* loop over the rows to compute the corresponding value in the unit vector */
         for( i = 0; i < nrows; ++i )
         {
            SCIP_CALL( SCIPlpiGetCoef(scip->lp->lpi, i, idx, &matrixval) );
            vecval += coef[i] * matrixval;
         }
      }
      else
      {
         assert( idx < 0 );

         /* retransform idx
          * - index i >= 0 corresponds to column i,
          * - index i < 0 to row -i-1
          */
         idx = -idx - 1;
         assert( idx >= 0 && idx < nrows );

         /* since idx < 0 we are in the case of a slack variable, i.e., the corresponding column
            is the idx-unit vector; note that some LP solver return a -idx-unit vector */
         /*   vecval = REALABS(coef[idx]);*/
         vecval = coef[idx];
      }

      /* check if vecval fits to the r-th unit vector */
      if( k == r && !SCIPisFeasEQ(scip, vecval, 1.0) )
      {
         /* we expected a 1.0 and found something different */
         SCIPmessagePrintWarning(SCIPgetMessagehdlr(scip), "checked SCIPgetLPBInvRow() found value <%g> expected 1.0\n", vecval);
      }
      else if( k != r && !SCIPisFeasZero(scip, vecval) )
      {
         /* we expected a 0.0 and found something different */
         SCIPmessagePrintWarning(SCIPgetMessagehdlr(scip), "checked SCIPgetLPBInvRow() found value <%g> expected 0.0\n", vecval);
      }
   }

   SCIPfreeBufferArray(scip, &basisind);

   return SCIP_OKAY;
}

#endif
