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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: debug.c,v 1.13 2005/10/18 15:21:27 bzfpfend Exp $"

/**@file   debug.c
 * @brief  methods for debugging
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "scip/def.h"
#include "scip/message.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/misc.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/debug.h"


#ifdef SCIP_DEBUG_SOLUTION

/** reads solution from given file into given arrays */
static
SCIP_RETCODE readSolfile(
   const char*           solfilename,        /**< solution filename to read */
   char***               names,              /**< pointer to store the array of variable names */
   SCIP_Real**           vals,               /**< pointer to store the array of solution values */
   int*                  nvals               /**< pointer to store the number of non-zero elements */
   )
{
   FILE* file;
   int solsize;

   assert(*names == NULL);
   assert(*vals == NULL);

   printf("***** debug: reading solution file <%s>\n", solfilename);

   /* open solution file */
   file = fopen(solfilename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open solution file <%s> specified in scip/debug.h\n", solfilename);
      return SCIP_NOFILE;
   }

   /* read data */
   solsize = 0;
   while( !feof(file) )
   {
      char buf[SCIP_MAXSTRLEN];
      char name[SCIP_MAXSTRLEN];
      char objstring[SCIP_MAXSTRLEN];
      SCIP_Real val;
      int nread;
      int i;

      fgets(buf, SCIP_MAXSTRLEN, file);

      /* the lines "solution status: ..." and "objective value: ..." may preceed the solution information */
      if( strncmp(buf, "solution", 8) == 0 || strncmp(buf, "objective", 9) == 0 )
         continue;

      nread = sscanf(buf, "%s %lf %s\n", name, &val, objstring);
      if( nread < 2 )
      {
         printf("invalid input line %d in solution file <%s>: <%s>\n", *nvals, SCIP_DEBUG_SOLUTION, name);
         fclose(file);
         return SCIP_READERROR;
      }

      /* allocate memory */
      if( *nvals >= solsize )
      {
         solsize *= 2;
         solsize = MAX(solsize, (*nvals)+1);
         SCIP_ALLOC( BMSreallocMemoryArray(names, solsize) );
         SCIP_ALLOC( BMSreallocMemoryArray(vals, solsize) );         
      }
      assert(*nvals < solsize);

      /* store solution value in sorted list */
      for( i = *nvals; i > 0 && strcmp(name, (*names)[i-1]) < 0; --i )
      {
         (*names)[i] = (*names)[i-1];
         (*vals)[i] = (*vals)[i-1];
      }
      SCIP_ALLOC( BMSduplicateMemoryArray(&(*names)[i], name, strlen(name)+1) );
      (*vals)[i] = val;
      (*nvals)++;
   }

   /* close file */
   fclose(file);

   printf("***** debug: read %d non-zero entries\n", *nvals);

   return SCIP_OKAY;
}

static char** solnames = NULL;
static SCIP_Real* solvals = NULL;
static int nsolvals = 0;

/** reads feasible solution to check from file */
static
SCIP_RETCODE readSolution(
   void
   )
{
   if( nsolvals > 0 )
      return SCIP_OKAY;

   SCIP_CALL( readSolfile(SCIP_DEBUG_SOLUTION, &solnames, &solvals, &nsolvals) );

   return SCIP_OKAY;
}

/** gets value of given variable in debugging solution */
static
SCIP_RETCODE getSolutionValue(
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

   assert(val != NULL);

   SCIP_CALL( readSolution() );

   /* ignore deleted variables */
   if( SCIPvarIsDeleted(var) )
   {
      SCIPdebugMessage("**** invalid solution value for deleted variable <%s>\n", SCIPvarGetName(var));
      *val = SCIP_INVALID;
      return SCIP_OKAY;
   }

   /* retransform variable onto orginal variable space */
   solvar = var;
   scalar = 1.0;
   constant = 0.0;
   if( SCIPvarIsNegated(solvar) )
   {
      scalar = -1.0;
      constant = SCIPvarGetNegationConstant(solvar);
      solvar = SCIPvarGetNegationVar(solvar);
   }
   if( SCIPvarIsTransformedOrigvar(solvar) )
   {
      SCIP_CALL( SCIPvarGetOrigvarSum(&solvar, &scalar, &constant) );
      if( solvar == NULL )
      {
         SCIPwarningMessage("variable <%s> has no original counterpart\n", SCIPvarGetName(var));
         *val = SCIP_INVALID;
         return SCIP_OKAY;
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
      SCIPwarningMessage("invalid solution value %g for variable <%s>[%g,%g]\n",
         *val, SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   return SCIP_OKAY;
}

static SCIP_HASHMAP* solinnode = NULL;       /**< maps nodes to bools, storing whether the solution is valid for the node */
static SCIP_Bool falseptr = FALSE;
static SCIP_Bool trueptr = TRUE;

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

   /* generate the hashmap */
   if( solinnode == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&solinnode, blkmem, 107377) );
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

   /* if the solution is not contained in the parent of the node, it cannot be containt in the current node */
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
         for( i = 0; i < domchgbound->nboundchgs && *solcontained; ++i )
         {
            SCIP_Real varsol;

            /* get solution value of variable */
            SCIP_CALL( getSolutionValue(boundchgs[i].var, &varsol) );

            /* compare the bound change with the solution value */
            if( boundchgs[i].boundtype == SCIP_BOUNDTYPE_LOWER )
               *solcontained = SCIPsetIsFeasGE(set, varsol, boundchgs[i].newbound);
            else
               *solcontained = SCIPsetIsFeasLE(set, varsol, boundchgs[i].newbound);
            if( !(*solcontained) && boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING )
            {
               SCIPerrorMessage("debugging solution was cut off in local node %p at depth %d by inference <%s>[%.8g] %s %.8g\n",
                  node, SCIPnodeGetDepth(node), SCIPvarGetName(boundchgs[i].var), varsol,
                  boundchgs[i].boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", boundchgs[i].newbound);
               SCIPABORT();
            }
         }
      }
   }

   /* remember the status of the current node */
   SCIP_CALL( SCIPhashmapSetImage(solinnode, (void*)node, *solcontained ? (void*)(&trueptr) : (void*)(&falseptr)) );
   
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
      SCIP_CALL( getSolutionValue(var, &solval) );
      if( solval != SCIP_INVALID )
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
      printf("***** debug: row <%s> violates debugging solution (lhs=%g, rhs=%g, activity=[%g,%g], local=%d)\n",
         SCIProwGetName(row), lhs, rhs, minactivity, maxactivity, SCIProwIsLocal(row));
      SCIProwPrint(row, NULL);

      /* output row with solution values */
      printf("\n\n");
      printf("***** debug: violated row <%s>:\n", SCIProwGetName(row));
      printf(" %g <= %g", lhs, SCIProwGetConstant(row));
      for( i = 0; i < nnonz; ++i )
      {
         /* get solution value of variable in debugging solution */
         SCIP_CALL( getSolutionValue(SCIPcolGetVar(cols[i]), &solval) );
         printf(" %+g<%s>[%g]", vals[i], SCIPvarGetName(SCIPcolGetVar(cols[i])), solval);
      }
      printf(" <= %g\n", rhs);

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

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(var, &varsol) );
   SCIPdebugMessage("debugging solution on lower bound of <%s>[%g] >= %g\n", SCIPvarGetName(var), varsol, lb);

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID && SCIPsetIsLT(set, varsol, lb) )
   {
      SCIPerrorMessage("invalid global lower bound: <%s>[%.8g] >= %.8g\n", SCIPvarGetName(var), varsol, lb);
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

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(var, &varsol) );
   SCIPdebugMessage("debugging solution on upper bound of <%s>[%g] <= %g\n", SCIPvarGetName(var), varsol, ub);

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID && SCIPsetIsGT(set, varsol, ub) )
   {
      SCIPerrorMessage("invalid global upper bound: <%s>[%.8g] <= %.8g\n", SCIPvarGetName(var), varsol, ub);
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

   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(var, &varsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID )
   {
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLT(set, varsol, newbound) )
      {
         SCIPerrorMessage("invalid local lower bound implication: <%s>[%.8g] >= %.8g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGT(set, varsol, newbound) )
      {
         SCIPerrorMessage("invalid local upper bound implication: <%s>[%.8g] <= %.8g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** informs solution debugger, that the given node will be freed */
SCIP_RETCODE SCIPdebugRemoveNode(
   SCIP_NODE*            node                /**< node that will be freed */
   )
{
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

   /* get solution value of variables */
   SCIP_CALL( getSolutionValue(var, &varsol) );
   SCIP_CALL( getSolutionValue(vbvar, &vbvarsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID && vbvarsol != SCIP_INVALID )
   {
      vb = vbcoef * vbvarsol + vbconstant;
      if( (vbtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLT(set, varsol, vb))
         || (vbtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGT(set, varsol, vb)) )
      {
         SCIPerrorMessage("invalid variable bound: <%s>[%.8g] %s %.8g<%s>[%g] %+.8g\n", 
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

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   /* get solution value of variable */
   SCIP_CALL( getSolutionValue(var, &solval) );
   if( solval == SCIP_INVALID )
      return SCIP_OKAY;
   assert(SCIPsetIsFeasEQ(set, solval, 0.0) || SCIPsetIsFeasEQ(set, solval, 1.0));

   /* check, whether the implication applies for the debugging solution */
   if( (solval > 0.5) != varfixing )
      return SCIP_OKAY;

   /* get solution value of implied variable */
   SCIP_CALL( getSolutionValue(implvar, &solval) );
   if( solval == SCIP_INVALID )
      return SCIP_OKAY;

   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      if( SCIPsetIsLT(set, solval, implbound) )
      {
         SCIPerrorMessage("invalid implication <%s> == %d -> <%s> >= %.8g (variable has value %g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }
   else
   {
      if( SCIPsetIsGT(set, solval, implbound) )
      {
         SCIPerrorMessage("invalid implication <%s> == %d -> <%s> <= %.8g (variable has value %g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** checks whether given conflict is valid for the debugging solution */
SCIP_RETCODE SCIPdebugCheckConflict(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node where the conflict clause is added */
   SCIP_VAR**            conflictset,        /**< variables in the conflict set */
   int                   nliterals           /**< number of literals in the conflict set */
   )
{
   SCIP_Real solval;
   SCIP_Bool solcontained;
   int i;

   assert(conflictset != NULL);
   
   /* check whether the debugging solution is contained in the local subproblem */
   SCIP_CALL( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* check, whether at least one literals is TRUE in the debugging solution */
   for( i = 0; i < nliterals; ++i )
   {
      if( SCIPvarGetType(conflictset[i]) != SCIP_VARTYPE_BINARY )
      {
         SCIPerrorMessage("non-binary variable <%s> in conflict set\n", SCIPvarGetName(conflictset[i]));
         SCIPABORT();
      }

      SCIP_CALL( getSolutionValue(conflictset[i], &solval) );
      if( solval == SCIP_INVALID )
         return SCIP_OKAY;
      if( solval > 0.5 )
         return SCIP_OKAY;
   }

   SCIPerrorMessage("invalid conflict clause:");
   for( i = 0; i < nliterals; ++i )
   {
      SCIP_CALL( getSolutionValue(conflictset[i], &solval) );
      printf(" <%s>[%g]", SCIPvarGetName(conflictset[i]), solval);
   }
   printf("\n");
   SCIPABORT();

   return SCIP_OKAY;
}


/** propagator to force finding the debugging solution */
static
SCIP_DECL_PROPEXEC(propExecDebug)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nvars;
   int i;

   *result = SCIP_DIDNOTFIND;
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
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

      SCIP_CALL( getSolutionValue(vars[i], &solval) );
      if( solval == SCIP_INVALID )
      {
         SCIPerrorMessage("original variable without debugging solution value\n");
         SCIPABORT();
      }

      lb = SCIPvarGetLbGlobal(vars[i]);
      ub = SCIPvarGetUbGlobal(vars[i]);
      if( SCIPisLT(scip, solval, lb) || SCIPisGT(scip, solval, ub) )
      {
         SCIPerrorMessage("solution value %g of <%s> outside bounds loc=[%.8g,%.8g], glb=[%.8g,%.8g]\n",
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
   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, "debug", "debugging propagator", 99999999, -1, FALSE,
         NULL, NULL, NULL, NULL, NULL, propExecDebug, NULL, NULL) );

   return SCIP_OKAY;
}

#endif
