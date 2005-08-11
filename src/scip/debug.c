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
#pragma ident "@(#) $Id: debug.c,v 1.8 2005/08/11 09:59:27 bzfpfend Exp $"

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
#include "scip/memory.h"
#include "scip/set.h"
#include "scip/misc.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/debug.h"


#ifdef DEBUG_SOLUTION

/** reads solution from given file into given arrays */
static
RETCODE readSolfile(
   const char*      solfilename,        /**< solution filename to read */
   char***          names,              /**< pointer to store the array of variable names */
   Real**           vals,               /**< pointer to store the array of solution values */
   int*             nvals               /**< pointer to store the number of non-zero elements */
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
      errorMessage("cannot open solution file <%s> specified in scip/debug.h\n", solfilename);
      return SCIP_NOFILE;
   }

   /* read data */
   solsize = 0;
   while( !feof(file) )
   {
      char buf[MAXSTRLEN];
      char name[MAXSTRLEN];
      char objstring[MAXSTRLEN];
      Real val;
      int nread;
      int i;

      fgets(buf, MAXSTRLEN, file);

      /* the lines "solution status: ..." and "objective value: ..." may preceed the solution information */
      if( strncmp(buf, "solution", 8) == 0 || strncmp(buf, "objective", 9) == 0 )
         continue;

      nread = sscanf(buf, "%s %lf %s\n", name, &val, objstring);
      if( nread < 2 )
      {
         printf("invalid input line %d in solution file <%s>: <%s>\n", *nvals, DEBUG_SOLUTION, name);
         fclose(file);
         return SCIP_READERROR;
      }

      /* allocate memory */
      if( *nvals >= solsize )
      {
         solsize *= 2;
         solsize = MAX(solsize, (*nvals)+1);
         ALLOC_OKAY( reallocMemoryArray(names, solsize) );
         ALLOC_OKAY( reallocMemoryArray(vals, solsize) );         
      }
      assert(*nvals < solsize);

      /* store solution value in sorted list */
      for( i = *nvals; i > 0 && strcmp(name, (*names)[i-1]) < 0; --i )
      {
         (*names)[i] = (*names)[i-1];
         (*vals)[i] = (*vals)[i-1];
      }
      ALLOC_OKAY( duplicateMemoryArray(&(*names)[i], name, strlen(name)+1) );
      (*vals)[i] = val;
      (*nvals)++;
   }

   /* close file */
   fclose(file);

   printf("***** debug: read %d non-zero entries\n", *nvals);

   return SCIP_OKAY;
}

static char** solnames = NULL;
static Real* solvals = NULL;
static int nsolvals = 0;

/** reads feasible solution to check from file */
static
RETCODE readSolution(
   void
   )
{
   if( nsolvals > 0 )
      return SCIP_OKAY;

   CHECK_OKAY( readSolfile(DEBUG_SOLUTION, &solnames, &solvals, &nsolvals) );

   return SCIP_OKAY;
}

/** gets value of given variable in debugging solution */
static
RETCODE getSolutionValue(
   VAR*             var,                /**< variable to get solution value for */
   Real*            val                 /**< pointer to store solution value */
   )
{
   VAR* solvar;
   Real scalar;
   Real constant;
   const char* name;
   int left;
   int right;
   int middle;
   int cmp;

   assert(val != NULL);

   CHECK_OKAY( readSolution() );

   /* ignore deleted variables */
   if( SCIPvarIsDeleted(var) )
   {
      debugMessage("**** invalid solution value for deleted variable <%s>\n", SCIPvarGetName(var));
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
      CHECK_OKAY( SCIPvarGetOrigvarSum(&solvar, &scalar, &constant) );
      if( solvar == NULL )
      {
         warningMessage("variable <%s> has no original counterpart\n", SCIPvarGetName(var));
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
      warningMessage("invalid solution value %g for variable <%s>[%g,%g]\n",
         *val, SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));
   }

   return SCIP_OKAY;
}

static HASHMAP* solinnode = NULL;       /**< maps nodes to bools, storing whether the solution is valid for the node */
static Bool falseptr = FALSE;
static Bool trueptr = TRUE;

/** returns whether the solution is contained in node's subproblem */
static
RETCODE isSolutionInNode(
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   NODE*            node,               /**< local node where this bound change was applied */
   Bool*            solcontained        /**< pointer to store whether the solution is contained in node's subproblem */
   )
{
   Bool* boolptr;

   /* generate the hashmap */
   if( solinnode == NULL )
   {
      CHECK_OKAY( SCIPhashmapCreate(&solinnode, blkmem, 107377) );
   }

   /* check, whether we know already whether the solution is contained in the given node */
   boolptr = (Bool*)SCIPhashmapGetImage(solinnode, (void*)node);
   if( boolptr != NULL )
   {
      if( boolptr != &falseptr && boolptr != &trueptr )
      {
         errorMessage("wrong value in node hashmap\n");
         SCIPABORT();
      }
      *solcontained = *boolptr;
      return SCIP_OKAY;
   }

   /* if the solution is not contained in the parent of the node, it cannot be containt in the current node */
   *solcontained = TRUE;
   if( node->parent != NULL )
   {
      CHECK_OKAY( isSolutionInNode(blkmem, set, node->parent, solcontained) );
   }

   if( *solcontained )
   {
      /* check whether the bound changes at the current node remove the debugging solution from the subproblem */
      if( node->domchg != NULL )
      {
         DOMCHGBOUND* domchgbound;
         BOUNDCHG* boundchgs;
         int i;

         domchgbound = &node->domchg->domchgbound;
         boundchgs = domchgbound->boundchgs;
         for( i = 0; i < domchgbound->nboundchgs && *solcontained; ++i )
         {
            Real varsol;

            /* get solution value of variable */
            CHECK_OKAY( getSolutionValue(boundchgs[i].var, &varsol) );

            /* compare the bound change with the solution value */
            if( boundchgs[i].boundtype == SCIP_BOUNDTYPE_LOWER )
               *solcontained = SCIPsetIsFeasGE(set, varsol, boundchgs[i].newbound);
            else
               *solcontained = SCIPsetIsFeasLE(set, varsol, boundchgs[i].newbound);
            if( !(*solcontained) && boundchgs[i].boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING )
            {
               errorMessage("debugging solution was cut off in local node %p at depth %d by inference <%s>[%.8g] %s %.8g\n",
                  node, SCIPnodeGetDepth(node), SCIPvarGetName(boundchgs[i].var), varsol,
                  boundchgs[i].boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", boundchgs[i].newbound);
               SCIPABORT();
            }
         }
      }
   }

   /* remember the status of the current node */
   CHECK_OKAY( SCIPhashmapSetImage(solinnode, (void*)node, *solcontained ? (void*)(&trueptr) : (void*)(&falseptr)) );
   
   return SCIP_OKAY;
}

/** checks whether given row is valid for the debugging solution */
RETCODE SCIPdebugCheckRow(
   SET*             set,                /**< global SCIP settings */
   ROW*             row                 /**< row to check for validity */
   )
{
   COL** cols;
   Real* vals;
   Real lhs;
   Real rhs;
   int nnonz;
   int i;
   Real minactivity;
   Real maxactivity;
   Real solval;

   /* if the row is only locally valid, check whether the debugging solution is contained in the local subproblem */
   if( SCIProwIsLocal(row) )
   {
      Bool solcontained;

      CHECK_OKAY( isSolutionInNode(SCIPblkmem(set->scip), set, SCIPgetCurrentNode(set->scip), &solcontained) );
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
      VAR* var;

      /* get solution value of variable in debugging solution */
      var = SCIPcolGetVar(cols[i]);
      CHECK_OKAY( getSolutionValue(var, &solval) );
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
   debugMessage("debugging solution on row <%s>: %g <= [%g,%g] <= %g\n", 
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
         CHECK_OKAY( getSolutionValue(SCIPcolGetVar(cols[i]), &solval) );
         printf(" %+g<%s>[%g]", vals[i], SCIPvarGetName(SCIPcolGetVar(cols[i])), solval);
      }
      printf(" <= %g\n", rhs);

      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given global lower bound is valid for the debugging solution */
RETCODE SCIPdebugCheckLbGlobal(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             lb                  /**< lower bound */
   )
{
   Real varsol;

   /* get solution value of variable */
   CHECK_OKAY( getSolutionValue(var, &varsol) );
   debugMessage("debugging solution on lower bound of <%s>[%g] >= %g\n", SCIPvarGetName(var), varsol, lb);

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID && SCIPsetIsLT(set, varsol, lb) )
   {
      errorMessage("invalid global lower bound: <%s>[%.8g] >= %.8g\n", SCIPvarGetName(var), varsol, lb);
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given global upper bound is valid for the debugging solution */
RETCODE SCIPdebugCheckUbGlobal(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Real             ub                  /**< upper bound */
   )
{
   Real varsol;

   /* get solution value of variable */
   CHECK_OKAY( getSolutionValue(var, &varsol) );
   debugMessage("debugging solution on upper bound of <%s>[%g] <= %g\n", SCIPvarGetName(var), varsol, ub);

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID && SCIPsetIsGT(set, varsol, ub) )
   {
      errorMessage("invalid global upper bound: <%s>[%.8g] <= %.8g\n", SCIPvarGetName(var), varsol, ub);
      SCIPABORT();
   }

   return SCIP_OKAY;
}

/** checks whether given local bound implication is valid for the debugging solution */
RETCODE SCIPdebugCheckInference(
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   NODE*            node,               /**< local node where this bound change was applied */
   VAR*             var,                /**< problem variable */
   Real             newbound,           /**< new value for bound */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   )
{
   Real varsol;
   Bool solcontained;

   /* check whether the debugging solution is contained in the local subproblem */
   CHECK_OKAY( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* get solution value of variable */
   CHECK_OKAY( getSolutionValue(var, &varsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID )
   {
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLT(set, varsol, newbound) )
      {
         errorMessage("invalid local lower bound implication: <%s>[%.8g] >= %.8g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGT(set, varsol, newbound) )
      {
         errorMessage("invalid local upper bound implication: <%s>[%.8g] <= %.8g\n", SCIPvarGetName(var), varsol, newbound);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** informs solution debugger, that the given node will be freed */
RETCODE SCIPdebugRemoveNode(
   NODE*            node                /**< node that will be freed */
   )
{
   /* remove node from the hash map */
   if( solinnode != NULL )
   {
      CHECK_OKAY( SCIPhashmapRemove(solinnode, (void*)node) );
   }

   return SCIP_OKAY;
}

/** checks whether given variable bound is valid for the debugging solution */
RETCODE SCIPdebugCheckVbound(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable x in x <= b*z + d  or  x >= b*z + d */
   BOUNDTYPE        vbtype,             /**< type of variable bound (LOWER or UPPER) */
   VAR*             vbvar,              /**< variable z    in x <= b*z + d  or  x >= b*z + d */
   Real             vbcoef,             /**< coefficient b in x <= b*z + d  or  x >= b*z + d */
   Real             vbconstant          /**< constant d    in x <= b*z + d  or  x >= b*z + d */
   )
{
   Real varsol;
   Real vbvarsol;
   Real vb;

   /* get solution value of variables */
   CHECK_OKAY( getSolutionValue(var, &varsol) );
   CHECK_OKAY( getSolutionValue(vbvar, &vbvarsol) );

   /* check validity of debugging solution */
   if( varsol != SCIP_INVALID && vbvarsol != SCIP_INVALID )
   {
      vb = vbcoef * vbvarsol + vbconstant;
      if( (vbtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsLT(set, varsol, vb))
         || (vbtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsGT(set, varsol, vb)) )
      {
         errorMessage("invalid variable bound: <%s>[%.8g] %s %.8g<%s>[%g] %+.8g\n", 
            SCIPvarGetName(var), varsol, vbtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", vbcoef,
            SCIPvarGetName(vbvar), vbvarsol, vbconstant);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** checks whether given implication is valid for the debugging solution */
RETCODE SCIPdebugCheckImplic(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< problem variable */
   Bool             varfixing,          /**< FALSE if y should be added in implications for x == 0, TRUE for x == 1 */
   VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER) or y >= b (SCIP_BOUNDTYPE_LOWER) */
   Real             implbound           /**< bound b    in implication y <= b or y >= b */
   )
{
   Real solval;

   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

   /* get solution value of variable */
   CHECK_OKAY( getSolutionValue(var, &solval) );
   if( solval == SCIP_INVALID )
      return SCIP_OKAY;
   assert(SCIPsetIsEQ(set, solval, 0.0) || SCIPsetIsEQ(set, solval, 1.0));

   /* check, whether the implication applies for the debugging solution */
   if( (solval > 0.5) != varfixing )
      return SCIP_OKAY;

   /* get solution value of implied variable */
   CHECK_OKAY( getSolutionValue(implvar, &solval) );
   if( solval == SCIP_INVALID )
      return SCIP_OKAY;

   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      if( SCIPsetIsLT(set, solval, implbound) )
      {
         errorMessage("invalid implication <%s> == %d -> <%s> >= %.8g (variable has value %g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }
   else
   {
      if( SCIPsetIsGT(set, solval, implbound) )
      {
         errorMessage("invalid implication <%s> == %d -> <%s> <= %.8g (variable has value %g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         SCIPABORT();
      }
   }

   return SCIP_OKAY;
}

/** checks whether given conflict is valid for the debugging solution */
RETCODE SCIPdebugCheckConflict(
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   NODE*            node,               /**< node where the conflict clause is added */
   VAR**            conflictset,        /**< variables in the conflict set */
   int              nliterals           /**< number of literals in the conflict set */
   )
{
   Real solval;
   Bool solcontained;
   int i;

   assert(conflictset != NULL);
   
   /* check whether the debugging solution is contained in the local subproblem */
   CHECK_OKAY( isSolutionInNode(blkmem, set, node, &solcontained) );
   if( !solcontained )
      return SCIP_OKAY;

   /* check, whether at least one literals is TRUE in the debugging solution */
   for( i = 0; i < nliterals; ++i )
   {
      if( SCIPvarGetType(conflictset[i]) != SCIP_VARTYPE_BINARY )
      {
         errorMessage("non-binary variable <%s> in conflict set\n", SCIPvarGetName(conflictset[i]));
         SCIPABORT();
      }

      CHECK_OKAY( getSolutionValue(conflictset[i], &solval) );
      if( solval == SCIP_INVALID )
         return SCIP_OKAY;
      if( solval > 0.5 )
         return SCIP_OKAY;
   }

   errorMessage("invalid conflict clause:");
   for( i = 0; i < nliterals; ++i )
   {
      CHECK_OKAY( getSolutionValue(conflictset[i], &solval) );
      printf(" <%s>[%g]", SCIPvarGetName(conflictset[i]), solval);
   }
   printf("\n");
   SCIPABORT();

   return SCIP_OKAY;
}


/** propagator to force finding the debugging solution */
static
DECL_PROPEXEC(propExecDebug)
{  /*lint --e{715}*/
   VAR** vars;
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
      Real solval;
      Real lb;
      Real ub;
      Bool infeasible;
      Bool fixed;

      CHECK_OKAY( getSolutionValue(vars[i], &solval) );
      if( solval == SCIP_INVALID )
      {
         errorMessage("original variable without debugging solution value\n");
         SCIPABORT();
      }

      lb = SCIPvarGetLbGlobal(vars[i]);
      ub = SCIPvarGetUbGlobal(vars[i]);
      if( SCIPisLT(scip, solval, lb) || SCIPisGT(scip, solval, ub) )
      {
         errorMessage("solution value %g of <%s> outside bounds loc=[%.8g,%.8g], glb=[%.8g,%.8g]\n",
            solval, SCIPvarGetName(vars[i]), lb, ub, SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]));
         SCIPABORT();
      }
      
      CHECK_OKAY( SCIPfixVar(scip, vars[i], solval, &infeasible, &fixed) );
      if( infeasible )
         *result = SCIP_CUTOFF;
      else if( fixed )
         *result = SCIP_REDUCEDDOM;
   }

   return SCIP_OKAY;
}

/** creates the debugging propagator and includes it in SCIP */
RETCODE SCIPdebugIncludeProp(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   /* include propagator */
   CHECK_OKAY( SCIPincludeProp(scip, "debug", "debugging propagator", 99999999, -1, FALSE,
         NULL, NULL, NULL, NULL, NULL, propExecDebug, NULL, NULL) );

   return SCIP_OKAY;
}

#endif
