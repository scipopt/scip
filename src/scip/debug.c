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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: debug.c,v 1.2 2005/05/02 11:42:55 bzfpfend Exp $"

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
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/scip.h"
#include "scip/debug.h"


#ifdef DEBUG_SOLUTION

static char** solnames = NULL;
static Real* solvals = NULL;
static int nsolvals = 0;

/** reads feasible solution to check from file */
static
RETCODE readSolution(
   void
   )
{
   FILE* file;
   int solsize;

   if( nsolvals > 0 )
      return SCIP_OKAY;

   assert(solnames == NULL);
   assert(solvals == NULL);

   printf("***** debug: reading solution file <%s>\n", DEBUG_SOLUTION);

   /* open solution file */
   file = fopen(DEBUG_SOLUTION, "r");
   if( file == NULL )
   {
      errorMessage("cannot open solution file <%s> specified in scip/debug.h\n", DEBUG_SOLUTION);
      return SCIP_NOFILE;
   }

   /* read data */
   solsize = 0;
   while( !feof(file) )
   {
      char name[MAXSTRLEN];
      Real val;
      int nread;
      int i;

      nread = fscanf(file, "%s %lf\n", name, &val);
      if( nread != 2 )
      {
         printf("invalid input line %d in solution file <%s>\n", nsolvals, DEBUG_SOLUTION);
         fclose(file);
         return SCIP_READERROR;
      }

      /* allocate memory */
      if( nsolvals >= solsize )
      {
         solsize *= 2;
         solsize = MAX(solsize, nsolvals+1);
         ALLOC_OKAY( reallocMemoryArray(&solnames, solsize) );
         ALLOC_OKAY( reallocMemoryArray(&solvals, solsize) );         
      }
      assert(nsolvals < solsize);

      /* store solution value in sorted list */
      for( i = nsolvals; i > 0 && strcmp(name, solnames[i-1]) < 0; --i )
      {
         solnames[i] = solnames[i-1];
         solvals[i] = solvals[i-1];
      }
      ALLOC_OKAY( duplicateMemoryArray(&solnames[i], name, strlen(name)+1) );
      solvals[i] = val;
      nsolvals++;
   }

   /* close file */
   fclose(file);

   printf("***** debug: read %d non-zero entries\n", nsolvals);

   return SCIP_OKAY;
}

/** gets value of given variable in debugging solution */
static
RETCODE getSolutionValue(
   VAR*             var,                /**< variable to get solution value for */
   Real*            val                 /**< pointer to store solution value */
   )
{
   VAR* origvar;
   Real scalar;
   Real constant;
   const char* name;
   int left;
   int right;
   int middle;
   int cmp;

   assert(val != NULL);

   CHECK_OKAY( readSolution() );

   /* retransform variable onto orginal variable space */
   origvar = var;
   scalar = 1.0;
   constant = 0.0;
   CHECK_OKAY( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );
   if( origvar == NULL )
   {
      warningMessage("variable <%s> has no original counterpart\n", SCIPvarGetName(var));
      *val = 0.0;
      return SCIP_OKAY;
   }
   
   /* perform a binary search for the variable */
   name = SCIPvarGetName(origvar);
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

   return SCIP_OKAY;
}

/** checks whether given row is valid for the debugging solution */
RETCODE SCIPdebugCheckRow(
   ROW*             row,                /**< row to check for validity */
   SET*             set                 /**< global SCIP settings */
   )
{
   COL** cols;
   Real* vals;
   Real lhs;
   Real rhs;
   int nnonz;
   int i;
   Real activity;
   Real solval;

   cols = SCIProwGetCols(row);
   vals = SCIProwGetVals(row);
   nnonz = SCIProwGetNNonz(row);
   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);

   /* calculate row's activity on debugging solution */
   activity = SCIProwGetConstant(row);;
   for( i = 0; i < nnonz; ++i )
   {
      /* get solution value of variable in debugging solution */
      CHECK_OKAY( getSolutionValue(SCIPcolGetVar(cols[i]), &solval) );
      activity += vals[i] * solval;
   }

   /* check row for violation */
   if( SCIPsetIsFeasLT(set, activity, lhs) || SCIPsetIsFeasGT(set, activity, rhs) )
   {
      printf("***** debug: row <%s> violates debugging solution (lhs=%g, rhs=%g, activity=%g)\n",
         SCIProwGetName(row), lhs, rhs, activity);
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

      abort();
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** checks whether given implication is valid for the debugging solution */
RETCODE SCIPdebugCheckImplic(
   VAR*             var,                /**< problem variable  */
   SET*             set,                /**< global SCIP settings */
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
   assert(SCIPsetIsEQ(set, solval, 0.0) || SCIPsetIsEQ(set, solval, 1.0));

   /* check, whether the implication applies for the debugging solution */
   if( (solval > 0.5) != varfixing )
      return SCIP_OKAY;

   /* get solution value of implied variable */
   CHECK_OKAY( getSolutionValue(implvar, &solval) );
   if( impltype == SCIP_BOUNDTYPE_LOWER )
   {
      if( SCIPsetIsLT(set, solval, implbound) )
      {
         errorMessage("invalid implication <%s> == %d -> <%s> >= %g (variable has value %g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         abort();
         return SCIP_ERROR;
      }
   }
   else
   {
      if( SCIPsetIsGT(set, solval, implbound) )
      {
         errorMessage("invalid implication <%s> == %d -> <%s> <= %g (variable has value %g in solution)\n",
            SCIPvarGetName(var), varfixing, SCIPvarGetName(implvar), implbound, solval);
         abort();
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

#endif
