/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  main file for C compilation
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "scip.h"
#include "cons_linear.h"
#include "nodesel_dfs.h"


static const int   nrows = 2;
static const int   nvars = 2;

static const char* var_name[] = { "x1"    , "x2"     };
static const Real  var_obj [] = { -1.0    , -1.0     };
static const Real  var_lb  [] = {  0.0    ,  0.0     };
static const Real  var_ub  [] = { 10.0    ,  5.0     };

static const char* row_name[] = { "lin1"  , "lin2"   };
static const int   row_len [] = {        2,        2 };
static const int   row_idx [] = {   0,   1,   0,   1 };
static const Real  row_val [] = { 2.0, 1.0, 1.0, 2.0 };
static const Real  row_rhs [] = {      2.0,      2.0 };

int
main(void)
{
   SCIP* scip = NULL;
   RETCODE retcode;
   VAR** vars;
   VAR** rowvars;
   Real* rowvals;
   int r;
   int v;
   int i;
   int pos;

   printf("SCIP version %g\n", SCIPversion());

   /***************************
    * Local Memory Allocation *
    ***************************/

   allocMemoryArray(vars, nvars);
   allocMemoryArray(rowvars, nvars);
   allocMemoryArray(rowvals, nvars);


   /*********
    * Setup *
    *********/

   printf("\nsetup SCIP\n");

   /* initialize SCIP */
   CHECK_SCIP( SCIPcreate(&scip) );

   /* change settings */
   CHECK_SCIP( SCIPsetVerbLevel(scip, SCIP_VERBLEVEL_FULL) );

   /* include user defined callbacks */
   CHECK_SCIP( SCIPincludeConsHdlrLinear(scip) );
   CHECK_SCIP( SCIPincludeNodeselDfs(scip) );


   /********************
    * Problem Creation *
    ********************/

   printf("\ncreate problem\n");

   /* create problem */
   CHECK_SCIP( SCIPcreateProb(scip, "test.lp") );

   /* create necessary variables */
   for( v = 0; v < nvars; ++v )
   {
      CHECK_SCIP( SCIPcreateVar(scip, &vars[v], var_name[v], var_lb[v], var_ub[v], var_obj[v], SCIP_VARTYPE_INTEGER) );
      CHECK_SCIP( SCIPaddVar(scip, vars[v]) );
   }

   /* create constraints and add them to the problem */
   pos = 0;
   for( r = 0; r < nrows; ++r )
   {
      CONS* cons;

      for( i = 0; i < row_len[r]; ++i )
      {
         rowvars[i] = vars[row_idx[pos]];
         rowvals[i] = row_val[pos];
         pos++;
      }
      CHECK_SCIP( SCIPcreateConsLinear(scip, &cons, row_name[r], row_len[r], rowvars, rowvals, 0.0, row_rhs[r], TRUE) );
      CHECK_SCIP( SCIPaddCons(scip, cons) ); /* add as a global constraint */
   }


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   printf("\nsolve problem\n");
   CHECK_SCIP( SCIPsolve(scip) );

   /* free solution process */
   printf("\nfree problem solution\n");
   CHECK_SCIP( SCIPfreeSolve(scip) );

   /* solve problem again */
   printf("\nsolve problem again\n");
   CHECK_SCIP( SCIPsolve(scip) );

#ifndef NDEBUG
   /*SCIPdebugMemory(scip);*/
#endif


   /********************
    * Deinitialization *
    ********************/

   printf("\nfree SCIP\n");

   /* free SCIP */
   CHECK_SCIP( SCIPfree(&scip) );


   /*****************************
    * Local Memory Deallocation *
    *****************************/

   freeMemoryArray(rowvals);
   freeMemoryArray(rowvars);
   freeMemoryArray(vars);

#ifndef NDEBUG
   memoryCheckEmpty();
#endif

   return 0;
}

