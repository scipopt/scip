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
#include "reader_mps.h"
#include "disp_default.h"
#include "cons_integral.h"
#include "cons_linear.h"
#include "nodesel_bfs.h"
#include "nodesel_dfs.h"
#include "branch_mostinf.h"
#include "branch_leastinf.h"


#if 0
static const int   nrows = 2;
static const int   nvars = 2;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MAXIMIZE;
static const char* var_name[] = { "var1"  , "var2"   };
static const Real  var_obj [] = {  1.0    ,  1.0     };     /* max  +x1 +x2 */
static const Real  var_lb  [] = {  1.0    ,  0.0     };     /*   1 <= x1 <= 10 integer */
static const Real  var_ub  [] = { 10.0    ,  5.0     };     /*   0 <= x2 <= 5  integer */

static const char* row_name[] = { "lin1"  , "lin2"   };     /* such that */
static const int   row_len [] = {        2,        2 };     /*   lin1: 0 <= +2x1 + x2 <= 4 */
static const int   row_idx [] = {   0,   1,   0,   1 };     /*   lin2: 0 <= + x1 +2x2 <= 3 */
static const Real  row_val [] = { 2.0, 1.0, 1.0, 2.0 };
static const Real  row_lhs [] = {      0.0,      0.0 };
static const Real  row_rhs [] = {      4.0,      3.0 };
#endif

#if 1
static const int   nrows = 3;
static const int   nvars = 2;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MAXIMIZE;
static const char* var_name[] = { "var1"  , "var2"   };
static const Real  var_obj [] = {  1.0    ,  1.0     };           /* max  +x1 +x2 */
static const Real  var_lb  [] = {  1.0    ,  0.0     };           /*   1 <= x1 <= 10 integer */
static const Real  var_ub  [] = { 10.0    ,  5.0     };           /*   0 <= x2 <= 5  integer */

static const char* row_name[] = { "lin1"  , "lin2"  , "lin3"   }; /* such that */
static const int   row_len [] = {        2,        2,        2 }; /*   lin1: 0 <= +2x1 + x2 <=  4 */
static const int   row_idx [] = {   0,   1,   0,   1,   0,   1 }; /*   lin2: 0 <= + x1 +2x2 <=  3 */
static const Real  row_val [] = { 2.0, 1.0, 1.0, 2.0, 1.0, 1.0 }; /*   lin3: 0 <= + x1 + x2 <= 20 */
static const Real  row_lhs [] = {      0.0,      0.0,      0.0 };
static const Real  row_rhs [] = {      4.0,      3.0,     20.0 };
#endif

#if 0
static const int   nrows = 2;
static const int   nvars = 2;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MINIMIZE;
static const char* var_name[] = { "var1"  , "var2"   };
static const Real  var_obj [] = { -2.0    , +1.0     };     /* min -2x1 +x2 */
static const Real  var_lb  [] = {  0.0    ,  0.0     };     /*   0 <= x1 <= 10 integer */
static const Real  var_ub  [] = { 10.0    , 10.0     };     /*   0 <= x2 <= 10 integer */

static const char* row_name[] = { "lin1"  , "lin2"   };     /* such that */
static const int   row_len [] = {        2,        1 };     /*   lin1: 1 <= -x1 +x2 */
static const int   row_idx [] = {   0,   1,        1 };     /*   lin2: 0 <=     -x2 */
static const Real  row_val [] = {-1.0, 1.0,     -1.0 };
static const Real  row_lhs [] = {      1.0,      0.0 };
static const Real  row_rhs [] = { 100000.0, 100000.0 };
#endif

#if 0
static const int   nrows = 3;
static const int   nvars = 3;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MINIMIZE;
static const char* var_name[] = { "var1"  , "var2",   "var3"  };       /* min -2x1 +x2 + 100000x3 */
static const Real  var_obj [] = { -2.0    ,  +1.0 , 100000.0  };       /*   0 <= x1 <=  10 integer */
static const Real  var_lb  [] = {  0.0    ,   0.0 ,      0.0  };       /*   0 <= x2 <=  10 integer */
static const Real  var_ub  [] = { 10.0    ,  10.0 ,    100.0  };       /*   0 <= x3 <= 100 integer */

static const char* row_name[] = { "lin1"       , "lin2"  , "lin3"   }; /* such that */
static const int   row_len [] = {             3,        1,        2 }; /*   lin1: 1 <= -x1 +x2 +x3       */
static const int   row_idx [] = {   0,   1,   2,        1,   0,   1 }; /*   lin2: 0 <=     -x2           */
static const Real  row_val [] = {-1.0, 1.0, 1.0,     -1.0, 1.0, 1.0 }; /*   lin3: 0 <=  x1 +x2     <= 30 */
static const Real  row_lhs [] = {           1.0,      0.0,      0.0 };
static const Real  row_rhs [] = {      100000.0, 100000.0,     30.0 };
#endif

int
main(int argc, char **argv)
{
   SCIP* scip = NULL;
   VAR** vars;
   VAR** rowvars;
   Real* rowvals;
   int r;
   int v;
   int i;
   int pos;

   SCIPprintVersion(NULL);

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
   CHECK_OKAY( SCIPincludeReaderMPS(scip) );
   CHECK_OKAY( SCIPincludeDispDefault(scip) );
   CHECK_SCIP( SCIPincludeConsHdlrIntegral(scip) );
   CHECK_SCIP( SCIPincludeConsHdlrLinear(scip) );
   CHECK_SCIP( SCIPincludeNodeselBfs(scip) );
   CHECK_SCIP( SCIPincludeNodeselDfs(scip) );
   CHECK_SCIP( SCIPincludeBranchruleMostinf(scip) );
   CHECK_SCIP( SCIPincludeBranchruleLeastinf(scip) );


   /********************
    * Problem Creation *
    ********************/

#if 0
   printf("\ncreate problem\n");

   /* create problem */
   CHECK_SCIP( SCIPcreateProb(scip, "test.lp") );
   CHECK_SCIP( SCIPsetObjsense(scip, objsen) );

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
      CHECK_SCIP( SCIPcreateConsLinear(scip, &cons, row_name[r], row_len[r], rowvars, rowvals,
                     row_lhs[r], row_rhs[r], TRUE, FALSE) );
      CHECK_SCIP( SCIPaddCons(scip, cons) ); /* add as a global constraint */
      CHECK_SCIP( SCIPreleaseCons(scip, &cons) );
   }

   /* release variables, because we don't need them any longer */
   for( v = 0; v < nvars; ++v )
   {
      CHECK_SCIP( SCIPreleaseVar(scip, &vars[v]) );
   }

#else
   if( argc < 2 )
   {
      printf("syntax: %s <problem>\n", argv[0]);
      return 0;
   }

   printf("\nread problem <%s>\n", argv[1]);
   CHECK_OKAY( SCIPreadProb(scip, argv[1]) );

   /*CHECK_OKAY( SCIPreadProb(scip, "IP/miplib/egout.mps") );*/
   /*CHECK_OKAY( SCIPreadProb(scip, "IP/miplib/markshare1.mps") );*/
   /*CHECK_OKAY( SCIPreadProb(scip, "IP/miplib/rentacar.mps") );*/
   /*CHECK_OKAY( SCIPreadProb(scip, "IP/miplib/dano3mip.mps") );*/

#endif


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   printf("\nsolve problem\n");
   CHECK_SCIP( SCIPsolve(scip) );

#if 0
   /* free solution process */
   printf("\nfree problem solution\n");
   CHECK_SCIP( SCIPfreeSolve(scip) );

   /* solve problem again */
   printf("\nsolve problem again\n");
   CHECK_SCIP( SCIPsolve(scip) );
#endif

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

