/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
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
#include "cons_setcover.h"
#include "cons_setpack.h"
#include "cons_setpart.h"
#include "nodesel_bfs.h"
#include "nodesel_dfs.h"
#include "branch_fullstrong.h"
#include "branch_mostinf.h"
#include "branch_leastinf.h"
#include "heur_diving.h"
#include "heur_rounding.h"
#include "sepa_gomory.h"



/*
 * Test Event Handler
 */

static
DECL_EVENTFREE(eventFreeTest)
{
   /*printf("free test event handler\n");*/

   return SCIP_OKAY;
}

static
DECL_EVENTINIT(eventInitTest)
{
   /*printf("init test event handler\n");*/

   CHECK_OKAY( SCIPcatchEvent(scip, 
                  /* SCIP_EVENTTYPE_NODEACTIVATED
                     | SCIP_EVENTTYPE_NODEFEASIBLE
                     | SCIP_EVENTTYPE_NODEINFEASIBLE
                     | SCIP_EVENTTYPE_NODEBRANCHED
                     | SCIP_EVENTTYPE_FIRSTLPSOLVED
                     | SCIP_EVENTTYPE_LPSOLVED
                     | SCIP_EVENTTYPE_POORSOLFOUND
                     | */
                  SCIP_EVENTTYPE_BESTSOLFOUND,
                  eventhdlr, NULL) );
   
   return SCIP_OKAY;
}

static
DECL_EVENTEXIT(eventExitTest)
{
   /*printf("exit test event handler\n");*/

   return SCIP_OKAY;
}

static
DECL_EVENTDELETE(eventDeleteTest)
{
   /*printf("delete test event handler\n");*/

   return SCIP_OKAY;
}

static
DECL_EVENTEXEC(eventExecTest)
{
   EVENTTYPE eventtype;

   eventtype = SCIPeventGetType(event);
   /*printf("exec test event handler: eventtype=0x%x\n", eventtype);*/

   /*???????????*/
   {
      char lpname[MAXSTRLEN];
      sprintf(lpname, "lp%lld.lp", SCIPgetNodenum(scip));
      CHECK_OKAY( SCIPwriteLP(scip, lpname) );
   }

   return SCIP_OKAY;
}

static
RETCODE includeTestEventHdlr(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CHECK_OKAY( SCIPincludeEventhdlr(scip, "testeventhdlr", "test event handler description",
                  eventFreeTest, eventInitTest, eventExitTest, eventDeleteTest, eventExecTest, NULL) );

   return SCIP_OKAY;
}



#if 1
static const int   nrows = 2;
static const int   nvars = 2;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MAXIMIZE;
static const char* var_name[] = { "var1"  , "var2"   };
static const Real  var_obj [] = {  1.0    ,  1.0     };     /* max  +x1 +x2 */
static const Real  var_lb  [] = {  0.0    ,  0.0     };     /*   0 <= x1 <= 10 integer */
static const Real  var_ub  [] = { 10.0    ,  5.0     };     /*   0 <= x2 <= 5  integer */

static const char* row_name[] = { "lin1"  , "lin2"   };     /* such that */
static const int   row_len [] = {        2,        2 };     /*   lin1:      +2x1 + x2 <= 2 */
static const int   row_idx [] = {   0,   1,   0,   1 };     /*   lin2:      + x1 +2x2 <= 2 */
static const Real  row_val [] = { 2.0, 1.0, 1.0, 2.0 };
static const Real  row_lhs [] = {   -1e+20,   -1e+20 };
static const Real  row_rhs [] = {      2.0,      2.0 };
#endif

#if 0
static const int   nrows = 2;
static const int   nvars = 2;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MAXIMIZE;
static const char* var_name[] = { "var1"  , "var2"   };
static const Real  var_obj [] = {  1.0    ,  1.0     };     /* max  +x1 +x2 */
static const Real  var_lb  [] = {  0.0    ,  0.0     };     /*   0 <= x1 <= 10 integer */
static const Real  var_ub  [] = { 10.0    ,  5.0     };     /*   0 <= x2 <= 5  integer */

static const char* row_name[] = { "lin1"  , "lin2"   };     /* such that */
static const int   row_len [] = {        2,        2 };     /*   lin1:       +2x1 + x2 <= 2 */
static const int   row_idx [] = {   0,   1,   0,   1 };     /*   lin2: -2 <= - x1 -2x2      */
static const Real  row_val [] = { 2.0, 1.0,-1.0,-2.0 };
static const Real  row_lhs [] = {   -1e+20,     -2.0 };
static const Real  row_rhs [] = {      2.0,    1e+20 };
#endif

#if 0
static const int   nrows = 2;
static const int   nvars = 2;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MAXIMIZE;
static const char* var_name[] = { "var1"  , "var2"   };
static const Real  var_obj [] = { -1.0    ,  1.0     };
static const Real  var_lb  [] = { -9.0    ,  0.0     };
static const Real  var_ub  [] = {  1.0    ,  5.0     };

static const char* row_name[] = { "lin1"  , "lin2"   };
static const int   row_len [] = {        2,        2 };
static const int   row_idx [] = {   0,   1,   0,   1 };
static const Real  row_val [] = {-2.0,+1.0, 1.0,-2.0 };
static const Real  row_lhs [] = {   -1e+20,     -1.0 };
static const Real  row_rhs [] = {      0.0,    1e+20 };
#endif

#if 0
static const int   nrows = 2;
static const int   nvars = 2;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MAXIMIZE;
static const char* var_name[] = { "var1"  , "var2"   };
static const Real  var_obj [] = {  1.0    ,  1.0     };     /* max  +x1 +x2 */
static const Real  var_lb  [] = {  1.0    ,  0.0     };     /*   1 <= x1 <= 10 integer */
static const Real  var_ub  [] = { 10.0    ,  5.0     };     /*   0 <= x2 <= 5  integer */

static const char* row_name[] = { "lin1"  , "lin2"   };     /* such that */
static const int   row_len [] = {        2,        2 };     /*   lin1:      +2x1 + x2 <= 4 */
static const int   row_idx [] = {   0,   1,   0,   1 };     /*   lin2:      + x1 +2x2 <= 3 */
static const Real  row_val [] = { 2.0, 1.0, 1.0, 2.0 };
static const Real  row_lhs [] = {   -1e+20,   -1e+20 };
static const Real  row_rhs [] = {      4.0,      3.0 };
#endif

#if 0
static const int   nrows = 4;
static const int   nvars = 4;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MAXIMIZE;
static const char* var_name[] = { "var1"  , "var2"  , "var3"  , "var4"  };
static const Real  var_obj [] = {  1.0    ,  1.0    ,  1.0    ,  1.0    };     /* max  +x1 +x2 +x3 +x4 */
static const Real  var_lb  [] = {  1.0    ,  0.0    ,  1.0    ,  0.0    };     /*   1 <= x1,x3 <= 10 integer */
static const Real  var_ub  [] = { 10.0    ,  5.0    , 10.0    ,  5.0    };     /*   0 <= x2,x4 <= 5  integer */

static const char* row_name[] = { "lin1"  , "lin2"  , "lin3"  , "lin4"  };     /* such that */
static const int   row_len [] = {        2,        2,        2,        2};     /*   lin1:      +2x1 + x2 <= 4 */
static const int   row_idx [] = {   0,   1,   0,   1,   2,   3,   2,   3};     /*   lin2:      + x1 +2x2 <= 3 */
static const Real  row_val [] = { 2.0, 1.0, 1.0, 2.0,-2.0,-1.0,-1.0,-2.0};     /*   lin3:-4 <= -2x1 - x2      */
static const Real  row_lhs [] = {   -1e+20,   -1e+20,     -4.0,     -3.0};     /*   lin4:-3 <= - x1 -2x2      */
static const Real  row_rhs [] = {      4.0,      3.0,    1e+20,    1e+20};
#endif

#if 0
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

#if 0
static const int   nrows = 3;
static const int   nvars = 3;

static const OBJSENSE objsen  = SCIP_OBJSENSE_MINIMIZE;
static const char* var_name[] = { "var1"  , "var2",   "var3"  };       /* min -2x1 +x2 + 3x3 */
static const Real  var_obj [] = { -2.0    ,  +1.0 ,     +3.0  };       /*   0 <= x1 <=   1 integer */
static const Real  var_lb  [] = {  0.0    ,   0.0 ,      0.0  };       /*   0 <= x2 <=   1 integer */
static const Real  var_ub  [] = {  1.0    ,   1.0 ,      1.0  };       /*   0 <= x3 <=   1 integer */

static const char* row_name[] = { "lin1"       , "lin2"  , "lin3"   }; /* such that */
static const int   row_len [] = {             3,        1,        2 }; /*   lin1: 1 <=  x1 +x2 +x3       */
static const int   row_idx [] = {   0,   1,   2,        1,   0,   1 }; /*   lin2: 0 <=     -x2           */
static const Real  row_val [] = { 1.0, 1.0, 1.0,     -1.0, 1.0, 1.0 }; /*   lin3: 0 <=  x1 +x2     <= 30 */
static const Real  row_lhs [] = {           1.0,      0.0,      0.0 };
static const Real  row_rhs [] = {        1e+100,   1e+100,     30.0 };
#endif


static
RETCODE runSCIP(
   int              argc,
   char**           argv
   )
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

   ALLOC_OKAY( allocMemoryArray(&vars, nvars) );
   ALLOC_OKAY( allocMemoryArray(&rowvars, nvars) );
   ALLOC_OKAY( allocMemoryArray(&rowvals, nvars) );


   /*********
    * Setup *
    *********/

   printf("\nsetup SCIP\n");
   printf("==========\n\n");

   /* initialize SCIP */
   CHECK_OKAY( SCIPcreate(&scip) );

   /* include user defined callbacks */
   CHECK_OKAY( SCIPincludeReaderMPS(scip) );
   CHECK_OKAY( SCIPincludeDispDefault(scip) );
   CHECK_OKAY( SCIPincludeConsHdlrIntegral(scip) );
   CHECK_OKAY( SCIPincludeConsHdlrLinear(scip) );
   CHECK_OKAY( SCIPincludeConsHdlrSetcover(scip) );
   CHECK_OKAY( SCIPincludeConsHdlrSetpack(scip) );
   CHECK_OKAY( SCIPincludeConsHdlrSetpart(scip) );
   CHECK_OKAY( SCIPincludeNodeselBfs(scip) );
   CHECK_OKAY( SCIPincludeNodeselDfs(scip) );
   /*CHECK_OKAY( SCIPincludeBranchruleFullstrong(scip) );*/
   CHECK_OKAY( SCIPincludeBranchruleMostinf(scip) );
   CHECK_OKAY( SCIPincludeBranchruleLeastinf(scip) );
   CHECK_OKAY( SCIPincludeHeurDiving(scip) );
   CHECK_OKAY( SCIPincludeHeurRounding(scip) );
   CHECK_OKAY( SCIPincludeSepaGomory(scip) );
   
   /*CHECK_OKAY( includeTestEventHdlr(scip) );*/


   /**************
    * Parameters *
    **************/

   if( SCIPfileExists("scip.set") )
   {
      printf("reading parameter file <scip.set>\n");
      CHECK_OKAY( SCIPreadParams(scip, "scip.set") );
   }
   else
      printf("parameter file <scip.set> not found - using default parameters\n");



   /********************
    * Problem Creation *
    ********************/

#if 0
   printf("\ncreate problem\n");
   printf("==============\n\n");

   /* create problem */
   CHECK_OKAY( SCIPcreateProb(scip, "test.lp", NULL, NULL, NULL) );
   CHECK_OKAY( SCIPsetObjsense(scip, objsen) );

   /* create necessary variables */
   for( v = 0; v < nvars; ++v )
   {
      CHECK_OKAY( SCIPcreateVar(scip, &vars[v], var_name[v], var_lb[v], var_ub[v], var_obj[v], 
                     SCIP_VARTYPE_INTEGER, FALSE) );
      CHECK_OKAY( SCIPaddVar(scip, vars[v]) );
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
      CHECK_OKAY( SCIPcreateConsLinear(scip, &cons, row_name[r], row_len[r], rowvars, rowvals,
                     row_lhs[r], row_rhs[r], TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
      CHECK_OKAY( SCIPaddCons(scip, cons) ); /* add as a global constraint */
      CHECK_OKAY( SCIPreleaseCons(scip, &cons) );
   }

   /* release variables, because we don't need them any longer */
   for( v = 0; v < nvars; ++v )
   {
      CHECK_OKAY( SCIPreleaseVar(scip, &vars[v]) );
   }

#else
   if( argc < 2 )
   {
      printf("syntax: %s <problem>\n", argv[0]);
      return SCIP_OKAY;
   }

   printf("\nread problem <%s>\n", argv[1]);
   printf("============\n\n");
   CHECK_OKAY( SCIPreadProb(scip, argv[1]) );
#endif


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   printf("\nsolve problem\n");
   printf("=============\n\n");
   CHECK_OKAY( SCIPsolve(scip) );

#if 0
   /* free solution process */
   printf("\nfree problem solution\n");
   printf("=====================\n\n");
   CHECK_OKAY( SCIPfreeSolve(scip) );

   /* solve problem again */
   printf("\nsolve problem again\n");
   printf("===================\n\n");
   CHECK_OKAY( SCIPsolve(scip) );
#endif

#if 1
   printf("\nprimal solution:\n");
   printf("================\n\n");
   CHECK_OKAY( SCIPprintBestSol(scip, NULL) );
#endif

#ifndef NDEBUG
   /*SCIPdebugMemory(scip);*/
#endif


   /**************
    * Statistics *
    **************/

   printf("\nStatistics\n");
   printf("==========\n\n");

   CHECK_OKAY( SCIPprintStatistics(scip, NULL) );


   /********************
    * Deinitialization *
    ********************/

   printf("\nfree SCIP\n");
   printf("=========\n\n");

   /* free SCIP */
   CHECK_OKAY( SCIPfree(&scip) );


   /*****************************
    * Local Memory Deallocation *
    *****************************/

   freeMemoryArray(&rowvals);
   freeMemoryArray(&rowvars);
   freeMemoryArray(&vars);

#ifndef NDEBUG
   memoryCheckEmpty();
#endif

   return SCIP_OKAY;
}

int
main(
   int              argc,
   char**           argv
   )
{
   RETCODE retcode;

   todoMessage("implement remaining events");
   todoMessage("avoid addition of identical rows");
   todoMessage("avoid addition of identical constraints");
   todoMessage("pricing for pseudo solutions");
   todoMessage("integrality check on objective function, abort if gap is below 1.0");
   todoMessage("numerical problems in tree->actpseudoobjval if variable's bounds are infinity");
   todoMessage("implement reduced cost fixing");
   todoMessage("statistics: count domain reductions and constraint additions of constraint handlers");
   todoMessage("it's a bit ugly, that user call backs may be called before the nodequeue was processed");
   todoMessage("information method if parameter changed");

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode, stderr);
      return -1;
   }

   return 0;
}
