/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dual.c
 * @brief  unit tests for checking dual events
 * @author João Dionísio
 * @author Ambros Gleixner
 * @author Dominik Kamp
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/scip_event.h"

#define EVENTHDLR_NAME         "dualboundimproved"
#define EVENTHDLR_DESC         "event handler for catching dual bound improvements"


/** EVENT HANDLER **/

/** structure to store event handler specific information */
struct SCIP_EventhdlrData
{
   SCIP_Bool ispseudoreached;                /**< whether pseudo objective is reached */
   SCIP_Bool islpreached;                    /**< whether LP value is reached */
   SCIP_Bool isoptreached;                   /**< whether optimal bound is reached */
   SCIP_Real lastdualbound;                  /**< last dual bound encountered before */
#ifdef SCIP_WITH_EXACTSOLVE
   SCIP_RATIONAL* lastdualboundexact;        /**< exact last dual bound encountered before */
#endif
};

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitDualBoundImproved)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

#ifdef SCIP_WITH_EXACTSOLVE
   if( SCIPisExact(scip) )
   {
      /* create last bound */
      SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &eventhdlrdata->lastdualboundexact) );

      /* initialize last bound */
      SCIPrationalSetNegInfinity(eventhdlrdata->lastdualboundexact);
   }
   else
#endif
   {
      /* initialize last bound */
      eventhdlrdata->lastdualbound = -SCIPinfinity(scip);
   }

   /* notify SCIP that your event handler wants to react on the event type best solution found */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_DUALBOUNDIMPROVED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitDualBoundImproved)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);

   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

#ifdef SCIP_WITH_EXACTSOLVE
   if( SCIPisExact(scip) )
   {
      SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
      assert(eventhdlrdata != NULL);

      /* free last bound */
      SCIPrationalFreeBlock(SCIPblkmem(scip), &eventhdlrdata->lastdualboundexact);
   }
#endif

   /* notify SCIP that your event handler wants to drop the event type best solution found */
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_DUALBOUNDIMPROVED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecDualBoundImproved)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(eventhdlr != NULL);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_DUALBOUNDIMPROVED);

   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

   SCIPdebugMsg(scip, "exec method of event handler for dual bound improvement\n");

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

#ifdef SCIP_WITH_EXACTSOLVE
   if( SCIPisExact(scip) )
   {
      SCIP_RATIONAL* dualboundexact;

      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &dualboundexact) );

      /* get dual bound */
      SCIPgetDualboundExact(scip, dualboundexact);

      /* print dual bound */
      SCIPinfoMessage(scip, NULL, "found new dual bound with value ");
      SCIPrationalMessage(SCIPgetMessagehdlr(scip), NULL, dualboundexact);
      SCIPinfoMessage(scip, NULL, " in problem <%s>\n", SCIPgetProbName(scip));

      /* ensure finite improvement */
      assert(!SCIPrationalIsAbsInfinity(dualboundexact));
      assert(SCIPrationalIsGT(dualboundexact, eventhdlrdata->lastdualboundexact));

      /* detect bound state */
      if( !eventhdlrdata->ispseudoreached )
      {
         if( SCIPrationalIsZero(dualboundexact) )
            eventhdlrdata->ispseudoreached = TRUE;
      }
      else if( !eventhdlrdata->islpreached )
      {
         if( SCIPrationalIsPositive(dualboundexact) && SCIPrationalIsLEReal(dualboundexact, 1.0) )
            eventhdlrdata->islpreached = TRUE;
      }
      else if( !eventhdlrdata->isoptreached )
      {
         if( SCIPrationalIsGTReal(dualboundexact, 1.0) )
            eventhdlrdata->isoptreached = TRUE;
      }

      /* save last bound */
      SCIPrationalSetRational(eventhdlrdata->lastdualboundexact, dualboundexact);

      SCIPrationalFreeBuffer(SCIPbuffer(scip), &dualboundexact);
   }
   else
#endif
   {
      SCIP_Real dualbound;

      /* get dual bound */
      dualbound = SCIPgetDualbound(scip);

      /* print dual bound */
      SCIPinfoMessage(scip, NULL, "found new dual bound with value %.15g in problem <%s>\n", dualbound,
         SCIPgetProbName(scip));

      /* ensure finite improvement */
      assert(!SCIPisInfinity(scip, ABS(dualbound)));
      assert(dualbound > eventhdlrdata->lastdualbound);

      /* detect bound state */
      if( !eventhdlrdata->ispseudoreached )
      {
         if( SCIPisZero(scip, dualbound) )
            eventhdlrdata->ispseudoreached = TRUE;
      }
      else if( !eventhdlrdata->islpreached )
      {
         if( SCIPisPositive(scip, dualbound) && SCIPisLE(scip, dualbound, 1.0) )
            eventhdlrdata->islpreached = TRUE;
      }
      else if( !eventhdlrdata->isoptreached )
      {
         if( SCIPisGT(scip, dualbound, 1.0) )
            eventhdlrdata->isoptreached = TRUE;
      }

      /* save last bound */
      eventhdlrdata->lastdualbound = dualbound;
   }

   return SCIP_OKAY;
}

/** includes event handler for best solution found */
static
SCIP_RETCODE SCIPincludeEventHdlrDualBoundImproved(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   assert(scip != NULL);
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata) );

   /* create event handler for dual bound improved */
   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecDualBoundImproved, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitDualBoundImproved) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitDualBoundImproved) );

   return SCIP_OKAY;
}

/** GLOBAL VARIABLES **/
static SCIP* scip_test = NULL;

/** TEST SUITES **/
static
void setup(void)
{
   cr_assert(scip_test == NULL);

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip_test) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip_test) );
   SCIP_CALL( SCIPincludeEventHdlrDualBoundImproved(scip_test) );
}

static
void teardown(void)
{
   SCIP_EVENTHDLR* eventhdlr = SCIPfindEventhdlr(scip_test, EVENTHDLR_NAME);
   cr_assert(eventhdlr != NULL);
   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   cr_assert(eventhdlrdata != NULL);

   /* free event handler data */
   SCIPfreeBlockMemory(scip_test, &eventhdlrdata);

   SCIP_CALL( SCIPfree(&scip_test) );

   cr_assert(scip_test == NULL);
   cr_assert(BMSgetMemoryUsed() == 0, "There is a memory leak!");
}

TestSuite(events, .init = setup, .fini = teardown);

/* TESTS */

Test(events, dualboundimprovedreal, .description = "tests real SCIP_EVENTTYPE_DUALBOUNDIMPROVED")
{
   /* create the following problem with LP dual bound 1 and optimal value 1 + mu with mu = 1e-5, such that the event
    * for the lowerbound improvement from 1 to 1 + mu should happen in real mode:
    *
    * min   x0
    * s.t.  x0 >= 1 - mu * x1 + mu * (1 - x1)  equivalent  x0 + 2mu * x1 >= 1 + mu
    *       x0 >= 1 + mu * x1 - mu * (1 - x1)  equivalent  x0 - 2mu * x1 >= 1 - mu
    *       x0 >= 0, x1 binary
    */
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIP_CALL( SCIPsetPresolving(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetSeparating(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIPinfoMessage(scip_test, NULL, "create test problem ...\n");

   SCIP_CALL( SCIPcreateProbBasic(scip_test, "real") );

   SCIPinfoMessage(scip_test, NULL, "... problem created\n");

   SCIP_CALL( SCIPcreateVarBasic(scip_test, &vars[0], "x0", 0.0, SCIPinfinity(scip_test), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[0]) );

   SCIP_CALL( SCIPcreateVarBasic(scip_test, &vars[1], "x1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[1]) );

   SCIPinfoMessage(scip_test, NULL, "... variables created\n");

   vals[0] = 1.0;
   rhs = SCIPinfinity(scip_test);

   vals[1] = 2e-5;
   lhs = 1.0+1e-5;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip_test, &cons, "c0", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   vals[1] = -2e-5;
   lhs = 1.0-1e-5;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip_test, &cons, "c1", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   SCIPinfoMessage(scip_test, NULL, "... constraint created\n");

   eventhdlr = SCIPfindEventhdlr(scip_test, EVENTHDLR_NAME);
   cr_assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   cr_assert(eventhdlrdata != NULL);

   eventhdlrdata->ispseudoreached = FALSE;
   eventhdlrdata->islpreached = FALSE;
   eventhdlrdata->isoptreached = FALSE;

   SCIP_CALL( SCIPsolve(scip_test) );

   /* expect correct optimum */
   lhs = 1.0+1e-5;
   rhs = SCIPgetSolVal(scip_test, SCIPgetBestSol(scip_test), vars[0]);
   cr_expect(SCIPisEQ(scip_test, rhs, lhs), "Solution value wrong");
   rhs = SCIPgetDualbound(scip_test);
   cr_expect(SCIPisEQ(scip_test, rhs, lhs), "Dual bound wrong");
   rhs = SCIPgetPrimalbound(scip_test);
   cr_expect(SCIPisEQ(scip_test, rhs, lhs), "Primal bound wrong");

   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[1]) );

   /* expect at least three improvements: 0 as pseudo bound, 1 as LP bound, and 1 + mu as optimal value */
   cr_assert(eventhdlrdata->ispseudoreached, "Pseudo objective missed");
   cr_assert(eventhdlrdata->islpreached, "LP value missed");
   cr_assert(eventhdlrdata->isoptreached, "Optimal bound missed");
}

Test(events, dualboundimprovedrealimplint, .description = "tests realimplint SCIP_EVENTTYPE_DUALBOUNDIMPROVED")
{
   /* create the following problem with LP dual bound 0.5 and implied integral optimal value 2, such that an event
    * for a lowerbound improvement from 1 to 2 should occur:
    *
    * min   x0
    * s.t.  x0 >= 0.5 - 1.5 * x1 + 1.5 * (1 - x1)  equivalent  x0 + 3 * x1 >= 2
    *       x0 >= 0.5 + 1.5 * x1 - 1.5 * (1 - x1)  equivalent  x0 - 3 * x1 >= -1
    *       x0 >= 0 contimplint, x1 binary
    */
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_Real lhs;
   SCIP_Real rhs;

   SCIP_CALL( SCIPsetPresolving(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetSeparating(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIPinfoMessage(scip_test, NULL, "create test problem ...\n");

   SCIP_CALL( SCIPcreateProbBasic(scip_test, "realimplint") );

   SCIPinfoMessage(scip_test, NULL, "... problem created\n");

   SCIP_CALL( SCIPcreateVarImpl(scip_test, &vars[0], "x0", 0.0, SCIPinfinity(scip_test), 1.0,
         SCIP_VARTYPE_CONTINUOUS, SCIP_IMPLINTTYPE_WEAK, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[0]) );

   SCIP_CALL( SCIPcreateVarBasic(scip_test, &vars[1], "x1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[1]) );

   SCIPinfoMessage(scip_test, NULL, "... variables created\n");

   vals[0] = 1.0;
   rhs = SCIPinfinity(scip_test);

   vals[1] = 3.0;
   lhs = 2.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip_test, &cons, "c0", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   vals[1] = -3.0;
   lhs = -1.0;
   SCIP_CALL( SCIPcreateConsBasicLinear(scip_test, &cons, "c1", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   SCIPinfoMessage(scip_test, NULL, "... constraint created\n");

   eventhdlr = SCIPfindEventhdlr(scip_test, EVENTHDLR_NAME);
   cr_assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   cr_assert(eventhdlrdata != NULL);

   eventhdlrdata->ispseudoreached = FALSE;
   eventhdlrdata->islpreached = FALSE;
   eventhdlrdata->isoptreached = FALSE;

   SCIP_CALL( SCIPsolve(scip_test) );

   /* expect correct optimum */
   lhs = 2.0;
   rhs = SCIPgetSolVal(scip_test, SCIPgetBestSol(scip_test), vars[0]);
   cr_expect(SCIPisEQ(scip_test, rhs, lhs), "Solution value wrong");
   rhs = SCIPgetDualbound(scip_test);
   cr_expect(SCIPisEQ(scip_test, rhs, lhs), "Dual bound wrong");
   rhs = SCIPgetPrimalbound(scip_test);
   cr_expect(SCIPisEQ(scip_test, rhs, lhs), "Primal bound wrong");

   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[1]) );

   /* expect at least three improvements: 0 as pseudo bound, 1 as rounded LP bound, and 2 as optimal value */
   cr_assert(eventhdlrdata->ispseudoreached, "Pseudo objective missed");
   cr_assert(eventhdlrdata->islpreached, "LP value missed");
   cr_assert(eventhdlrdata->isoptreached, "Optimal bound missed");
}

#ifdef SCIP_WITH_EXACTSOLVE
Test(events, dualboundimprovedexact, .description = "tests exact SCIP_EVENTTYPE_DUALBOUNDIMPROVED")
{
   /* create the following problem with LP dual bound 1 and optimal value 1 + mu with mu = 1e-16, such that the event
    * for the lowerbound improvement from 1 to 1 + mu can only occur in exact mode:
    *
    * min   x0
    * s.t.  x0 >= 1 - mu * x1 + mu * (1 - x1)  equivalent  x0 + 2mu * x1 >= 1 + mu
    *       x0 >= 1 + mu * x1 - mu * (1 - x1)  equivalent  x0 - 2mu * x1 >= 1 - mu
    *       x0 >= 0, x1 binary
    */
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_RATIONAL* vals[2];
   SCIP_RATIONAL* lhs;
   SCIP_RATIONAL* rhs;

   SCIP_CALL( SCIPenableExactSolving(scip_test, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetSeparating(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIPinfoMessage(scip_test, NULL, "create test problem ...\n");

   SCIP_CALL( SCIPcreateProbBasic(scip_test, "exact") );

   SCIPinfoMessage(scip_test, NULL, "... problem created\n");

   SCIP_CALL( SCIPcreateVarBasic(scip_test, &vars[0], "x0", 0.0, SCIPinfinity(scip_test), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVarExactData(scip_test, vars[0], NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[0]) );

   SCIP_CALL( SCIPcreateVarBasic(scip_test, &vars[1], "x1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVarExactData(scip_test, vars[1], NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[1]) );

   SCIPinfoMessage(scip_test, NULL, "... variables created\n");

   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &vals[0]) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &vals[1]) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &lhs) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &rhs) );

   SCIPrationalSetString(vals[0], "1");
   SCIPrationalSetInfinity(rhs);

   SCIPrationalSetString(vals[1], "1/5000000000000000");
   SCIPrationalSetString(lhs, "10000000000000001/10000000000000000");
   SCIP_CALL( SCIPcreateConsBasicExactLinear(scip_test, &cons, "c0", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   SCIPrationalSetString(vals[1], "-1/5000000000000000");
   SCIPrationalSetString(lhs, "9999999999999999/10000000000000000");
   SCIP_CALL( SCIPcreateConsBasicExactLinear(scip_test, &cons, "c1", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   SCIPinfoMessage(scip_test, NULL, "... constraint created\n");

   eventhdlr = SCIPfindEventhdlr(scip_test, EVENTHDLR_NAME);
   cr_assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   cr_assert(eventhdlrdata != NULL);

   eventhdlrdata->ispseudoreached = FALSE;
   eventhdlrdata->islpreached = FALSE;
   eventhdlrdata->isoptreached = FALSE;

   SCIP_CALL( SCIPsolve(scip_test) );

   /* expect correct optimum */
   SCIPrationalSetString(lhs, "10000000000000001/10000000000000000");
   SCIPgetSolValExact(scip_test, SCIPgetBestSol(scip_test), vars[0], rhs);
   cr_expect(SCIPrationalIsEQ(rhs, lhs), "Solution value wrong");
   SCIPgetDualboundExact(scip_test, rhs);
   cr_expect(SCIPrationalIsEQ(rhs, lhs), "Dual bound wrong");
   SCIPgetPrimalboundExact(scip_test, rhs);
   cr_expect(SCIPrationalIsEQ(rhs, lhs), "Primal bound wrong");

   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &rhs);
   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &lhs);
   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &vals[1]);
   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &vals[0]);

   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[1]) );

   /* expect at least three improvements: 0 as pseudo bound, 1 as LP bound, and 1 + mu as optimal value */
   cr_assert(eventhdlrdata->ispseudoreached, "Pseudo objective missed");
   cr_assert(eventhdlrdata->islpreached, "LP value missed");
   cr_assert(eventhdlrdata->isoptreached, "Optimal bound missed");
}

Test(events, dualboundimprovedexactimplint, .description = "tests exactimplint SCIP_EVENTTYPE_DUALBOUNDIMPROVED")
{
   /* create the following problem with LP dual bound 1/2 and implied integral optimal value 2, such that an event
    * for a lowerbound improvement from 1 to 2 should occur:
    *
    * min   x0
    * s.t.  x0 >= 1/2 - 3/2 * x1 + 3/2 * (1 - x1)  equivalent  x0 + 3 * x1 >= 2
    *       x0 >= 1/2 + 3/2 * x1 - 3/2 * (1 - x1)  equivalent  x0 - 3 * x1 >= -1
    *       x0 >= 0 contimplint, x1 binary
    */
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_CONS* cons;
   SCIP_VAR* vars[2];
   SCIP_RATIONAL* vals[2];
   SCIP_RATIONAL* lhs;
   SCIP_RATIONAL* rhs;

   SCIP_CALL( SCIPenableExactSolving(scip_test, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetSeparating(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip_test, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIPinfoMessage(scip_test, NULL, "create test problem ...\n");

   SCIP_CALL( SCIPcreateProbBasic(scip_test, "exactimplint") );

   SCIPinfoMessage(scip_test, NULL, "... problem created\n");

   SCIP_CALL( SCIPcreateVarImpl(scip_test, &vars[0], "x0", 0.0, SCIPinfinity(scip_test), 1.0,
         SCIP_VARTYPE_CONTINUOUS, SCIP_IMPLINTTYPE_WEAK, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVarExactData(scip_test, vars[0], NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[0]) );

   SCIP_CALL( SCIPcreateVarBasic(scip_test, &vars[1], "x1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVarExactData(scip_test, vars[1], NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip_test, vars[1]) );

   SCIPinfoMessage(scip_test, NULL, "... variables created\n");

   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &vals[0]) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &vals[1]) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &lhs) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip_test), &rhs) );

   SCIPrationalSetString(vals[0], "1");
   SCIPrationalSetInfinity(rhs);

   SCIPrationalSetString(vals[1], "3");
   SCIPrationalSetString(lhs, "2");
   SCIP_CALL( SCIPcreateConsBasicExactLinear(scip_test, &cons, "c0", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   SCIPrationalSetString(vals[1], "-3");
   SCIPrationalSetString(lhs, "-1");
   SCIP_CALL( SCIPcreateConsBasicExactLinear(scip_test, &cons, "c1", 2, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip_test, cons) );
   SCIP_CALL( SCIPreleaseCons(scip_test, &cons) );

   SCIPinfoMessage(scip_test, NULL, "... constraint created\n");

   eventhdlr = SCIPfindEventhdlr(scip_test, EVENTHDLR_NAME);
   cr_assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   cr_assert(eventhdlrdata != NULL);

   eventhdlrdata->ispseudoreached = FALSE;
   eventhdlrdata->islpreached = FALSE;
   eventhdlrdata->isoptreached = FALSE;

   SCIP_CALL( SCIPsolve(scip_test) );

   /* expect correct optimum */
   SCIPrationalSetString(lhs, "2");
   SCIPgetSolValExact(scip_test, SCIPgetBestSol(scip_test), vars[0], rhs);
   cr_expect(SCIPrationalIsEQ(rhs, lhs), "Solution value wrong");
   SCIPgetDualboundExact(scip_test, rhs);
   cr_expect(SCIPrationalIsEQ(rhs, lhs), "Dual bound wrong");
   SCIPgetPrimalboundExact(scip_test, rhs);
   cr_expect(SCIPrationalIsEQ(rhs, lhs), "Primal bound wrong");

   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &rhs);
   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &lhs);
   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &vals[1]);
   SCIPrationalFreeBlock(SCIPblkmem(scip_test), &vals[0]);

   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[0]) );
   SCIP_CALL( SCIPreleaseVar(scip_test, &vars[1]) );

   /* expect at least three improvements: 0 as pseudo bound, 1 as rounded LP bound, and 2 as optimal value */
   cr_assert(eventhdlrdata->ispseudoreached, "Pseudo objective missed");
   cr_assert(eventhdlrdata->islpreached, "LP value missed");
   cr_assert(eventhdlrdata->isoptreached, "Optimal bound missed");
}
#endif
