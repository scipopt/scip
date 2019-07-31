/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disable_upgrades.c
 * @brief  unit tests for printing and disable_upgrades linear constraints
 * @author Gregor Hendel
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_or.h"
#include "scip/pub_cons.h"
#include <stdio.h>

#define MEMSIZE 10
#define FNAME ".disable_upgrades-or-%s.cip"


/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;

/* TEST SUITE */

/** the setup creates the necessary source and target SCIPs, initializes linear data */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
}

/** free all allocated memory */
static
void teardown(void)
{
   SCIPfree(&scip);

/*
    remove the cip-file for this test
   if( strncmp(filename, FNAME, 5)  == 0 )
      remove(filename);
*/

}


TestSuite(disable_upgrades, .init = setup, .fini = teardown);

/* TESTS  */
Test(disable_upgrades, create_and_free)
{
   /* calls setup and teardown */
}

Test(disable_upgrades, disable_upgrades_or, .description="disable upgrades of or-constraints to and-constraints to keep or-constraints during solution process")
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   int i;

   SCIP_CALL( SCIPreadProb(scip, "../check/instances/Or/or_constraint.cip", "cip") );

   conshdlr = SCIPfindConshdlr(scip, "or");

   nconss = SCIPgetNOrigConss(scip);
   conss = SCIPgetOrigConss(scip);

   cr_assert_eq(nconss, 8);

   /* todo set or constraints to be modifiable */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsGetHdlr(conss[i]) == conshdlr )
         SCIPsetConsModifiable(scip, conss[i], TRUE);
   }


   SCIP_CALL( SCIPsolve(scip) );
}

Test(disable_upgrades, write_problem, .description="test that CIP write method works for or constraints")
{
   SCIP_CALL( SCIPreadProb(scip, "../check/instances/Or/or_constraint.cip", "cip") );

   /* author bzfhende
    *
    * TODO @Helena
    */

}


