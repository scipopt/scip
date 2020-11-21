/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   activity.c
 * @brief  unittest for the activity datastructure in misc.c
 * @author Merlin Viernickel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/pub_misc.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_RESOURCEACTIVITY* activity;
static SCIP_VAR* var;
static int duration;
static int demand;

static
void setup(void)
{
   /* create scip */
   SCIP_CALL( SCIPcreate(&scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create variable */
   SCIP_CALL( SCIPcreateVarBasic(scip, &var, "var", -5, 5, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, var) );

   /* create resource activity */
   duration = 5;
   demand = 10;
   SCIP_CALL( SCIPactivityCreate(&activity, var, duration, demand) );
}


static
void teardown(void)
{
   /* free activity */
   SCIPactivityFree(&activity);

   /* free variable */
   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   /* free scip */
   SCIP_CALL( SCIPfree(&scip) );
}


TestSuite(activity, .init = setup, .fini = teardown);

Test(activity, setup_and_teardown, .description = "test that setup and teardown work correctly")
{
}

Test(activity, test_activity_getters, .description = "test that the resource activity returns entries correctly.")
{
   int energy;

   energy = duration * demand;

   cr_assert_eq(demand, SCIPactivityGetDemand(activity));
   cr_assert_eq(duration, SCIPactivityGetDuration(activity));
   cr_assert_eq(energy, SCIPactivityGetEnergy(activity));
   cr_assert_eq(var, SCIPactivityGetVar(activity));
}
