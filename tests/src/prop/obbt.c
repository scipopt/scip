/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   obbt.c
 * @brief  unit test for prop_obbt methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/prop_obbt.c"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* origx;
static SCIP_VAR* origy;
static SCIP_VAR* origz;

/* creates scip, problem, and includes quadratic and nonlinear constraint handlers */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include quadratic and nonlinear conshdlr */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );
   SCIP_CALL( SCIPincludeConshdlrQuadratic(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &origx, "x", -10.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &origy, "y", -5.0, 5.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &origz, "z", 1.0, 3.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, origx) );
   SCIP_CALL( SCIPaddVar(scip, origy) );
   SCIP_CALL( SCIPaddVar(scip, origz) );

   TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE);
   x = SCIPvarGetTransVar(origx);
   y = SCIPvarGetTransVar(origy);
   z = SCIPvarGetTransVar(origz);
}

/* frees scip */
static
void teardown(void)
{
   /* free variables */
   SCIP_CALL( SCIPreleaseVar(scip, &origz) );
   SCIP_CALL( SCIPreleaseVar(scip, &origy) );
   SCIP_CALL( SCIPreleaseVar(scip, &origx) );

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

/* helper function to add a row containing two variables via SCIPaddRowProbing() */
static
SCIP_RETCODE addRowProbing(SCIP* _scip, SCIP_VAR* _x, SCIP_VAR* _y, SCIP_Real coefx, SCIP_Real coefy, SCIP_Real lhs, SCIP_Real rhs)
{
   SCIP_ROW* row;

   SCIP_CALL( SCIPcreateEmptyRowUnspec(_scip, &row, "row", lhs, rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarToRow(_scip, row, _x, coefx) );
   SCIP_CALL( SCIPaddVarToRow(_scip, row, _y, coefy) );
   SCIP_CALL( SCIPaddRowProbing(_scip, row) );
   SCIP_CALL( SCIPreleaseRow(_scip, &row) );

   return SCIP_OKAY;
}

Test(projection, onerow, .init = setup, .fini=teardown,
   .description = "single row example"
   )
{
   SCIP_Bool infeasible;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real constant;
   SCIP_Bool lperror;

   /* construct LP */
   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );
   assert(!infeasible);

   SCIP_CALL( SCIPstartProbing(scip) );

   /* x <= 0.4 y + 10 */
   SCIP_CALL( addRowProbing(scip, x, y, 1.0, -0.4, -SCIPinfinity(scip), 10.0) );

   /* necessary to solve to probing LP before creating new probing nodes in solveBilinearLP() */
   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
   assert(!lperror);

   SCIP_CALL( solveBilinearLP(scip, x, y, -10, 5, 10, -5, &xcoef, &ycoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, xcoef, 1.0));
   cr_expect(SCIPisEQ(scip, ycoef, 0.4));
   cr_expect(SCIPisEQ(scip, constant, 10.0));

   SCIP_CALL( SCIPendProbing(scip) );
}

Test(projection, tworows, .init = setup, .fini=teardown,
   .description = "two row example"
   )
{
   SCIP_Bool infeasible;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real constant;
   SCIP_Bool lperror;

   /* construct LP */
   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );
   assert(!infeasible);

   SCIP_CALL( SCIPstartProbing(scip) );

   /* x <= 0.4 y + 10 */
   SCIP_CALL( addRowProbing(scip, x, y, 1.0, -0.4, -SCIPinfinity(scip), 10.0) );

   /* x <= 0.7 y + 5 */
   SCIP_CALL( addRowProbing(scip, x, y, 1.0, -0.7, -SCIPinfinity(scip), 5.0) );

   /* necessary to solve to probing LP before creating new probing nodes in solveBilinearLP() */
   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
   assert(!lperror);

   SCIP_CALL( solveBilinearLP(scip, x, y, -10, 5, 10, -5, &xcoef, &ycoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, xcoef, 1.0));
   cr_expect(SCIPisEQ(scip, ycoef, 0.7));
   cr_expect(SCIPisEQ(scip, constant, 5.0));

   SCIP_CALL( SCIPendProbing(scip) );
}

Test(projection, threevars, .init = setup, .fini=teardown,
   .description = "three variable example"
   )
{
   SCIP_Bool infeasible;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real zcoef;
   SCIP_Real constant;
   SCIP_Bool lperror;

   /* construct LP */
   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );
   assert(!infeasible);

   SCIP_CALL( SCIPsetIntParam(scip, "lp/scaling", 0) );
   SCIP_CALL( SCIPstartProbing(scip) );

   /* 2.2x <= 1.2 z + 1 */
   SCIP_CALL( addRowProbing(scip, x, z, 2.2, -1.2, -SCIPinfinity(scip), 1.0) );

   /* -0.4z >= -0.8 y + 2  <=> 2 <= 0.8y - 0.4z */
   SCIP_CALL( addRowProbing(scip, z, y, -0.4, 0.8, 2.0, SCIPinfinity(scip)) );

   /* necessary to solve to probing LP before creating new probing nodes in solveBilinearLP() */
   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
   assert(!lperror);

   SCIP_CALL( solveBilinearLP(scip, x, z, -10, 3, 10, 1, &xcoef, &zcoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, xcoef, 1.0));
   cr_expect(SCIPisEQ(scip, zcoef, 1.2/2.2));
   cr_expect(SCIPisEQ(scip, constant, 1.0/2.2));

   SCIP_CALL( solveBilinearLP(scip, z, y, 1, 5, 3, -5, &zcoef, &ycoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, zcoef, 1.0));
   cr_expect(SCIPisEQ(scip, ycoef, 2.0));
   cr_expect(SCIPisEQ(scip, constant, -5.0));

   /* the resulting inequality y >= 3 results in no useful inequality */
   SCIP_CALL( solveBilinearLP(scip, x, y, -10, 5, 10, -5, &xcoef, &ycoef, &constant, -1) );
   cr_expect(xcoef == SCIP_INVALID);
   cr_expect(ycoef == SCIP_INVALID);
   cr_expect(constant == SCIP_INVALID);

   SCIP_CALL( SCIPendProbing(scip) );
}

Test(projection, fourrows, .init = setup, .fini=teardown,
   .description = "four rows example"
   )
{
   SCIP_Bool infeasible;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real constant;
   SCIP_Bool lperror;

   /* construct LP */
   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );
   assert(!infeasible);

   SCIP_CALL( SCIPsetIntParam(scip, "lp/scaling", 0) );
   SCIP_CALL( SCIPstartProbing(scip) );

   /* TR: 3 x <= -5 y + 30  <=> -3x >= 5y - 30 */
   SCIP_CALL( addRowProbing(scip, x, y, -3.0, -5.0, -30.0, SCIPinfinity(scip)) );

   /* TL: -10 x <= -5 y + 100 */
   SCIP_CALL( addRowProbing(scip, x, y, -10.0, 5.0, -SCIPinfinity(scip), 100.0) );

   /* BR: x <= y + 8 */
   SCIP_CALL( addRowProbing(scip, x, y, 1.0, -1.0, -SCIPinfinity(scip), 8.0) );

   /* BL: -3 x <= 5 y + 40  <=> 3x >= -5y - 40 */
   SCIP_CALL( addRowProbing(scip, x, y, 3.0, 5.0, -40.0, SCIPinfinity(scip)) );

   /* necessary to solve to probing LP before creating new probing nodes in solveBilinearLP() */
   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, NULL) );
   assert(!lperror);

   /* TR */
   SCIP_CALL( solveBilinearLP(scip, x, y, -10, -5, 10, 5, &xcoef, &ycoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, xcoef, 1.0));
   cr_expect(SCIPisEQ(scip, ycoef, -5.0/3.0));
   cr_expect(SCIPisEQ(scip, constant, 10.0));

   /* TL */
   SCIP_CALL( solveBilinearLP(scip, x, y, 10, -5, -10, 5, &xcoef, &ycoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, xcoef, -1.0));
   cr_expect(SCIPisEQ(scip, ycoef, -0.5));
   cr_expect(SCIPisEQ(scip, constant, 10.0));

   /* BR */
   SCIP_CALL( solveBilinearLP(scip, x, y, -10, 5, 10, -5, &xcoef, &ycoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, xcoef, 1.0));
   cr_expect(SCIPisEQ(scip, ycoef, 1.0));
   cr_expect(SCIPisEQ(scip, constant, 8.0));

   /* BL */
   SCIP_CALL( solveBilinearLP(scip, x, y, 10, 5, -10, -5, &xcoef, &ycoef, &constant, -1) );
   cr_expect(SCIPisEQ(scip, xcoef, -1.0));
   cr_expect(SCIPisEQ(scip, ycoef, 5.0/3.0));
   cr_expect(SCIPisEQ(scip, constant, 40.0/3.0));

   SCIP_CALL( SCIPendProbing(scip) );
}
