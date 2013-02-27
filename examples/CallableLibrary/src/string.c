/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   string.c
 * @brief  Coil Compression String Design model
 * @author Stefan Vigerske
 *
 * This example shows how to setup quadratic and nonlinear constraints in SCIP when using SCIP as callable library.
 * The example implements a model for the design of a coil compression string as it can be found in the GAMS model library:
 * http://www.gams.com/modlib/libhtml/spring.htm
 *
 * The task is to find a minimum volume of a wire for the production of a coil compression spring.
 *
 * Original model source:
 * @par
 *    E. Sangren@n
 *    Nonlinear Integer and Discrete Programming in Mechanical Design Optimization@n
 *    Journal of Mechanical Design, Trans. ASME 112 (1990), 223-229
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <math.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif


/* Model parameters */

/** number of possible wire types */
#define nwires 11

/** diameters of available diameters (in) */
static const SCIP_Real diameters[] = { 0.207, 0.225, 0.244, 0.263, 0.283, 0.307, 0.331, 0.362, 0.394, 0.4375, 0.500 };

/** preload (lb) */
static const SCIP_Real preload = 300;

/** maximal working load (lb) */
static const SCIP_Real maxworkload = 1000;

/** maximal deflection (in) */
static const SCIP_Real maxdeflect = 6;

/** deflection from preload (in) */
static const SCIP_Real deflectpreload = 1.25;

/** maximal free length of spring (in) */
static const SCIP_Real maxfreelen = 14.0;

/** maximal coil diameter (in) */
static const SCIP_Real maxcoildiam = 3.0;

/** maximal shear stress */
static const SCIP_Real maxshearstress = 189000.0;

/** shear modulus of material */
static const SCIP_Real shearmod = 11500000.0;


/** sets up problem */
static
SCIP_RETCODE setupProblem(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR* coil;        /* coil diameter */
   SCIP_VAR* wire;        /* wire diameter */
   SCIP_VAR* defl;        /* deflection */
   SCIP_VAR* ncoils;      /* number of coils (integer) */
   SCIP_VAR* const1;      /* a constant */
   SCIP_VAR* const2;      /* another constant */
   SCIP_VAR* volume;      /* total volume */
   SCIP_VAR* y[nwires];   /* wire choice (binary) */

   SCIP_CONS* voldef;
   SCIP_CONS* defconst1;
   SCIP_CONS* defconst2;
   SCIP_CONS* shear;
   SCIP_CONS* defdefl;
   SCIP_CONS* freel;
   SCIP_CONS* coilwidth;
   SCIP_CONS* defwire;
   SCIP_CONS* selectwire;

   char name[SCIP_MAXSTRLEN];
   int i;

   /* create empty problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "string") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &coil, "coildiam", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &wire, "wirediam", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &defl, "deflection", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &ncoils, "ncoils", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &const1, "const1", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &const2, "const2", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &volume, "volume", 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
   for( i = 0; i < nwires; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "wire%d", i+1);
      SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], name, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   }

   /* set nonstandard variable bounds */
   SCIP_CALL( SCIPchgVarLb(scip, defl, deflectpreload / (maxworkload - preload)) );
   SCIP_CALL( SCIPchgVarUb(scip, defl, maxdeflect / preload) );

   /* add variables to problem */
   SCIP_CALL( SCIPaddVar(scip, coil) );
   SCIP_CALL( SCIPaddVar(scip, wire) );
   SCIP_CALL( SCIPaddVar(scip, defl) );
   SCIP_CALL( SCIPaddVar(scip, ncoils) );
   SCIP_CALL( SCIPaddVar(scip, const1) );
   SCIP_CALL( SCIPaddVar(scip, const2) );
   SCIP_CALL( SCIPaddVar(scip, volume) );
   for( i = 0; i < nwires; ++i )
   {
      SCIP_CALL( SCIPaddVar(scip, y[i]) );
   }

   /* nonlinear constraint voldef: PI/2 * (ncoils+2)*coil*wire^2 - volume == 0 */
   {
      SCIP_EXPR* ncoilsplus2;
      SCIP_EXPR* coilexpr;
      SCIP_EXPR* wireexpr;
      SCIP_EXPR* expr;
      SCIP_EXPR* children[3];
      SCIP_Real exponents[3] = { 1.0, 1.0, 2.0 };
      SCIP_VAR* vars[3];
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_EXPRTREE* exprtree;
      SCIP_Real one;
      SCIP_Real minusone;

      one = 1.0;
      minusone = -1.0;

      /* setup expression tree for PI/2 * (N+2)*coil*wire^2
       * in the expression tree, we relate the variable indices as follows:
       *   0: ncoils
       *   1: coil
       *   2: wire
       */

      /* setup expression for ncoils+2 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &ncoilsplus2, SCIP_EXPR_VARIDX, 0) );
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &ncoilsplus2, 1, &ncoilsplus2, &one, 2.0) );

      /* setup expression for variable coil */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &coilexpr, SCIP_EXPR_VARIDX, 1) );

      /* setup expression for variable wire */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &wireexpr, SCIP_EXPR_VARIDX, 2) );

      /* setup monomial for PI/2 * ncoilsplus2 * coilexpr * wireexpr^2 */
      SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, M_PI / 2.0, 3, NULL, exponents) );

      /* setup polynomial expression for only one monomial
       * (FALSE as last argument indicates that the polynomial assumes ownership of monomial)
       */
      children[0] = ncoilsplus2;
      children[1] = coilexpr;
      children[2] = wireexpr;
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), &expr, 3, children, 1, &monomial, 0.0, FALSE) );

      /* setup expression tree with expr as root expression, the tree is defined w.r.t. 3 variables */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 3, 0, NULL) );

      /* assign SCIP variables to tree */
      vars[0] = ncoils;
      vars[1] = coil;
      vars[2] = wire;
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, 3, vars) );

      /* create nonlinear constraint for exprtree - volume = 0.0 */
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &voldef, "voldef", 1, &volume, &minusone, 1, &exprtree, &one, 0.0, 0.0) );

      /* free expression tree, because it was copied by the constraint */
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   /* nonlinear constraint defconst1: coil / wire - const1 == 0.0 */
   {
      SCIP_EXPR* coilexpr;
      SCIP_EXPR* wireexpr;
      SCIP_EXPR* expr;
      SCIP_EXPRTREE* exprtree;
      SCIP_VAR* vars[2];

      SCIP_Real one;
      SCIP_Real minusone;

      one = 1.0;
      minusone = -1.0;

      /* expression for variables coilexpr (index 0) and wireexpr (index 1) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &coilexpr, SCIP_EXPR_VARIDX, 0) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &wireexpr, SCIP_EXPR_VARIDX, 1) );

      /* expression for coil / wire */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, coilexpr, wireexpr) );

      /* expression tree from expr */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 2, 0, NULL) );

      /* set variables in expression tree */
      vars[0] = coil;
      vars[1] = wire;
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, 2, vars) );

      /* create nonlinear constraint for exprtree - const1 = 0.0 */
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &defconst1, "defconst1", 1, &const1, &minusone, 1, &exprtree, &one, 0.0, 0.0) );

      /* free expression tree, because it was copied by the constraint */
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   /* nonlinear constraint defconst2: (4.0 * const1 - 1.0) / (4.0 * const1 - 4.0) + 0.615 / const1 - const2 == 0.0 */
   {
      SCIP_EXPR* const1expr;
      SCIP_EXPR* denom;
      SCIP_EXPR* nomin;
      SCIP_EXPR* expr;
      SCIP_EXPRTREE* exprtrees[2];

      SCIP_Real coef;
      SCIP_Real coefs[2];
      SCIP_Real minusone;

      minusone = -1.0;

      /* expression for denominator 4.0 * const1 - 1.0 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &const1expr, SCIP_EXPR_VARIDX, 0) );
      coef = 4.0;
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &denom, 1, &const1expr, &coef, -1.0) );

      /* expression for nominator 4.0 * const1 - 4.0
       * (we cannot reuse const1expr a second time in the expression tree, thus we create a new expression for this variable)
       */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &const1expr, SCIP_EXPR_VARIDX, 0) );
      coef = 4.0;
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &nomin, 1, &const1expr, &coef, -4.0) );

      /* expression for quotient (4.0 * const1 - 1.0) / (4.0 * const1 - 4.0) */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, denom, nomin) );

      /* expression tree from expr */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtrees[0], expr, 1, 0, NULL) );
      SCIP_CALL( SCIPexprtreeSetVars(exprtrees[0], 1, &const1) );


      /* expression for 1.0 / const1 */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &const1expr, SCIP_EXPR_VARIDX, 0) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_INTPOWER, const1expr, -1) );

      /* expression tree from expr */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtrees[1], expr, 1, 0, NULL) );
      SCIP_CALL( SCIPexprtreeSetVars(exprtrees[1], 1, &const1) );


      /* create nonlinear constraint for exprtree[0] + 0.615 * exprtree[1] - const2 = 0.0 */
      coefs[0] = 1.0;
      coefs[1] = 0.615;
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &defconst2, "defconst2", 1, &const2, &minusone, 2, exprtrees, coefs, 0.0, 0.0) );


      /* free expression trees, because they were copied by the constraint */
      SCIP_CALL( SCIPexprtreeFree(&exprtrees[0]) );
      SCIP_CALL( SCIPexprtreeFree(&exprtrees[1]) );
   }

   /* quadratic constraint shear: 8.0*maxworkload/PI * const1 * const2 - maxshearstress * wire^2 <= 0.0 */
   {
      /* create empty quadratic constraint with right-hand-side 0.0 */
      SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &shear, "shear", 0, NULL, NULL, 0, NULL, NULL, NULL, -SCIPinfinity(scip), 0.0) );

      /* add bilinear term 8.0*maxworkload/PI * const1 * const2 */
      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, shear, const1, const2, 8.0 * maxworkload / M_PI) );

      /* add square term -maxshearstress * wire^2 */
      SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, shear, wire, -maxshearstress) );
   }

   /* nonlinear constraint defdefl: 8.0/shearmod * ncoils * const1^3 / wire - defl == 0.0 */
   {
      SCIP_EXPR* expr;
      SCIP_EXPR* children[3];
      SCIP_Real exponents[3] = { 1.0, 3.0, -1.0 };
      SCIP_VAR* vars[3];
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_EXPRTREE* exprtree;
      SCIP_Real one;
      SCIP_Real minusone;

      one = 1.0;
      minusone = -1.0;

      /* we relate the variable indices as follows:
       *  0: ncoils
       *  1: const1
       *  2: wire
       */

      /* setup expressions for ncoils, const1, and wire */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[0], SCIP_EXPR_VARIDX, 0) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[1], SCIP_EXPR_VARIDX, 1) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &children[2], SCIP_EXPR_VARIDX, 2) );

      /* setup monomial for 8.0/shearmod * ncoils * const1^3 / wire */
      SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, 8.0 / shearmod, 3, NULL, exponents) );

      /* setup polynomial expression for only one monomial */
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), &expr, 3, children, 1, &monomial, 0.0, FALSE) );

      /* setup expression tree with expr as root expression */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 3, 0, NULL) );

      /* assign SCIP variables to tree */
      vars[0] = ncoils;
      vars[1] = const1;
      vars[2] = wire;
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, 3, vars) );

      /* create nonlinear constraint for exprtree - defl = 0.0 */
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &defdefl, "defdefl", 1, &defl, &minusone, 1, &exprtree, &one, 0.0, 0.0) );

      /* free expression tree, because it was copied by the constraint */
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   /* quadratic constraint freel: maxworkload * defl + 1.05 * ncoils * wire + 2.1 * wire <= maxfreelen */
   {
      SCIP_Real one05;

      /* create quadratic constraint maxworkload * defl + 1.05 * ncoils * wire <= maxfreelen */
      one05 = 1.05;
      SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &freel, "freel",
         1, &defl, (SCIP_Real*)&maxworkload,
         1, &ncoils, &wire, &one05, -SCIPinfinity(scip), maxfreelen) );

      /* add linear term 2.1 * wire for variable wire in quadratic part of constraint */
      SCIP_CALL( SCIPaddQuadVarLinearCoefQuadratic(scip, freel, wire, 2.1) );
   }

   /* linear constraint coilwidth: coil + wire <= maxcoildiam */
   {
      /* create empty linear constraint with right-hand-side maxcoildiam */
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &coilwidth, "coilwidth", 0, NULL, NULL, -SCIPinfinity(scip), maxcoildiam) );

      /* add linear term coil + wire */
      SCIP_CALL( SCIPaddCoefLinear(scip, coilwidth, coil, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, coilwidth, wire, 1.0) );
   }

   /* linear constraint defwire: sum_i b[i]*y[i] - wire == 0.0 */
   {
      /* create linear constraint sum_i b[i]*y[i] == 0.0 */
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &defwire, "defwire", nwires, y, (SCIP_Real*)diameters, 0.0, 0.0) );

      /* add term -wire */
      SCIP_CALL( SCIPaddCoefLinear(scip, defwire, wire, -1.0) );
   }

   /* specialized linear constraint selectwire: sum_i y[i] == 1.0 */
   {
      SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &selectwire, "selectwire", nwires, y) );
   }


   /* add constraints to problem */
   SCIP_CALL( SCIPaddCons(scip, voldef) );
   SCIP_CALL( SCIPaddCons(scip, defconst1) );
   SCIP_CALL( SCIPaddCons(scip, defconst2) );
   SCIP_CALL( SCIPaddCons(scip, shear) );
   SCIP_CALL( SCIPaddCons(scip, defdefl) );
   SCIP_CALL( SCIPaddCons(scip, freel) );
   SCIP_CALL( SCIPaddCons(scip, coilwidth) );
   SCIP_CALL( SCIPaddCons(scip, defwire) );
   SCIP_CALL( SCIPaddCons(scip, selectwire) );


   /* release variables and constraints
    * the problem has them captured, and we do not require them anymore
    */
   SCIP_CALL( SCIPreleaseVar(scip, &coil) );
   SCIP_CALL( SCIPreleaseVar(scip, &wire) );
   SCIP_CALL( SCIPreleaseVar(scip, &defl) );
   SCIP_CALL( SCIPreleaseVar(scip, &ncoils) );
   SCIP_CALL( SCIPreleaseVar(scip, &const1) );
   SCIP_CALL( SCIPreleaseVar(scip, &const2) );
   SCIP_CALL( SCIPreleaseVar(scip, &volume) );
   for( i = 0; i < nwires; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
   }

   SCIP_CALL( SCIPreleaseCons(scip, &voldef) );
   SCIP_CALL( SCIPreleaseCons(scip, &defconst1) );
   SCIP_CALL( SCIPreleaseCons(scip, &defconst2) );
   SCIP_CALL( SCIPreleaseCons(scip, &shear) );
   SCIP_CALL( SCIPreleaseCons(scip, &defdefl) );
   SCIP_CALL( SCIPreleaseCons(scip, &freel) );
   SCIP_CALL( SCIPreleaseCons(scip, &coilwidth) );
   SCIP_CALL( SCIPreleaseCons(scip, &defwire) );
   SCIP_CALL( SCIPreleaseCons(scip, &selectwire) );

   return SCIP_OKAY;
}

/* runs string example */
static
SCIP_RETCODE runString(void)
{
   SCIP* scip;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "************************************************\n");
   SCIPinfoMessage(scip, NULL, "* Running Coil Compression String Design Model *\n");
   SCIPinfoMessage(scip, NULL, "************************************************\n");
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( setupProblem(scip) );

   SCIPinfoMessage(scip, NULL, "Original problem:\n");
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPpresolve(scip) );

   /* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
   */

   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPfreeTransform(scip) );

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
      SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
   }

   SCIP_CALL( SCIPfree(&scip) );

   return SCIP_OKAY;
}


/** main method starting SCIP */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )  /*lint --e{715}*/
{
   SCIP_RETCODE retcode;

   retcode = runString();

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
