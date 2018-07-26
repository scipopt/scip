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

/**@file   parsing.c
 * @brief  unit tests for printing and parsing linear constraints
 * @author Gregor Hendel
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include <stdio.h>

#define MEMSIZE 10
#define FNAME ".parsing-prob%s.cip"

/** structure to hold all necessary data to define a linear constraint for the unit tests */
struct LinData
{
   const char            names[MEMSIZE][SCIP_MAXSTRLEN]; /**< name string for every variable */
   SCIP_Real             coefs[MEMSIZE];     /**< array of linear coefficients */
   SCIP_Real             lhs;                /**< left hand side of the linear constraint */
   SCIP_Real             rhs;                /**< right hand side of the constraint */
   int                   nnonz;              /**< number of nonzero entries of the linear constraint */
};
typedef struct LinData LINDATA;

/** GLOBAL VARIABLES **/
static LINDATA lindata;
static SCIP* srcscip = NULL;
static SCIP_CONS* srccons = NULL;
static SCIP* targetscip = NULL;
static SCIP_CONS* targetcons = NULL;
static char filename[1024];

static SCIP_VAR** vars;


/* helper methods */

/** initializes an empty linear constraint with (negative) infinite (left) right hand side */
static
void initLindata(void)
{
   assert(srcscip != NULL);

   lindata.rhs = +SCIPinfinity(srcscip);
   lindata.lhs = -SCIPinfinity(srcscip);
   lindata.nnonz = 0;
}

/* TEST SUITE */

/** the setup creates the necessary source and target SCIPs, initializes linear data */
static
void setup(void)
{
   int i;
   SCIP* scips[2];
   /* set up a file pointer to read from, and the src and target SCIP */

   /* initialization loop over SCIPs */
   for( i = 0; i < 2; ++i )
   {
      SCIP_CALL( SCIPcreate(&(scips[i])) );
      SCIP_CALL( SCIPincludeConshdlrLinear(scips[i]) );
      SCIP_CALL( SCIPincludeReaderCip(scips[i]) );
   }

   srcscip = scips[0];
   targetscip = scips[1];
   SCIP_CALL (SCIPcreateProbBasic(srcscip, "parsing-prob"));

   SCIP_CALL( SCIPallocMemoryArray(srcscip, &vars, MEMSIZE) );

   initLindata();
}

/** free all allocated memory */
static
void teardown(void)
{
   int v;
   /** free the SCIPs */
   SCIPfree(&targetscip);


   if( srccons != NULL )
   {
      SCIP_CALL( SCIPreleaseCons(srcscip, &srccons) );
   }

   /* release the variables, otherwise, bad things happen */
   for( v = 0; v < lindata.nnonz; ++v )
   {
      SCIP_CALL( SCIPreleaseVar(srcscip, &(vars[v])) );
   }
   SCIPfreeMemoryArray(srcscip, &vars);

   SCIPfree(&srcscip);

   /* remove the cip-file for this test */
   if( strncmp(filename, FNAME, 5)  == 0 )
      remove(filename);

}

/** compares a linear data structure to a given constraint */
static
void compareConstraint(
   SCIP*                 scip,               /**< source or target scip */
   SCIP_CONS*            cons                /**< linear constraint of the source or target SCIP instance */
   )
{
   int i;

   cr_assert_eq(lindata.nnonz, SCIPgetNVarsLinear(scip, cons), "Constraint number of variables != %d", lindata.nnonz);

   /* test variable names, coefficients, etc */
   for( i = 0; i < lindata.nnonz; ++i )
   {
      /* test for epsilon equality of every coefficient */
      cr_assert_float_eq(lindata.coefs[i], SCIPgetValsLinear(scip, cons)[i], SCIPepsilon(scip), "Unequal coefficients %.15g != %.15g",
            lindata.coefs[i], SCIPgetValsLinear(scip, cons)[i]);

      /* test all variable names */
      cr_assert_str_eq(lindata.names[i], SCIPvarGetName(SCIPgetVarsLinear(scip, cons)[i]), "Variable names are different: %s != %s",
            lindata.names[i], SCIPvarGetName(SCIPgetVarsLinear(scip, cons)[i]));
   }


}

/** create */
static
void createSourceProblemCons(void)
{
   int i;
   /* create variables and add them to SCIP */
   for( i = 0; i < lindata.nnonz; ++i )
   {
      SCIP_CALL( SCIPcreateVarBasic(srcscip, &vars[i], lindata.names[i], 0, SCIPinfinity(srcscip), 0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(srcscip, vars[i]) );
   }

   SCIP_CALL( SCIPcreateConsBasicLinear(srcscip, &srccons, "source-constraint", lindata.nnonz,
         vars, lindata.coefs, lindata.lhs, lindata.rhs) );

   SCIP_CALL( SCIPaddCons(srcscip, srccons) );
   compareConstraint(srcscip, srccons);
}

static
void writeSourceProblem(void)
{
   SCIP_CALL( SCIPwriteOrigProblem(srcscip, filename, "cip", FALSE) );
}

static
void readIntoTargetProblem(void)
{
   SCIP_CALL( SCIPreadProb(targetscip, filename, "cip") );

   cr_assert_eq(SCIPgetNConss(targetscip), 1, "Number %d of target SCIP constraints should be 1", SCIPgetNConss(targetscip));
   targetcons = SCIPgetConss(targetscip)[0];

   compareConstraint(targetscip, targetcons);
}

static
void performTestOnLinData(void)
{
   createSourceProblemCons();
   writeSourceProblem();
   readIntoTargetProblem();
}

TestSuite(parsing_linear, .init = setup, .fini = teardown);

/* TESTS  */
Test(parsing_linear, create_and_free)
{
   /* calls setup and teardown */
}

Test(parsing_linear, equation_0, .description = "test if an equation 0 == 0 is correctly handled")
{
   lindata.lhs = lindata.rhs = 0.0;
   sprintf(filename, FNAME, "equation_0");

   performTestOnLinData();
}

Test(parsing_linear, equation_const, .description = "test if a ranged row -const <= const is correctly handled")
{
   lindata.lhs = -5.0;
   lindata.rhs = 5.0;
   sprintf(filename, FNAME, "equation_const");

   performTestOnLinData();
}

Test(parsing_linear, inf_lhs, .description = "test if an infinite left hand side is correctly handled")
{
   lindata.rhs = 5.0;
   sprintf(filename, FNAME, "inf_lhs");

   performTestOnLinData();
}

Test(parsing_linear, inf_rhs, .description = "test if an infinite right hand side is correctly handled")
{
   lindata.lhs = 0.0;
   sprintf(filename, FNAME, "inf_rhs");

   performTestOnLinData();
}

Test(parsing_linear, equation_onevariable, .description = "test if an equation const == const * variable is correctly handled")
{
   lindata.lhs = 5.0;
   lindata.rhs = 5.0;
   lindata.coefs[0] = 5.0;
   BMScopyMemoryArray(lindata.names[0], "v1", 3);
   lindata.nnonz = 1;

   sprintf(filename, FNAME, "equation_onevariable");

   performTestOnLinData();
}

Test(parsing_linear, free_constraint, .description = "test if a free row is correctly handled")
{
   /* left and right hand side are initialized to infinity */
   lindata.coefs[0] = 5.0;
   BMScopyMemoryArray(lindata.names[0], "v1", 3);
   lindata.nnonz = 1;

   sprintf(filename, FNAME, "free_constraint");

   performTestOnLinData();
}

Test(parsing_linear, infeasible_constraint, .description = "test infeasible constraint with left hand side larger than right hand side")
{
   lindata.lhs = 10.0;
   lindata.rhs = 0.0;
   lindata.coefs[0] = 5.0;
   BMScopyMemoryArray(lindata.names[0], "v1", 3);
   lindata.nnonz = 1;

   sprintf(filename, FNAME, "infeasible_constraint");

   performTestOnLinData();
}
