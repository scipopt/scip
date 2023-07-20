/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   creation.c
 * @brief  unit test that checks that creation of bound disjunctions works correctly
 * @author Marc Pfetsch
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include <stdio.h>

/* TESTS  */

/* Test1: Redundancy */
Test(bounddisjunction, redundant)
{
   char name[SCIP_MAXSTRLEN];
   SCIP* scip;
   SCIP_VAR* vars[2];
   SCIP_BOUNDTYPE boundtypes[4];
   SCIP_VAR* boundvars[4];
   SCIP_Real bounds[4];
   SCIP_CONS* cons;
   int nconsvars;
   SCIP_BOUNDTYPE* consboundtypes;
   SCIP_VAR** consvars;
   int i;

   /* create an empty problem */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "bounddisjunction_redundant") );

   /* create variables */
   for (i = 0; i < 2; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
   }

   /* prepare bounds */
   boundvars[0] = vars[0];
   bounds[0] = 1.0;
   boundtypes[0] = SCIP_BOUNDTYPE_LOWER;

   boundvars[1] = vars[1];
   bounds[1] = 1.0;
   boundtypes[1] = SCIP_BOUNDTYPE_LOWER;

   /* redundant bound */
   boundvars[2] = vars[1];
   bounds[2] = 2.0;
   boundtypes[2] = SCIP_BOUNDTYPE_LOWER;

   boundvars[3] = vars[1];
   bounds[3] = 3.0;
   boundtypes[3] = SCIP_BOUNDTYPE_UPPER;

   /* create bound disjunction constraint */
   SCIP_CALL( SCIPcreateConsBasicBounddisjunctionRedundant(scip, &cons, "bounddisjunction", 4, boundvars, boundtypes, bounds) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* check that redundant literal has been eliminated */
   nconsvars = SCIPgetNVarsBounddisjunction(scip, cons);
   consvars = SCIPgetVarsBounddisjunction(scip, cons);
   consboundtypes = SCIPgetBoundtypesBounddisjunction(scip, cons);

   cr_assert( nconsvars == 3 );
   cr_assert( consvars[0] == vars[0] );
   cr_assert( consvars[1] == vars[1] );
   cr_assert( consvars[2] == vars[1] );
   cr_assert( consboundtypes[0] == SCIP_BOUNDTYPE_LOWER );
   cr_assert( consboundtypes[1] == SCIP_BOUNDTYPE_LOWER );
   cr_assert( consboundtypes[2] == SCIP_BOUNDTYPE_UPPER );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* release variables */
   for( i = 0; i < 2; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip , &vars[i]) );
   }

   /* release SCIP */
   SCIP_CALL( SCIPfree(&scip ) );
}


/* Test2: Compression 1 - all variables are fixed */
Test(bounddisjunction, compression1)
{
   SCIP* scip;
   SCIP_VAR* vars[2];
   SCIP_BOUNDTYPE boundtypes[4];
   SCIP_VAR* boundvars[4];
   SCIP_Real bounds[4];
   SCIP_CONS* cons;
   int nconsvars;
   SCIP_BOUNDTYPE* consboundtypes;
   SCIP_VAR** consvars;
   int i;

   /* create an empty problem */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "bounddisjunction_redundant") );

   /* create variables: both fixed */
   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[0], "x0", 1.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, vars[0]) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[1], "x1", 0.0, 0.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, vars[1]) );

   /* prepare bounds */
   boundvars[0] = vars[0];
   bounds[0] = 1.0;
   boundtypes[0] = SCIP_BOUNDTYPE_LOWER;

   boundvars[1] = vars[1];
   bounds[1] = 1.0;
   boundtypes[1] = SCIP_BOUNDTYPE_LOWER;

   /* redundant bound */
   boundvars[2] = vars[1];
   bounds[2] = 2.0;
   boundtypes[2] = SCIP_BOUNDTYPE_LOWER;

   boundvars[3] = vars[1];
   bounds[3] = 3.0;
   boundtypes[3] = SCIP_BOUNDTYPE_UPPER;

   /* enable compression */
   SCIPenableConsCompression(scip);

   /* create bound disjunction constraint */
   SCIP_CALL( SCIPcreateConsBasicBounddisjunctionRedundant(scip, &cons, "bounddisjunction", 4, boundvars, boundtypes, bounds) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* check that redundant literal has been eliminated */
   nconsvars = SCIPgetNVarsBounddisjunction(scip, cons);
   consvars = SCIPgetVarsBounddisjunction(scip, cons);
   consboundtypes = SCIPgetBoundtypesBounddisjunction(scip, cons);

   cr_assert( nconsvars == 1 );
   cr_assert( consvars[0] == vars[0] );
   cr_assert( consboundtypes[0] == SCIP_BOUNDTYPE_LOWER );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* release variables */
   for( i = 0; i < 2; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip , &vars[i]) );
   }

   /* release SCIP */
   SCIP_CALL( SCIPfree(&scip ) );
}


/* Test3: Compression 2 - only one variable is fixed */
Test(bounddisjunction, compression2)
{
   SCIP* scip;
   SCIP_VAR* vars[2];
   SCIP_BOUNDTYPE boundtypes[5];
   SCIP_VAR* boundvars[5];
   SCIP_Real bounds[5];
   SCIP_CONS* cons;
   int nconsvars;
   SCIP_BOUNDTYPE* consboundtypes;
   SCIP_VAR** consvars;
   int i;

   /* create an empty problem */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "bounddisjunction_redundant") );

   /* create variables: first variable fixed */
   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[0], "x0", 1.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, vars[0]) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &vars[1], "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, vars[1]) );

   /* prepare bounds */
   boundvars[0] = vars[0];
   bounds[0] = 2.0;
   boundtypes[0] = SCIP_BOUNDTYPE_LOWER;

   boundvars[1] = vars[0];
   bounds[1] = -1.0;
   boundtypes[1] = SCIP_BOUNDTYPE_UPPER;

   boundvars[2] = vars[1];
   bounds[2] = 1.0;
   boundtypes[2] = SCIP_BOUNDTYPE_LOWER;

   /* redundant bound */
   boundvars[3] = vars[1];
   bounds[3] = 2.0;
   boundtypes[3] = SCIP_BOUNDTYPE_LOWER;

   boundvars[4] = vars[1];
   bounds[4] = 3.0;
   boundtypes[4] = SCIP_BOUNDTYPE_UPPER;

   /* enable compression */
   SCIPenableConsCompression(scip);

   /* create bound disjunction constraint */
   SCIP_CALL( SCIPcreateConsBasicBounddisjunctionRedundant(scip, &cons, "bounddisjunction", 5, boundvars, boundtypes, bounds) );
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* check that redundant literal has been eliminated */
   nconsvars = SCIPgetNVarsBounddisjunction(scip, cons);
   consvars = SCIPgetVarsBounddisjunction(scip, cons);
   consboundtypes = SCIPgetBoundtypesBounddisjunction(scip, cons);

   cr_assert( nconsvars == 2 );
   cr_assert( consvars[0] == vars[1] );
   cr_assert( consvars[1] == vars[1] );
   cr_assert( consboundtypes[0] == SCIP_BOUNDTYPE_LOWER );
   cr_assert( consboundtypes[1] == SCIP_BOUNDTYPE_UPPER );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* release variables */
   for( i = 0; i < 2; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip , &vars[i]) );
   }

   /* release SCIP */
   SCIP_CALL( SCIPfree(&scip ) );
}
