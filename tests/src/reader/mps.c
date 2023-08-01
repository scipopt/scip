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

/**@file   mps.c
 * @brief  Unittest for MPS reader
 * @author Mathieu Besançon, Rolf van der Hulst
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/reader_mps.h"
#include "scip/scipdefplugins.h"

#include "include/scip_test.h"

static SCIP* scip;

static
void setup(void)
{
   /* create SCIP instance */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include reader */
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

}

static
void teardown(void)
{
   /* free variable and SCIP */
   SCIP_CALL( SCIPfree(&scip) );
}

/* TEST SUITE */
TestSuite(readermps, .init = setup, .fini = teardown);

Test(readermps, read1, .description = "check the function for reading a *.mps file")
{
   /*The file tested below originally failed because 'OBJSENSE MAX' is on one line. */
   char filename[SCIP_MAXSTRLEN];
   TESTsetTestfilename(filename, __FILE__, "oc5.mps");

   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );
   cr_expect( SCIPgetNVars(scip) == 5 );
}

Test(readermps, read2, .description = "check the function for reading a *.mps file")
{
    /*A 'normal' mps file with OBJSENSE MAX on two lines */
    char filename[SCIP_MAXSTRLEN];
    TESTsetTestfilename(filename, __FILE__, "oc5_1.mps");

    SCIP_CALL( SCIPreadProb(scip, filename, NULL) );
    cr_expect( SCIPgetNVars(scip) == 5 );
}
