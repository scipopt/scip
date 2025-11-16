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

/**@file   sym_symresack.c
 * @brief  symmetry handler for symresack constraints
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_symresack.h"
#include "scip/scip_sym.h"
#include "scip/sym_symresack.h"

/* symmetry handler properties */
#define SYM_NAME            "symresack"
#define SYM_DESC            "symmetry handler for symresack constraint"
#define SYM_PRIORITY          -1000000           /**< propagator priority */
#define SYM_FREQ                     1           /**< propagator frequency */


/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symhdlrTryaddSymresack)
{  /*lint --e{715}*/
   int s;

   if( symtype != SYM_SYMTYPE_PERM )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   *success = TRUE;

   for( s = 0; s < nsymmetries; ++s )
   {
      SCIP_CONS* cons;

      SCIP_CALL( SCIPcreateSymbreakCons(scip, &cons, "cons", symmetries[s], symvars, nsymvars,
            FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}


/** include symmetry handler for symresack constraints */
SCIP_RETCODE SCIPincludeSymhdlrSymresack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLR* symhdlr;

   assert(scip != NULL);

   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, &symhdlr, SYM_NAME, SYM_DESC,
         1, 1, 1, 1, -1, -1, FALSE, FALSE, -1, SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_FAST,
         symhdlrTryaddSymresack, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}
