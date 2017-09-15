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

/**@file   separation_sin.c
 * @brief  tests separation of sin()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_sin.c"
#include "separation.h"

Test(separation, sinus, .init = setup, .fini = teardown,
   .description = "test separation for a sinus expression"
   )
{
#if 0
   SCIP_CONSEXPR_EXPR* expr = NULL;
   SCIP_ROW* cut = NULL;
   int i;

   /*
    * overestimate sin(x)
    */

   /* TODO create sinus expression */

   /* add the auxiliary variable to the expression; variable will be released in CONSEXITSOL */
   SCIP_CALL( SCIPcaptureVar(scip, auxvar) );
   SCIP_CALL( SCIPaddVarLocks(scip, auxvar, 1, 1) );
   expr->auxvar = auxvar;

   /* TODO solution value of x and auxvar */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 0.0) );

   /* TODO call separation method here */
   cr_assert(cut != NULL);

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];

      /* TODO check nonzeros */
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   /*
    * underestimate sin(x)
    */

   /* TODO do the same as above, but now underestimate sin(x) */

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
#endif
}
