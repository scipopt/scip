/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_default.c
 * @brief  default nonlinear handler that calls expression handler methods
 * @author Stefan Vigerske
 *
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_default.h"
#include "scip/cons_expr.h"
#include "scip/struct_cons_expr.h"  // FIXME to get exprhdlr->sepa

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "default"
#define NLHDLR_DESC         "default handler for expressions"
#define NLHDLR_PRIORITY     0

static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectDefault)
{ /*lint --e{715}*/
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* return interval evaluation possibility if exprhdlr for expr has a inteval callback */
   if( SCIPgetConsExprExprHdlr(expr)->inteval != NULL )
   {
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_INTEVAL;
      *success = TRUE;
   }

   /* return reverse propagation possibility if exprhdlr for expr has a reverseprop callback */
   if( SCIPgetConsExprExprHdlr(expr)->reverseprop != NULL )
   {
      /* one could claim that reverse propagation is sufficient for enforcement, but separation is probably stronger
       * so, not setting enforcedbelow/above to TRUE here for now
       */
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_REVERSEPROP;
      *success = TRUE;
   }

   /* return sepa possibility if exprhdlr for expr has a sepa callback and enforcement is not ensured already */
   if( SCIPgetConsExprExprHdlr(expr)->sepa != NULL && (!*enforcedbelow || !*enforcedabove) )
   {
      /* make sure that an (auxiliary) variable exists for every child */
      for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
      {
         /* todo skip this for value-expressions? */
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], NULL) );
      }

      /* communicate back what the nlhdlr will do
       * - it will enforce via separation
       * - it will enforce from both below and above
       * - it needs to be called for this expression
       */
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABOTH;
      *enforcedbelow = TRUE;
      *enforcedabove = TRUE;
      *success = TRUE;
   }

   if( *success )
   {
      /* remember in the nlhdlr exprdata (pointer) which enforcemethod we advertised */
      *nlhdlrexprdata = (SCIP_CONSEXPR_NLHDLREXPRDATA*)(size_t)(*enforcemethods);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaDefault)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   assert(scip != NULL);
   assert(expr != NULL);

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   assert(exprhdlr != NULL);

   /* if we will not separate, then don't call initsepa */
   if( ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) == 0 )
      return SCIP_OKAY;

   if( exprhdlr->initsepa == NULL )
      return SCIP_OKAY;

   /* call the separation initialization callback of the expression handler */
   SCIP_CALL( exprhdlr->initsepa(scip, conshdlr, expr, infeasible) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrSepaDefault)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(ncuts != NULL);

   /* if we did not say that we will separate, then stay to it */
   if( ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) == 0 )
      return SCIP_OKAY;

   if( separated )
   {
      /* don't do anything if someone already separated */
      *result = SCIP_DIDNOTFIND;
      *ncuts = 0;

      return SCIP_OKAY;
   }

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   assert(exprhdlr != NULL);
   assert(exprhdlr->sepa != NULL);

   /* call the separation callback of the expression handler */
   SCIP_CALL( exprhdlr->sepa(scip, conshdlr, expr, sol, minviolation, result, ncuts) );

   return SCIP_OKAY;
}


static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaDefault)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   assert(scip != NULL);
   assert(expr != NULL);

   /* if we have not separated, then don't call exitsepa */
   if( ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) == 0 )
      return SCIP_OKAY;

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   assert(exprhdlr != NULL);

   if( exprhdlr->exitsepa == NULL )
      return SCIP_OKAY;

   /* call the separation deinitialization callback of the expression handler */
   SCIP_CALL( exprhdlr->exitsepa(scip, expr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalDefault)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   assert(scip != NULL);
   assert(expr != NULL);

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   assert(exprhdlr != NULL);

   if( exprhdlr->inteval == NULL )
      return SCIP_OKAY;

   assert((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_INTEVAL);

   /* call the interval evaluation callback of the expression handler */
   SCIP_CALL( exprhdlr->inteval(scip, expr, interval, varboundrelax) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropDefault)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   assert(scip != NULL);
   assert(expr != NULL);

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   assert(exprhdlr != NULL);

   if( exprhdlr->reverseprop == NULL )
      return SCIP_OKAY;

   assert((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_REVERSEPROP);

   /* call the reverse propagation callback of the expression handler */
   SCIP_CALL( exprhdlr->reverseprop(scip, expr, reversepropqueue, infeasible, nreductions, force) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrDefault)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrDefault(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}

/** includes default nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectDefault, NULL) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrDefault);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, nlhdlrInitSepaDefault, nlhdlrSepaDefault, nlhdlrExitSepaDefault);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalDefault, nlhdlrReversepropDefault);

   return SCIP_OKAY;
}
