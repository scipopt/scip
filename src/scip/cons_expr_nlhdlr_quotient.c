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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_quotient.h
 * @brief  quotient nonlinear handler
 * @author Benjamin Mueller
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_quotient.h"
#include "scip/cons_expr.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "quotient"
#define NLHDLR_DESC         "quotient handler for expressions"
#define NLHDLR_PRIORITY     0

/*
 * Data structures
 */

/* TODO: fill in the necessary nonlinear handler data */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
};

/*
 * Local methods
 */

/* TODO: put your local methods here, and declare them static */

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrCopyhdlrQuotient NULL
#endif

/** callback to free data of handler */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrFreehdlrdataQuotient NULL
#endif


/** callback to free expression specific data */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataQuotient)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrFreeExprDataQuotient NULL
#endif


/** callback to be called in initialization */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINIT(nlhdlrInitQuotient)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitQuotient NULL
#endif


/** callback to be called in deinitialization */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLREXIT(nlhdlrExitQuotient)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrExitQuotient NULL
#endif


/** callback to detect structure in expression tree */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectQuotient)
{ /*lint --e{715}*/
   *success = FALSE;

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalauxQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSINITLP) */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitSepaQuotient NULL
#endif


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrExitSepaQuotient NULL
#endif


/** nonlinear handler under/overestimation callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrIntevalQuotient NULL
#endif


/** nonlinear handler callback for reverse propagation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReversepropQuotient NULL
#endif


/** nonlinear handler callback for reformulation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(nlhdlrReformulateQuotient)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of quotient nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReformulateQuotient NULL
#endif

/*
 * nonlinear handler specific interface methods
 */

/** includes Quotient nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrQuotient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   /* create nonlinear handler data */
   nlhdlrdata = NULL;

   /* TODO: create and store nonlinear handler specific data here */

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectQuotient, nlhdlrEvalauxQuotient, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrQuotient);
   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrFreehdlrdataQuotient);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrFreeExprDataQuotient);
   SCIPsetConsExprNlhdlrInitExit(scip, nlhdlr, nlhdlrInitQuotient, nlhdlrExitQuotient);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, nlhdlrInitSepaQuotient, NULL, nlhdlrEstimateQuotient, nlhdlrExitSepaQuotient);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalQuotient, nlhdlrReversepropQuotient);
   SCIPsetConsExprNlhdlrReformulate(scip, nlhdlr, nlhdlrReformulateQuotient);

   return SCIP_OKAY;
}
