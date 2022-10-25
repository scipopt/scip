/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_xyz.h
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  xyz nonlinear handler
 * @author Jane Doe
 */

#include <string.h>

#include "scip/nlhdlr_xyz.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc_rowprep.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "xyz"
#define NLHDLR_DESC               "xyz handler for expressions"
#define NLHDLR_DETECTPRIORITY     0
#define NLHDLR_ENFOPRIORITY       0

/*
 * Data structures
 */

/* TODO: fill in the necessary nonlinear handler data */

/** nonlinear handler data */
struct SCIP_NlhdlrData
{
};

/** nonlinear handler expression data */
struct SCIP_NlhdlrExprData
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
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrCopyhdlrXyz NULL
#endif

/** callback to free data of handler */
#if 0
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrFreehdlrdataXyz NULL
#endif


/** callback to free expression specific data */
#if 0
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrFreeExprDataXyz NULL
#endif


/** callback to be called in initialization */
#if 0
static
SCIP_DECL_NLHDLRINIT(nlhdlrInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitXyz NULL
#endif


/** callback to be called in deinitialization */
#if 0
static
SCIP_DECL_NLHDLREXIT(nlhdlrExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrExitXyz NULL
#endif


/** callback to detect structure in expression tree */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectXyz)
{ /*lint --e{715}*/
   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** separation initialization method of a nonlinear handler (called during CONSINITLP) */
#if 0
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitSepaXyz NULL
#endif


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
#if 0
static
SCIP_DECL_NLHDLREXITSEPA(nlhdlrExitSepaXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrExitSepaXyz NULL
#endif


/** nonlinear handler enforcement callback */
#if 0
static
SCIP_DECL_NLHDLRENFO(nlhdlrEnfoXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrEnfoXyz NULL
#endif


/** nonlinear handler under/overestimation callback */
#if 0
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrEstimateXyz NULL
#endif


/** nonlinear handler solution linearization callback */
#if !1
static
SCIP_DECL_NLHDLRSOLLINEARIZE(nlhdlrSollinearizeXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrSollinearizeXyz NULL
#endif


/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_NLHDLRINTEVAL(nlhdlrIntevalXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrIntevalXyz NULL
#endif


/** nonlinear handler callback for reverse propagation */
#if 0
static
SCIP_DECL_NLHDLRREVERSEPROP(nlhdlrReversepropXyz)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReversepropXyz NULL
#endif


/*
 * nonlinear handler specific interface methods
 */

/** includes Xyz nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeNlhdlrXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   /* create nonlinear handler data */
   nlhdlrdata = NULL;

   /* TODO: create and store nonlinear handler specific data here */

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectXyz, nlhdlrEvalauxXyz, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrXyz);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrdataXyz);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataXyz);
   SCIPnlhdlrSetInitExit(nlhdlr, nlhdlrInitXyz, nlhdlrExitXyz);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaXyz, nlhdlrEnfoXyz, nlhdlrEstimateXyz, nlhdlrExitSepaXyz);
   SCIPnlhdlrSetProp(nlhdlr, nlhdlrIntevalXyz, nlhdlrReversepropXyz);
   SCIPnlhdlrSetSollinearize(nlhdlr, nlhdlrSollinearizeXyz);

   return SCIP_OKAY;
}
