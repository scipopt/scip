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

/**@file   nlhdlr_signomial.h
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  signomial nonlinear handler
 * @author Liding Xu
 */

#define SCIP_DEBUG
#include <string.h>

#include "nlhdlr_signomial.h"
#include "scip/cons_nonlinear.h"
#include "scip/pub_misc_rowprep.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "signomial"
#define NLHDLR_DESC               "signomial handler for expressions"
#define NLHDLR_DETECTPRIORITY     -15  /**< it is important that the nlhdlr runs after the default nldhlr */
#define NLHDLR_ENFOPRIORITY       -10

#define MIN_INTERIORITY           0.01 /**< minimum interiority for a reference point for applying separation */
#define MIN_ABSBOUNDSIZE          0.1  /**< minimum size of variable bounds for applying separation */


/*
 * Data structures
 */

/** nonlinear handler signomial expression data. A signomial expression is cx^e associated with an auxiliary variable y. We assume x^e  = t = y/c. After reformulation,
u^f = v^g, the last of right v is t*/
struct SCIP_NlhdlrExprData{
   SCIP_EXPR*            expr;              /**< expression */
   SCIP_Real            coef;              /**< coeffcient c */
   SCIP_VAR*             auxvar;              /**< auxiliary variable y */
   SCIP_EXPR**            factors;              /**< factors x */
   int                  nfactors;             /**< number of factors x */
   SCIP_Real*           exponents;         /**< exponents a */
   SCIP_EXPR**            lfactors;             /**< left factors u */
   int                   nlfactors;            /**< number of left factors u */   
   SCIP_Real*           lexponents;        /**< left exponents f */
   SCIP_EXPR**            rfactors;             /**< right factors v */
   int                   nrfactors;            /**< number of right factors v */   
   SCIP_Real*           rexponents;        /**< right exponents g */
};


/** nonlinear handler data */
struct SCIP_NlhdlrData
{
   int nsignomials;
   /* parameters */
};


/*
 * Local methods
 */




/* TODO: put your local methods here, and declare them static */

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrSignomial)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrSignomial(targetscip) );

   return SCIP_OKAY;
}



/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataSignomial)
{  /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   int pos;

   assert(expr != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);


   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->factors, (*nlhdlrexprdata)->nfactors);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->exponents, (*nlhdlrexprdata)->nfactors);

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->lfactors, (*nlhdlrexprdata)->nlfactors);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->lexponents, (*nlhdlrexprdata)->nlfactors);

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->rfactors, (*nlhdlrexprdata)->nrfactors);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->rexponents, (*nlhdlrexprdata)->nrfactors);


   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}




/** callback to detect structure in expression tree */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectSignomial)
{ /*lint --e{715}*/

   SCIP_NLHDLRDATA* nlhdlrdata;


   assert(expr != NULL);
   assert(participating != NULL);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);


   /* check for product expressions with more than one children */
   if(SCIPisExprProduct(scip, expr) && SCIPexprGetNChildren(expr) >= 2)
   {

      int c;
      int nc;

      nc = SCIPexprGetNChildren(expr);

      /* create expression data for the nonlinear handler */
      SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );

      SCIPallocBlockMemoryArray(scip, & (*nlhdlrexprdata)->factors, nc);
      SCIPallocBlockMemoryArray(scip, & (*nlhdlrexprdata)->exponents, nc);


      SCIP_CALL( SCIPgetExprMonomialData(scip, expr, &((*nlhdlrexprdata)->coef), (*nlhdlrexprdata)->exponents, (*nlhdlrexprdata)->factors));

      /* add expression to nlhdlrdata and capture it */
      (*nlhdlrexprdata)->expr = expr;
      SCIPcaptureExpr(expr);

      /* auxuliary variable will be detected later */
      (*nlhdlrexprdata)->auxvar = NULL;
      (*nlhdlrexprdata)->nfactors = nc;

      int nlexponents = 0;
      SCIP_Real sumlexponents = 0;
      SCIP_Real sumrexponents = 1;

      /* count positive exponent terms  */
      for(c = 0; c < nc; c++){
         if(SCIPisPositive(scip, (*nlhdlrexprdata)->exponents[c]) ){
            sumlexponents += (*nlhdlrexprdata)->exponents[c];
            nlexponents++;
         }
         else{
            sumrexponents += (*nlhdlrexprdata)->exponents[c];
         }
      }

      (*nlhdlrexprdata)->nlfactors = nlexponents;
      SCIPallocBlockMemoryArray(scip, & (*nlhdlrexprdata)->lfactors, (*nlhdlrexprdata)->nlfactors);
      SCIPallocBlockMemoryArray(scip, & (*nlhdlrexprdata)->lexponents, (*nlhdlrexprdata)->nlfactors);

      (*nlhdlrexprdata)->nrfactors = nc + 1 - nlexponents;
      SCIPallocBlockMemoryArray(scip, & (*nlhdlrexprdata)->rfactors, (*nlhdlrexprdata)->nrfactors - 1);
      SCIPallocBlockMemoryArray(scip, & (*nlhdlrexprdata)->rexponents, (*nlhdlrexprdata)->nrfactors);

      /* compute and capture left and right exponent terms */
      int lc = 0, rc = 0;
      for(c = 0; c < nc; c++){
         SCIPcaptureExpr((*nlhdlrexprdata)->factors[c]);
         if(SCIPisPositive(scip, (*nlhdlrexprdata)->exponents[c]) ){
            (*nlhdlrexprdata)->lfactors[lc] = (*nlhdlrexprdata)->factors[c];
            (*nlhdlrexprdata)->lexponents[lc] = (*nlhdlrexprdata)->exponents[c] / sumlexponents;
            lc++;
         }
         else{
            (*nlhdlrexprdata)->rfactors[rc] = (*nlhdlrexprdata)->factors[c];
            (*nlhdlrexprdata)->rexponents[rc] = (*nlhdlrexprdata)->exponents[c] / sumrexponents;
            rc++;
         }
      }

      /* tell children that we will use their auxvar and use its activity for both estimate and domain propagation */
      for(c = 0; c < nc; c++ )
      {
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, (*nlhdlrexprdata)->factors[c], TRUE, FALSE, TRUE, TRUE) );
      }

   }

   if( *nlhdlrexprdata != NULL )
   {
      /* we want to join separation, if not disabled by parameter */
      *participating = SCIP_NLHDLR_METHOD_SEPABOTH;
   }

#ifdef SCIP_DEBUG
   if( *participating )
   {
      SCIPdebugMsg(scip, "detected expr: total vars (exps) %d, left vars (exps) %d, right vars (exps) %d,  ", (*nlhdlrexprdata)->nfactors, (*nlhdlrexprdata)->nlfactors, (*nlhdlrexprdata)->nrfactors);
      SCIPprintExpr(scip, expr, NULL);
      SCIPinfoMessage(scip, NULL, " participating: %d\n", *participating);
   }
#endif

   return SCIP_OKAY;
}

/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxSignomial)
{ /*lint --e{715}*/
   SCIP_VAR * var;
   SCIP_Real val;
   int c;

   *auxvalue = nlhdlrexprdata->coef;
   for( c = 0; c < nlhdlrexprdata->nfactors; ++c )
   {
      assert(nlhdlrexprdata->factors[c] != NULL);
      var =  SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->factors[c]);
      assert(var != NULL);
      val = SCIPgetSolVal(scip, sol, var);
      *auxvalue *= pow(val, nlhdlrexprdata->exponents[c]);
   }

   return SCIP_OKAY;
}


/** nonlinear handler to free data of handler */
#if 0
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataSignomial)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);
   assert((*nlhdlrdata)->nsignomials == 0);
   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}
#else
#define nlhdlrFreehdlrdataSignomial NULL
#endif


/** callback to free data of handler */



/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_NLHDLRINTEVAL(nlhdlrIntevalSignomial)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrIntevalSignomial NULL
#endif



/** nonlinear handler callback for reverse propagation */
#if 0
static
SCIP_DECL_NLHDLRREVERSEPROP(nlhdlrReversepropSignomial)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReversepropSignomial NULL
#endif



/** nonlinear handler enforcement callback */
#if 0
static
SCIP_DECL_NLHDLRENFO(nlhdlrEnfoSignomial)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrEnfoSignomial NULL
#endif


/** nonlinear handler under/overestimation callback */
#if 0
static
SCIP_DECL_NLHDLRESTIMATE(nlhdlrEstimateSignomial)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of xyz nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrEstimateSignomial NULL
#endif


/*
 * nonlinear handler specific interface methods
 */

/** includes Signomial nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeNlhdlrSignomial(
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
      NLHDLR_ENFOPRIORITY, nlhdlrDetectSignomial, nlhdlrEvalauxSignomial, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrSignomial);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrdataSignomial);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataSignomial);
   SCIPnlhdlrSetSepa(nlhdlr, NULL, nlhdlrEnfoSignomial, nlhdlrEstimateSignomial, NULL);
   SCIPnlhdlrSetProp(nlhdlr, nlhdlrIntevalSignomial, nlhdlrReversepropSignomial);

   return SCIP_OKAY;
}
