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

/**@file   nlhdlr.h
 * @brief  private functions of nonlinear handlers of nonlinear constraints
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_NLHDLR_H_
#define SCIP_NLHDLR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "scip/type_scip.h"
#include "scip/type_nlhdlr.h"

/**@name Nonlinear Handler Methods */
/**@{ */

/** set the copy handler callback of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrCopyHdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRCOPYHDLR((*copy)) /**< copy callback (can be NULL) */
);

/** set the nonlinear handler callback to free the nonlinear handler data */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrFreeHdlrData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEHDLRDATA((*freehdlrdata)) /**< handler free callback (can be NULL) */
);

/** set the expression handler callback to free expression specific data of nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrFreeExprData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEEXPRDATA((*freeexprdata)) /**< nonlinear handler expression data free callback (can be NULL if data does not need to be freed) */
);

/** set the initialization and deinitialization callback of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrInitExit(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRINIT((*init)),   /**< initialization callback (can be NULL) */
   SCIP_DECL_NLHDLREXIT((*exit))    /**< deinitialization callback (can be NULL) */
);

/** set the propagation callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRINTEVAL((*inteval)), /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
);

/** set the enforcement callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPsetConsExprNlhdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_DECL_NLHDLRINITSEPA((*initsepa)), /**< separation initialization callback (can be NULL) */
   SCIP_DECL_NLHDLRENFO((*enfo)),         /**< enforcement callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_NLHDLRESTIMATE((*estimate)), /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_NLHDLREXITSEPA((*exitsepa))  /**< separation deinitialization callback (can be NULL) */
);

/** gives name of nonlinear handler */
SCIP_EXPORT
const char* SCIPgetConsExprNlhdlrName(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives description of nonlinear handler, can be NULL */
SCIP_EXPORT
const char* SCIPgetConsExprNlhdlrDesc(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives detection priority of nonlinear handler */
SCIP_EXPORT
int SCIPgetConsExprNlhdlrDetectPriority(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives enforcement priority of nonlinear handler */
SCIP_EXPORT
int SCIPgetConsExprNlhdlrEnfoPriority(
        SCIP_NLHDLR*      nlhdlr    /**< nonlinear handler */
);

/** returns a nonlinear handler of a given name (or NULL if not found) */
SCIP_EXPORT
SCIP_NLHDLR* SCIPfindConsExprNlhdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of nonlinear handler */
   );

/** gives handler data of nonlinear handler */
SCIP_EXPORT
SCIP_NLHDLRDATA* SCIPgetConsExprNlhdlrData(
   SCIP_NLHDLR*      nlhdlr         /**< nonlinear handler */
);

/** gives nonlinear handler expression data
 *
 * @return NULL if expr has not been detected by nlhdlr or nlhdlr did not store data
 */
SCIP_EXPORT
SCIP_NLHDLREXPRDATA* SCIPgetConsExprNlhdlrExprData(
   SCIP_NLHDLR*      nlhdlr,        /**< nonlinear handler */
   SCIP_EXPR*        expr           /**< expression */
);

/** returns whether nonlinear handler implements the interval evaluation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrInteval(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the reverse propagation callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrReverseProp(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the separation initialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrInitSepa(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the separation deinitialization callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrExitSepa(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the enforcement callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrEnfo(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the estimator callback */
SCIP_EXPORT
SCIP_Bool SCIPhasConsExprNlhdlrEstimate(
   SCIP_NLHDLR* nlhdlr              /**< nonlinear handler */
);

/** call the detect callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRDETECT(SCIPdetectConsExprNlhdlr);

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLREVALAUX(SCIPevalauxConsExprNlhdlr);

/** calls the interval evaluation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRINTEVAL(SCIPintevalConsExprNlhdlr);

/** calls the reverse propagation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRREVERSEPROP(SCIPreversepropConsExprNlhdlr);

/** calls the separation initialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRINITSEPA(SCIPinitsepaConsExprNlhdlr);

/** calls the separation deinitialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLREXITSEPA(SCIPexitsepaConsExprNlhdlr);

/** calls the enforcement callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRENFO(SCIPenfoConsExprNlhdlr);

/** calls the estimator callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRESTIMATE(SCIPestimateConsExprNlhdlr);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_NLHDLR_H_ */
