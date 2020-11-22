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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_nonlinear.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for nonlinear functions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_NONLINEAR_H__
#define __SCIP_SCIP_NONLINEAR_H__


#include "nlpi/type_nlpi.h"
#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_misc.h"
#include "scip/type_nlp.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicNonlinearMethods
 *
 * @{
 */

/** computes coefficients of linearization of a square term in a reference point */
SCIP_EXPORT
void SCIPaddSquareLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             refpoint,           /**< point where to linearize */
   SCIP_Bool             isint,              /**< whether corresponding variable is a discrete variable, and thus linearization could be moved */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   );

/** computes coefficients of secant of a square term */
SCIP_EXPORT
void SCIPaddSquareSecant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             sqrcoef,            /**< coefficient of square term */
   SCIP_Real             lb,                 /**< lower bound on variable */
   SCIP_Real             ub,                 /**< upper bound on variable */
   SCIP_Real             refpoint,           /**< point for which to compute value of linearization */
   SCIP_Real*            lincoef,            /**< buffer to add coefficient of secant */
   SCIP_Real*            linconstant,        /**< buffer to add constant of secant */
   SCIP_Bool*            success             /**< buffer to set to FALSE if secant has failed due to large numbers or unboundedness */
   );

/** computes coefficients of linearization of a bilinear term in a reference point */
SCIP_EXPORT
void SCIPaddBilinLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             refpointx,          /**< point where to linearize first  variable */
   SCIP_Real             refpointy,          /**< point where to linearize second variable */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   );

/** computes coefficients of McCormick under- or overestimation of a bilinear term */
SCIP_EXPORT
void SCIPaddBilinMcCormick(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real*            lincoefx,           /**< buffer to add coefficient of first  variable in linearization */
   SCIP_Real*            lincoefy,           /**< buffer to add coefficient of second variable in linearization */
   SCIP_Real*            linconstant,        /**< buffer to add constant of linearization */
   SCIP_Bool*            success             /**< buffer to set to FALSE if linearization has failed due to large numbers */
   );

/** computes coefficients of linearization of a bilinear term in a reference point when given a linear inequality
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 */
SCIP_EXPORT
void SCIPcomputeBilinEnvelope1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             xcoef,              /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             ycoef,              /**< y coefficient of linear inequality */
   SCIP_Real             constant,           /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   );

/** computes coefficients of linearization of a bilinear term in a reference point when given two linear inequality
 *  involving only the variables of the bilinear term
 *
 *  @note the formulas are extracted from "Convex envelopes of bivariate functions through the solution of KKT systems"
 *        by Marco Locatelli
 *
 */
SCIP_EXPORT
void SCIPcomputeBilinEnvelope2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             bilincoef,          /**< coefficient of bilinear term */
   SCIP_Real             lbx,                /**< lower bound on first variable */
   SCIP_Real             ubx,                /**< upper bound on first variable */
   SCIP_Real             refpointx,          /**< reference point for first variable */
   SCIP_Real             lby,                /**< lower bound on second variable */
   SCIP_Real             uby,                /**< upper bound on second variable */
   SCIP_Real             refpointy,          /**< reference point for second variable */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_Real             alpha1,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             beta1,              /**< y coefficient of linear inequality */
   SCIP_Real             gamma1,             /**< constant of linear inequality */
   SCIP_Real             alpha2,             /**< x coefficient of linear inequality; must be in {-1,0,1} */
   SCIP_Real             beta2,              /**< y coefficient of linear inequality */
   SCIP_Real             gamma2,             /**< constant of linear inequality */
   SCIP_Real* RESTRICT   lincoefx,           /**< buffer to store coefficient of first  variable in linearization */
   SCIP_Real* RESTRICT   lincoefy,           /**< buffer to store coefficient of second variable in linearization */
   SCIP_Real* RESTRICT   linconstant,        /**< buffer to store constant of linearization */
   SCIP_Bool* RESTRICT   success             /**< buffer to store whether linearization was successful */
   );

/** creates an NLP relaxation and stores it in a given NLPI problem; the function computes for each variable which the
 *  number of non-linearly occurrences and stores it in the nlscore array
 *
 *  @note the first row corresponds always to the cutoff row (even if cutoffbound is SCIPinfinity(scip))
 **/
SCIP_EXPORT
SCIP_RETCODE SCIPcreateNlpiProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLROW**          nlrows,             /**< nonlinear rows */
   int                   nnlrows,            /**< total number of nonlinear rows */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< empty nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   SCIP_HASHMAP*         nlrow2idx,          /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem, can be NULL */
   SCIP_Real*            nlscore,            /**< array to store the score of each nonlinear variable (NULL if not
                                              *   needed) */
   SCIP_Real             cutoffbound,        /**< cutoff bound */
   SCIP_Bool             setobj,             /**< should the objective function be set? */
   SCIP_Bool             onlyconvex          /**< filter only for convex constraints */
   );

/** updates bounds of each variable and the cutoff row in the nlpiproblem */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateNlpiProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_HASHMAP*         var2nlpiidx,        /**< mapping between variables and nlpi indices */
   SCIP_VAR**            nlpivars,           /**< array containing all variables of the nlpi */
   int                   nlpinvars,          /**< total number of nlpi variables */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   );

/** adds linear rows to the NLP relaxation */
SCIP_EXPORT
SCIP_RETCODE SCIPaddNlpiProbRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< nlpi problem */
   SCIP_HASHMAP*         var2idx,            /**< empty hash map to store mapping between variables and indices in nlpi
                                              *   problem */
   SCIP_ROW**            rows,               /**< rows to add */
   int                   nrows               /**< total number of rows to add */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
