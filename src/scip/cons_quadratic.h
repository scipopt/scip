/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_quadratic.h,v 1.3 2009/07/22 20:04:48 bzfviger Exp $"

/**@file   cons_quadratic.h
 * @brief  constraint handler for quadratic constraints
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_QUADRATIC_H__
#define __SCIP_CONS_QUADRATIC_H__

#include "scip/scip.h"
#ifdef WITH_NLPI
#include "type_nlpi.h"
#endif

/** creates the handler for quadratic constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a quadratic constraint */
extern
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< variables in linear part */
   SCIP_Real*            lincoeff,           /**< coefficients of variables in linear part */
   int                   nquadterm,          /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< index of first variable in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< index of second variable in quadratic terms */
   SCIP_Real*            quadcoeff,          /**< coefficients of quadratic terms */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** creates and captures a quadratic constraint */
extern
SCIP_RETCODE SCIPcreateConsQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n_linvar,           /**< number of linear terms */
   SCIP_VAR**            linvar,             /**< variables in linear part */ 
   SCIP_Real*            lincoeff,           /**< coefficients of variables in linear part */ 
   int                   n_quadvar,          /**< number of quadratic terms */
   SCIP_VAR**            quadvar,            /**< variables in quadratic terms */
   SCIP_Real*            quadlincoeff,       /**< linear coefficients of quadratic variables */
   SCIP_Real*            quadsqrcoeff,       /**< coefficients of square terms of quadratic variables */
   int*                  n_adjbilin,         /**< number of bilinear terms where the variable is involved */
   int**                 adjbilin,           /**< indices of bilinear terms in which variable is involved */
   int                   n_bilin,            /**< number of bilinear terms */
   SCIP_VAR**            bilinvar1,          /**< first variable in bilinear term */
   SCIP_VAR**            bilinvar2,          /**< second variable in bilinear term */
   SCIP_Real*            bilincoeff,         /**< coefficient of bilinear term */
   SCIP_Real             lhs,                /**< constraint  left hand side */
   SCIP_Real             rhs,                /**< constraint right hand side */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint dynamic? */
   SCIP_Bool             removable           /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

/** Gets the number of variables in the linear term of a quadratic constraint.
 */
extern
int SCIPgetNLinearVarsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the variables in the linear part of a quadratic constraint.
 * Length is given by SCIPgetNLinearVarsQuadratic.
 */
extern
SCIP_VAR** SCIPgetLinearVarsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the number of variables in the quadratic term of a quadratic constraint.
 */
extern
int SCIPgetNQuadVarsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the variables in the quadratic part of a quadratic constraint.
 * Length is given by SCIPgetNQuadVarsQuadratic.
 */
extern
SCIP_VAR** SCIPgetQuadVarsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the coefficients in the linear part of a quadratic constraint.
 * Length is given by SCIPgetNLinearVarsQuadratic.
 */
extern
SCIP_Real* SCIPgetCoeffLinearVarsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the linear coefficients in the quadratic part of a quadratic constraint.
 * Length is given by SCIPgetNQuadVarsQuadratic.
 */
extern
SCIP_Real* SCIPgetLinearCoeffQuadVarsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the square coefficients in the quadratic part of a quadratic constraint.
 * Length is given by SCIPgetNQuadVarsQuadratic.
 */
extern
SCIP_Real* SCIPgetSqrCoeffQuadVarsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the number of bilinear terms in a quadratic constraint.
 */
extern
int SCIPgetNBilinTermQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the first variables in the bilinear terms in a quadratic constraint.
 * Length is given by SCIPgetNBilinTermQuadratic.
 */
extern
SCIP_VAR** SCIPgetBilinVar1Quadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the second variables in the bilinear terms in a quadratic constraint.
 * Length is given by SCIPgetNBilinTermQuadratic.
 */
extern
SCIP_VAR** SCIPgetBilinVar2Quadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the coefficients of the bilinear terms in a quadratic constraint.
 * Length is given by SCIPgetNBilinTermQuadratic.
 */
extern
SCIP_Real* SCIPgetBilinCoeffQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the left hand side of a quadratic constraint.
 */
extern
SCIP_Real SCIPgetLhsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

/** Gets the right hand side of a quadratic constraint.
 */
extern
SCIP_Real SCIPgetRhsQuadratic(
   SCIP_CONS*            cons                /**< pointer to hold the created constraint */
   );

#ifdef WITH_NLPI
/** NLPI initialization method of constraint handler
 * 
 * The constraint handler should create an NLPI representation of the constraints in the provided NLPI.
 */
extern
SCIP_RETCODE SCIPconsInitnlpiQuadratic(
   SCIP*           scip,     /**< SCIP data structure */
   SCIP_CONSHDLR*  conshdlr, /**< quadratic constraint handler - it's C wrapper*/
   SCIP_NLPI*      nlpi,     /**< NLPI where to add constraints */
   int             nconss,   /**< number of constraints */
   SCIP_CONS**     conss,    /**< quadratic constraints */
   SCIP_HASHMAP*   var_scip2nlp /**< mapping from SCIP variables to variable indices in NLPI */
   );
#endif

#endif
