/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_quadratic.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for quadratic constraints
 * @author Stefan Vigerske
 * 
 * This constraint handler handles constraints of the form
 * \f[
 *   \ell \leq \sum_{i,j=1}^n a_{i,j} x_ix_j + \sum_{i=1}^n b_i x_i \leq u  
 * \f]
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_QUADRATIC_H__
#define __SCIP_CONS_QUADRATIC_H__

#include "scip/scip.h"
#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif
   
typedef struct SCIP_QuadConsUpgrade SCIP_QUADCONSUPGRADE; /**< quadratic constraint update method */

/** upgrading method for quadratic constraints into more specific constraints
 * 
 * the method might also upgrade only one side of a quadratic constraint
 * if both sides are upgraded into the same constraint, then upgdconslhs and upgdconsrhs should be set to the same pointer 
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cons            : the quadratic constraint to upgrade
 *  - nbinlin         : number of binary variables in linear part
 *  - nbinquad        : number of binary variables in quadratic part
 *  - nintlin         : number of integer variables in linear part
 *  - nintquad        : number of integer variables in quadratic part
 *  - nimpllin        : number of implicit integer variables in linear part
 *  - nimplquad       : number of implicit integer variables in quadratic part
 *  - ncontlin        : number of continuous variables in linear part
 *  - ncontquad       : number of continuous variables in quadratic part
 *  - integral        : TRUE iff constraints activity value is always integral
 *  - upgdconslhs     : pointer to store the upgrade for the quadratic constraint w.r.t. the left  hand side 
 *  - upgdconsrhs     : pointer to store the upgrade for the quadratic constraint w.r.t. the right hand side
 */
#define SCIP_DECL_QUADCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, \
   int nbinlin, int nbinquad, int nintlin, int nintquad, int nimpllin, int nimplquad, int ncontlin, int ncontquad, \
   SCIP_Bool integral, SCIP_CONS** upgdconslhs, SCIP_CONS** upgdconsrhs)

/** creates the handler for quadratic constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrQuadratic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes a quadratic constraint update method into the quadratic constraint handler */
extern
SCIP_RETCODE SCIPincludeQuadconsUpgrade(
   SCIP*                   scip,               /**< SCIP data structure */
   SCIP_DECL_QUADCONSUPGD((*quadconsupgd)),    /**< method to call for upgrading quadratic constraint */
   int                     priority,           /**< priority of upgrading method */
   const char*             conshdlrname        /**< name of the constraint handler */
   );

/** Creates and captures a quadratic constraint.
 * 
 *  The constraint should be given in the form
 *  \f[
 *  \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m a_j y_jz_j \leq u,
 *  \f]
 *  where \f$x_i = y_j = z_k\f$ is possible.
 */
extern
SCIP_RETCODE SCIPcreateConsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< variables in linear part (x_i) or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< coefficients of variables in linear part (b_i) or NULL if nlinvars == 0 */
   int                   nquadterms,         /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms (y_j) or NULL if nquadterms == 0 */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms (z_j) or NULL if nquadterms == 0 */
   SCIP_Real*            quadcoeffs,         /**< array with coefficients of quadratic terms (a_j) or NULL if nquadterms == 0 */
   SCIP_Real             lhs,                /**< left hand side of quadratic equation (l) */
   SCIP_Real             rhs,                /**< right hand side of quadratic equation (u) */
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

/** Creates and captures a quadratic constraint.
 * 
 * The constraint should be given in the form
 * \f[
 * \ell \leq \sum_{i=1}^n b_i x_i + \sum_{j=1}^m (a_j y_j^2 + b_j y_j) + \sum_{k=1}^p c_kv_kw_k \leq u.
 * \f]
 */
extern
SCIP_RETCODE SCIPcreateConsQuadratic2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms (n) */
   SCIP_VAR**            linvars,            /**< array with variables in linear part (x_i) or NULL if nlinvars == 0 */ 
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part (b_i) or NULL if nlinvars == 0 */ 
   int                   nquadvars,          /**< number of quadratic terms (m) */
   SCIP_VAR**            quadvars,           /**< array with variables in quadratic terms (y_j) or NULL if nquadvars == 0 */
   SCIP_Real*            quadlincoefs,       /**< array with linear coefficients of quadratic variables (b_j) or NULL if nquadvars == 0 */
   SCIP_Real*            quadsqrcoefs,       /**< array with coefficients of square terms of quadratic variables (a_j) or NULL if nquadterms == 0 */
   int*                  nadjbilin,          /**< number of bilinear terms where the variable is involved or NULL if nquadterms == 0 */
   int**                 adjbilin,           /**< indices of bilinear terms in which variable is involved or NULL if nquadterms == 0 */
   int                   nbilin,             /**< number of bilinear terms (p) */
   SCIP_VAR**            bilinvars1,         /**< array with first variables in bilinear term (v_k) or NULL if nbilin == 0 */
   SCIP_VAR**            bilinvars2,         /**< array with second variables in bilinear term (w_k) or NULL if nbilin == 0 */
   SCIP_Real*            bilincoefs,         /**< array with coefficients of bilinear term (c_k) or NULL if nbilin == 0 */
   SCIP_Real             lhs,                /**< constraint left hand side (l) */
   SCIP_Real             rhs,                /**< constraint right hand side (u) */
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

/** Gets the number of variables in the linear part of a quadratic constraint.
 */
extern
int SCIPgetNLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the variables in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNLinearVarsQuadratic.
 */
extern
SCIP_VAR** SCIPgetLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the number of variables in the quadratic part of a quadratic constraint.
 */
extern
int SCIPgetNQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the variables in the quadratic part of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic.
 */
extern
SCIP_VAR** SCIPgetQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the coefficients in the linear part of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic.
 */
extern
SCIP_Real* SCIPgetCoefsLinearVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the linear coefficients in the quadratic part of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic.
 */
extern
SCIP_Real* SCIPgetLinearCoefsQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the square coefficients in the quadratic part of a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic.
 */
extern
SCIP_Real* SCIPgetSqrCoefsQuadVarsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the number of bilinear terms in a quadratic constraint.
 */
extern
int SCIPgetNBilinTermsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the first variables in the bilinear terms in a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
extern
SCIP_VAR** SCIPgetBilinVars1Quadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the second variables in the bilinear terms in a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
extern
SCIP_VAR** SCIPgetBilinVars2Quadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the coefficients of the bilinear terms in a quadratic constraint.
 *  Length is given by SCIPgetNBilinTermQuadratic.
 */
extern
SCIP_Real* SCIPgetBilinCoefsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets for each quadratic variable the number of bilinear terms in which the variable is involved in a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic
 */
extern
int* SCIPgetNAdjBilinQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets for each quadratic variable the indices of bilinear terms in which the variable is involved in a quadratic constraint.
 *  Length is given by SCIPgetNQuadVarsQuadratic, length of each entry is given by SCIPgetNAdjBilinQuadratic.
 */
extern
int** SCIPgetAdjBilinQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the left hand side of a quadratic constraint.
 */
extern
SCIP_Real SCIPgetLhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the right hand side of a quadratic constraint.
 */
extern
SCIP_Real SCIPgetRhsQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) convex.
 */
extern
SCIP_Bool SCIPisConvexQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Indicates whether the quadratic function of a quadratic constraint is (known to be) concave.
 */
extern
SCIP_Bool SCIPisConcaveQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   );

/** Gets the violation of a constraint by a solution */
extern
SCIP_RETCODE SCIPgetViolationQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< buffer to store violation of constraint */
   );
   

/** NLPI initialization method of constraint handler
 * 
 *  The constraint handler should create an NLPI representation of the constraints in the provided NLPI.
 */
extern
SCIP_RETCODE SCIPconsInitNlpiQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for quadratic constraints */
   SCIP_NLPI*            nlpi,               /**< NLPI where to add constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_CONS**           conss,              /**< quadratic constraints */
   SCIP_HASHMAP*         var_scip2nlp        /**< mapping from SCIP variables to variable indices in NLPI */
   );

#ifdef __cplusplus
}
#endif

#endif
