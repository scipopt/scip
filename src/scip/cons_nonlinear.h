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

/**@file   cons_nonlinear.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for nonlinear constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_NONLINEAR_H__
#define __SCIP_CONS_NONLINEAR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** linear auxiliary expression of the form xy {<=,>=,==} coefs[0]w + coefs[1]x + coefs[2]y + cst */
struct SCIP_ConsExpr_Auxexpr
{
   SCIP_Real             coefs[3];           /**< coefficients in the expression */
   SCIP_Real             cst;                /**< constant */
   SCIP_VAR*             auxvar;             /**< auxiliary variable w in xy {<=,>=,==} auxexpr(w, x, y) */
   SCIP_Bool             underestimate;      /**< whether the auxexpr underestimates the product */
   SCIP_Bool             overestimate;       /**< whether the auxexpr overestimates the product */
};
typedef struct SCIP_ConsExpr_Auxexpr SCIP_CONSEXPR_AUXEXPR;


/** bilinear term structure
 *
 * This can represent a product which
 * - explicitly exists in the problem and is under- and/or overestimated by a single auxiliary variable
 * stored as auxvar in the union (case nauxexprs == 0) or
 * - is involved in bilinear relations implicitly given by linear constraints with binary variables, and
 * is under- and/or overestimated by linear expression(s) stored as auxexprs in the union (case nauxexprs > 0).
 *
 * An explicitly existing product can also be involved in implicit relations, then it will be stored as in
 * the second case.
 */
struct SCIP_ConsExpr_BilinTerm
{
   SCIP_VAR*             x;                  /**< first variable */
   SCIP_VAR*             y;                  /**< second variable */
   union
   {
      SCIP_CONSEXPR_AUXEXPR** exprs;         /**< auxiliary expressions for the implicit product of x and y */
      SCIP_VAR*          var;                /**< auxiliary variable for the explicit product of x and y */
   } aux;
   int                   nauxexprs;          /**< number of auxexprs (0 for products without implicit relations) */
   int                   sauxexprs;          /**< size of the auxexprs array */
   int                   nlockspos;          /**< number of positive expression locks */
   int                   nlocksneg;          /**< number of negative expression locks */
   SCIP_Bool             existing;           /**< does the product exist explicitly in the problem? */
};
typedef struct SCIP_ConsExpr_BilinTerm SCIP_CONSEXPR_BILINTERM; /**< bilinear term structure */


/** storage for a linear row in preparation
 *
 * Uses to assemble data that could eventually make a SCIP_ROW.
 * @note Only one-sided rows are allowed here.
 */
struct SCIP_RowPrep
{
   SCIP_VAR**            vars;               /**< variables */
   SCIP_Real*            coefs;              /**< coefficients of variables */
   int                   nvars;              /**< number of variables (= number of coefficients) */
   int                   varssize;           /**< length of variables array (= lengths of coefficients array) */
   SCIP_Real             side;               /**< side */
   SCIP_SIDETYPE         sidetype;           /**< type of side */
   SCIP_Bool             local;              /**< whether the row is only locally valid (i.e., for the current node) */
   char                  name[SCIP_MAXSTRLEN]; /**< row name */

   SCIP_Bool             recordmodifications;/**< whether to remember variables which coefficients were modified during cleanup */
   SCIP_VAR**            modifiedvars;       /**< variables which coefficient were modified by cleanup */
   int                   nmodifiedvars;      /**< number of variables which coefficient was modified */
   int                   modifiedvarssize;   /**< length of modifiedvars array */
   SCIP_Bool             modifiedside;       /**< whether the side was modified (relaxed) by cleanup */
};
typedef struct SCIP_RowPrep SCIP_ROWPREP;


/** evaluation callback for (vertex-polyhedral) functions used as input for facet computation of its envelopes
 *
 * input:
 * - args the point to be evaluated
 * - nargs the number of arguments of the function (length of array args)
 * - funcdata user-data of function evaluation callback
 *
 * return:
 * - value of function in point x or SCIP_INVALID if could not be evaluated
 */
#define SCIP_DECL_VERTEXPOLYFUN(f) SCIP_Real f (SCIP_Real* args, int nargs, void* funcdata)

/** maximum dimension of vertex-polyhedral function for which we can try to compute a facet of its convex or concave envelope */
#define SCIP_MAXVERTEXPOLYDIM 14

/** upgrading method for expression constraints into more specific constraints
 *
 * the method might upgrade an expression constraint into a set of upgrade constraints
 * the caller provided an array upgdconss to store upgrade constraints
 * the length of upgdconss is given by upgdconsssize
 * if an upgrade is not possible, set *nupgdconss to zero
 * if more than upgdconsssize many constraints shall replace cons, the function
 * should return the required number as negated value in *nupgdconss
 * i.e., if cons should be replaced by 3 constraints, the function should set
 * *nupgdconss to -3 and return with SCIP_OKAY
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cons            : the nonlinear constraint to upgrade
 *  - nvarexprs       : total number of variable expressions in the expression constraint
 *  - nupgdconss      : pointer to store number of constraints that replace this constraint
 *  - upgdconss       : array to store constraints that replace this constraint
 *  - upgdconsssize   : length of the provided upgdconss array
 */
#define SCIP_DECL_EXPRCONSUPGD(x) SCIP_RETCODE x (SCIP* scip, SCIP_CONS* cons, int nvarexprs, \
      int* nupgdconss, SCIP_CONS** upgdconss, int upgdconsssize)

/** creates the handler for nonlinear constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 * */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@name Expression Constraint Handler Methods */
/**@{ */


/** gets tag indicating current local variable bounds */
SCIP_EXPORT
unsigned int SCIPgetConsExprCurBoundsTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** gets the curboundstag at the last time where variable bounds were relaxed */
SCIP_EXPORT
unsigned int SCIPgetConsExprLastBoundRelaxTag(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** returns the hashmap that is internally used to map variables to their corresponding variable expressions */
SCIP_EXPORT
SCIP_HASHMAP* SCIPgetConsExprVarHashmap(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** notifies conshdlr that a variable expression is to be freed
 *
 * the conshdlr will then update its var2expr hashmap
 *
 * @note To be called only by var-exprhdlr.
 * @note Temporary method that will be replaced by ownerdata-free
 */
SCIP_EXPORT
SCIP_RETCODE SCIPnotifyConsExprExprVarFreed(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_EXPR*        varexpr         /**< variable expression to be freed */
   );

/** collects all bilinear terms for a given set of constraints
 *
 * @note This method should only be used for unit tests that depend on SCIPgetConsExprBilinTerms(),
 *       SCIPgetConsExprBilinTerm() or SCIPgetConsExprBilinTermIdx().
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcollectConsExprBilinTerms(
   SCIP*                      scip,           /**< SCIP data structure */
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_CONS**                conss,          /**< expression constraints */
   int                        nconss          /**< total number of expression constraints */
   );

/** returns the total number of bilinear terms that are contained in all expression constraints
 *
 *  @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 */
SCIP_EXPORT
int SCIPgetConsExprNBilinTerms(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** returns all bilinear terms that are contained in all expression constraints
 *
 * @note This method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @note The value of the auxiliary variable of a bilinear term might be NULL, which indicates that the term does not have an auxiliary variable.
 */
SCIP_EXPORT
SCIP_CONSEXPR_BILINTERM* SCIPgetConsExprBilinTerms(
   SCIP_CONSHDLR*             consexprhdlr    /**< expression constraint handler */
   );

/** returns the index of the bilinear term representing the product of the two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns -1 if the variables do not appear bilinearly.
 */
SCIP_EXPORT
int SCIPgetConsExprBilinTermIdx(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y               /**< second variable */
   );

/** returns the bilinear term that representing the product of two given variables
 *
 * @note The method should only be used after auxiliary variables have been created, i.e., after CONSINITLP.
 * @return The method returns NULL if the variables do not appear bilinearly.
 */
SCIP_EXPORT
SCIP_CONSEXPR_BILINTERM* SCIPgetConsExprBilinTerm(
   SCIP_CONSHDLR*             consexprhdlr,   /**< expression constraint handler */
   SCIP_VAR*                  x,              /**< first variable */
   SCIP_VAR*                  y               /**< second variable */
   );

/** evaluates an auxiliary expression for a bilinear term */
SCIP_EXPORT
SCIP_Real SCIPevalConsExprBilinAuxExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             x,                  /**< first variable of the bilinear term */
   SCIP_VAR*             y,                  /**< second variable of the bilinear term */
   SCIP_CONSEXPR_AUXEXPR* auxexpr,           /**< auxiliary expression */
   SCIP_SOL*             sol                 /**< solution at which to evaluate (can be NULL) */
   );

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertBilinearTermExisting(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   int                   nlockspos,          /**< number of positive expression locks */
   int                   nlocksneg           /**< number of negative expression locks */
   );

/** stores the variables of a bilinear term in the data of the constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertBilinearTermImplicit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_VAR*             x,                  /**< first variable */
   SCIP_VAR*             y,                  /**< second variable */
   SCIP_VAR*             auxvar,             /**< auxiliary variable (might be NULL) */
   SCIP_Real             coefx,              /**< coefficient of x in the auxiliary expression */
   SCIP_Real             coefy,              /**< coefficient of y in the auxiliary expression */
   SCIP_Real             coefaux,            /**< coefficient of auxvar in the auxiliary expression */
   SCIP_Real             cst,                /**< constant of the auxiliary expression */
   SCIP_Bool             overestimate        /**< whether the auxiliary expression overestimates the bilinear product */
   );

/** returns the number of enforcements for an expression */
SCIP_EXPORT
int SCIPgetConsExprExprNEnfos(
   SCIP_EXPR*   expr                /**< expression */
   );

/** returns the data for one of the enforcements of an expression */
SCIP_EXPORT
void SCIPgetConsExprExprEnfoData(
   SCIP_EXPR*   expr,                         /**< expression */
   int                   idx,                          /**< position of enforcement in enfos array */
   SCIP_CONSEXPR_NLHDLR** nlhdlr,                      /**< buffer to store nlhldr */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata,      /**< buffer to store nlhdlr data for expression, or NULL */
   SCIP_CONSEXPR_EXPRENFO_METHOD* nlhdlrparticipation, /**< buffer to store methods where nonlinear handler participates, or NULL */
   SCIP_Bool*            sepabelowusesactivity,        /**< buffer to store whether sepabelow uses activity of some expression, or NULL */
   SCIP_Bool*            sepaaboveusesactivity,        /**< buffer to store whether sepaabove uses activity of some expression, or NULL */
   SCIP_Real*            auxvalue                      /**< buffer to store current auxvalue, or NULL */
   );

/** sets the auxiliary value of expression for one of the enforcements of an expression */
SCIP_EXPORT
void SCIPsetConsExprExprEnfoAuxValue(
   SCIP_EXPR*   expr,               /**< expression */
   int                   idx,                /**< position of enforcement in enfos array */
   SCIP_Real             auxvalue            /**< the new value of auxval */
   );


/** creates the handler for expr constraints and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrExpr(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes an expression constraint upgrade method into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprconsUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECL_EXPRCONSUPGD((*exprconsupgd)),  /**< method to call for upgrading expression constraint, or NULL */
   int                   priority,           /**< priority of upgrading method */
   SCIP_Bool             active,             /**< should the upgrading method by active by default? */
   const char*           conshdlrname        /**< name of the constraint handler */
   );

/** creates and captures a expr constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** creates and captures a expr constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** creates and captures a quadratic expression constraint */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs,          /**< array with coefficients of quadratic terms */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   );

/** @} */


/**@addtogroup CONSHDLRS
 *
 * @{
 *
 * @name Nonlinear Constraints
 *
 * @{
 */

/**@name Expression Constraint Methods */
/**@{ */

/** returns the expression of the given expression constraint */
SCIP_EXPORT
SCIP_EXPR* SCIPgetExprConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the root curvature of the given expression constraint
 *
 * @note The curvature information are computed during CONSINITSOL.
 */
SCIP_EXPORT
SCIP_EXPRCURV SCIPgetCurvatureConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns representation of the expression of the given expression constraint as quadratic form, if possible
 *
 * Only sets *quaddata to non-NULL if the whole expression is quadratic (in the non-extended formulation) and non-linear.
 * That is, the expr in each SCIP_QUADEXPR_QUADTERM will be a variable expressions and
 * \ref SCIPgetConsExprExprVarVar() can be used to retrieve the variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetQuadExprConsExpr(
   SCIP*                    scip,               /**< SCIP data structure */
   SCIP_CONS*               cons,               /**< constraint data */
   SCIP_QUADEXPR** quaddata            /**< buffer to store pointer to quaddata, if quadratic; stores NULL, otherwise */
   );

/** gets the expr constraint as a nonlinear row representation. */
SCIP_EXPORT
SCIP_RETCODE SCIPgetNlRowConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< pointer to store nonlinear row */
   );

/** gets the left hand side of an expression constraint */
SCIP_EXPORT
SCIP_Real SCIPgetLhsConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the right hand side of an expression constraint */
SCIP_EXPORT
SCIP_Real SCIPgetRhsConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** adds coef * var to expression constraint
 *
 * @attention This method can only be called in the problem stage.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddLinearTermConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_Real             coef,               /**< coefficient */
   SCIP_VAR*             var                 /**< variable */
   );

/** gets absolute violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * If this value is at most SCIPfeastol(scip), the constraint would be considered feasible.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetAbsViolationConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   );

/** gets scaled violation of expression constraint
 *
 * This function evaluates the constraints in the given solution.
 *
 * The scaling that is applied to the absolute violation of the constraint
 * depends on the setting of parameter constraints/expr/violscale.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetRelViolationConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to check */
   SCIP_Real*            viol                /**< buffer to store computed violation */
   );

/** gives the unique index of an expression constraint
 *
 * Each expression constraint gets an index assigned when it is created.
 * This index never changes and is unique among all expression constraints
 * within the same SCIP instance.
 * Thus, it can be used to sort a set of expression constraints.
 */
SCIP_EXPORT
int SCIPgetConsExprIndex(
   SCIP_CONS*            cons                /**< constraint data */
   );

/** compares two expression constraints by its index
 *
 * Usable as compare operator in array sort functions.
 */
SCIP_EXPORT
int SCIPcompareConsExprIndex(
   void*                 cons1,
   void*                 cons2
   );

/** returns an equivalent linear constraint if possible */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLinearConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_CONS**           lincons             /**< buffer to store linear constraint data */
   );

/** returns a variable that appears linearly that may be decreased without making any other constraint infeasible */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLinvarMayDecreaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   );

/** returns a variable that appears linearly that may be increased without making any other constraint infeasible */
SCIP_EXPORT
SCIP_RETCODE SCIPgetLinvarMayIncreaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_VAR**            var,                /**< pointer to store the variable */
   SCIP_Real*            coef                /**< pointer to store the coefficient */
   );

/** detects nonlinear handlers that can handle the expressions and creates needed auxiliary variables
 *
 *  @note this method is only used for testing purposes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPdetectConsExprNlhdlrs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss              /**< total number of constraints */
   );

/** add the cut and maybe report branchscores */
SCIP_EXPORT
SCIP_RETCODE SCIPprocessConsExprRowprep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_NLHDLR* nlhdlr,             /**< nonlinear handler which provided the estimator */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_ROWPREP*         rowprep,            /**< cut to be added */
   SCIP_Bool             overestimate,       /**< whether the expression needs to be over- or underestimated */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             auxvalue,           /**< current value of expression w.r.t. auxiliary variables as obtained from EVALAUX */
   SCIP_Bool             allowweakcuts,      /**< whether we should only look for "strong" cuts, or anything that separates is fine */
   SCIP_Bool             branchscoresuccess, /**< buffer to store whether the branching score callback of the estimator was successful */
   SCIP_Bool             inenforcement,      /**< whether we are in enforcement, or only in separation */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result */
   );

/** checks whether an expression is quadratic and returns the corresponding coefficients
 *
 * An expression is quadratic if it is either a square (of some expression), a product (of two expressions),
 * or a sum of terms where at least one is a square or a product.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetConsExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_EXPR*   expr,               /**< expression */
   SCIP_QUADEXPR** quaddata         /**< buffer to store pointer to quadratic representation of expression, if it is quadratic, otherwise stores NULL */
   );

/** @} */

/** creates and captures a nonlinear constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

/** creates and captures a nonlinear constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   );

/** @} */

/** @} */

/**@name Nonlinear Handler Methods */
/**@{ */

/** creates the nonlinearity handler and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprNlhdlrBasic(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLR**      nlhdlr,             /**< buffer where to store nonlinear handler */
   const char*                 name,               /**< name of nonlinear handler (must not be NULL) */
   const char*                 desc,               /**< description of nonlinear handler (can be NULL) */
   int                         detectpriority,     /**< detection priority of nonlinear handler */
   int                         enfopriority,       /**< enforcement priority of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLRDETECT((*detect)),  /**< structure detection callback of nonlinear handler */
   SCIP_DECL_CONSEXPR_NLHDLREVALAUX((*evalaux)), /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_CONSEXPR_NLHDLRDATA*   data                /**< data of nonlinear handler (can be NULL) */
   );


/** computes a facet of the convex or concave envelope of a vertex polyhedral function
 *
 * If \f$ f(x) \f$ is vertex-polyhedral, then \f$ g \f$ is a convex underestimator if and only if
 * \f$ g(v^i) \leq f(v^i), \forall i \f$, where \f$ \{ v^i \}_{i = 1}^{2^n} \subseteq \mathbb R^n \f$ are the vertices
 * of the domain of \f$ x \f$, \f$ [\ell,u] \f$. Hence, we can compute a linear underestimator by solving the following
 * LP (we don't necessarily get a facet of the convex envelope, see Technical detail below):
 *
 * \f{align*}{
 *              \max \, & \alpha^T x^* + \beta \\
 *     s.t. \; & \alpha^T v^i + \beta \le f(v^i), \, \forall i = 1, \ldots, 2^n
 * \f}
 *
 * In principle, one would need to update the LP whenever the domain changes. However, \f$ [\ell,u] = T([0, 1]^n) \f$,
 * where \f$ T \f$ is an affine linear invertible transformation given by \f$ T(y)_i = (u_i - \ell_i) y_i + \ell_i \f$.
 * Working with the change of variables \f$ x = T(y) \f$ allows us to keep the constraints of the LP, even if the domain
 * changes. Indeed, after the change of variables, the problem is: find an affine underestimator \f$ g \f$ such that \f$
 * g(T(y)) \le f(T(y)) \f$, for all \f$ y \in [0, 1]^n \f$. Now \f$ f(T(y)) \f$ is componentwise affine, but still
 * satisfies that \f$ g \f$ is a valid underestimator if and only if \f$ g(T(u)) \leq f(T(u)), \forall u \in \{0, 1\}^n
 * \f$. So we now look for \f$ \bar g(y) := g(T(y)) = g(((u_i - \ell_i) y_i + \ell_i)_i) = \bar \alpha^T y + \bar \beta
 * \f$, where \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i \f$ and \f$ \bar \beta = \sum_i \alpha_i \ell_i + \beta \f$. So
 * we find \f$ \bar g \f$ by solving the LP:
 *
 * \f{align*}{
 *              \max \, & \bar \alpha^T T^{-1}(x^*) + \bar \beta \\
 *     s.t. \; & \bar \alpha^T u + \bar \beta \le f(T(u)), \, \forall u \in \{0, 1\}^n
 * \f}
 *
 * and recover \f$ g \f$ by solving \f$ \bar \alpha_i = (u_i - \ell_i) \alpha_i, \bar \beta = \sum_i \alpha_i \ell_i +
 * \beta \f$. Notice that \f$ f(T(u^i)) = f(v^i) \f$ so the right hand side doesn't change after the change of variables.
 *
 * Furthermore, the LP has more constraints than variables, so we solve its dual:
 * \f{align*}{
 *              \min \, & \sum_i \lambda_i f(v^i) \\
 *     s.t. \; & \sum_i \lambda_i u^i = T^{-1}(x^*) \\
 *             & \sum_i \lambda_i = 1 \\
 *             & \forall i, \, \lambda_i \geq 0
 * \f}
 *
 * In case we look for an overestimate, we do exactly the same, but have to maximize in the dual LP instead
 * of minimize.
 *
 * #### Technical and implementation details
 * -# \f$ U \f$ has exponentially many variables, so we only apply this separator for \f$ n \leq 14 \f$.
 * -# If the bounds are not finite, there is no underestimator. Also, \f$ T^{-1}(x^*) \f$ must be in the domain,
 * otherwise the dual is infeasible.
 * -# After a facet is computed, we check whether it is a valid facet (i.e. we check \f$ \alpha^T v + \beta \le f(v) \f$
 *  for every vertex \f$ v \f$). If we find a violation of at most ADJUSTFACETFACTOR * SCIPlpfeastol, then we weaken \f$
 *  \beta \f$ by this amount, otherwise, we discard the cut.
 * -# If a variable is fixed within tolerances, we replace it with its value and compute the facet of the remaining
 * expression. Note that since we are checking the cut for validity, this will never produce wrong result.
 * -# If \f$ x^* \f$ is in the boundary of the domain, then the LP has infinitely many solutions, some of which might
 * have very bad numerical properties. For this reason, we perturb \f$ x^* \f$ to be in the interior of the region.
 * Furthermore, for some interior points, there might also be infinitely many solutions (e.g. for \f$ x y \f$ in \f$
 * [0,1]^2 \f$ any point \f$ (x^*, y^*) \f$ such that \f$ y^* = 1 - x^* \f$ has infinitely many solutions). For this
 * reason, we perturb any given \f$ x^* \f$. The idea is to try to get a facet of the convex/concave envelope. This only
 * happens when the solution has \f$ n + 1 \f$ non zero \f$ \lambda \f$'s (i.e. the primal has a unique solution).
 * -# We need to compute \f$ f(v^i) \f$ for every vertex of \f$ [\ell,u] \f$. A vertex is encoded by a number between 0
 * and \f$ 2^n - 1 \f$, via its binary representation (0 bit is lower bound, 1 bit is upper bound), so we can compute
 * all these values by iterating between 0 and \f$ 2^n - 1 \f$.
 * -# To check that the computed cut is valid we do the following: we use a gray code to loop over the vertices
 * of the box domain w.r.t. unfixed variables in order to evaluate the underestimator. To ensure the validity of the
 * underestimator, we check whether \f$ \alpha v^i + \beta \le f(v^i) \f$ for every vertex \f$ v^i \f$ and adjust
 * \f$ \beta \f$ if the maximal violation is small.
 *
 * @todo the solution is a facet if all variables of the primal have positive reduced costs (i.e. the solution is
 * unique). In the dual, this means that there are \f$ n + 1 \f$ variables with positive value. Can we use this or some
 * other information to handle any of both cases (point in the boundary or point in the intersection of polytopes
 * defining different pieces of the convex envelope)? In the case where the point is in the boundary, can we use that
 * information to maybe solve another to find a facet? How do the polytopes defining the pieces where the convex
 * envelope is linear looks like, i.e, given a point in the interior of a facet of the domain, does the midpoint of the
 * segment joining \f$ x^* \f$ with the center of the domain, always belongs to the interior of one of those polytopes?
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeFacetVertexPolyhedral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_Bool             overestimate,       /**< whether to compute facet of concave (TRUE) or convex (FALSE) envelope */
   SCIP_DECL_VERTEXPOLYFUN((*function)),     /**< pointer to vertex polyhedral function */
   void*                 fundata,            /**< data for function evaluation (can be NULL) */
   SCIP_Real*            xstar,              /**< point to be separated */
   SCIP_Real*            box,                /**< box where to compute facet: should be lb_1, ub_1, lb_2, ub_2... */
   int                   nallvars,           /**< half of the length of box */
   SCIP_Real             targetvalue,        /**< target value: no need to compute facet if value in xstar would be worse than this value */
   SCIP_Bool*            success,            /**< buffer to store whether a facet could be computed successfully */
   SCIP_Real*            facetcoefs,         /**< buffer to store coefficients of facet defining inequality; must be an array of length at least nallvars */
   SCIP_Real*            facetconstant       /**< buffer to store constant part of facet defining inequality */
);


/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeHyperplaneThreePoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             a1,                 /**< first coordinate of a */
   SCIP_Real             a2,                 /**< second coordinate of a */
   SCIP_Real             a3,                 /**< third coordinate of a */
   SCIP_Real             b1,                 /**< first coordinate of b */
   SCIP_Real             b2,                 /**< second coordinate of b */
   SCIP_Real             b3,                 /**< third coordinate of b */
   SCIP_Real             c1,                 /**< first coordinate of c */
   SCIP_Real             c2,                 /**< second coordinate of c */
   SCIP_Real             c3,                 /**< third coordinate of c */
   SCIP_Real*            alpha,              /**< coefficient of first coordinate */
   SCIP_Real*            beta,               /**< coefficient of second coordinate */
   SCIP_Real*            gamma_,             /**< coefficient of third coordinate */
   SCIP_Real*            delta               /**< constant right-hand side */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
