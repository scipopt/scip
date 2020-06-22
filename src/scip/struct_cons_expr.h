/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_cons_expr.h
 * @brief  (public) data structures of expression constraints
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 *
 * These are in particular data structures to manage the expressions in cons_expr.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONS_EXPR_H__
#define __SCIP_STRUCT_CONS_EXPR_H__

#include "scip/type_cons_expr.h"
#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** generic data and callback methods of an expression handler */
struct SCIP_ConsExpr_ExprHdlr
{
   char*                   name;          /**< expression handler name */
   char*                   desc;          /**< expression handler description (can be NULL) */
   SCIP_CONSEXPR_EXPRHDLRDATA*   data;    /**< data of handler */
   unsigned int            precedence;    /**< precedence of expression operation relative to other expression (used for printing) */

   SCIP_Longint            nestimatecalls;/**< number of times, the estimation callback were called */
   SCIP_Longint            nintevalcalls; /**< number of times, the interval evaluation callback was called */
   SCIP_Longint            npropcalls;    /**< number of times, the propagation callback was called */
   SCIP_Longint            ncutsfound;    /**< number of cuts added by this expression handler */
   SCIP_Longint            ncutoffs;      /**< number of cutoffs found so far by this expression handler */
   SCIP_Longint            ndomreds;      /**< number of domain reductions found so far by this expression handler */
   SCIP_Longint            nsimplifycalls; /**< number of times, the simplification callback was called */
   SCIP_Longint            nsimplified;   /**< number of times the simplification callback was succesful */
   SCIP_Longint            nbranchscores; /**< number of times, branching scores were added by (or for) this expression handler */

   SCIP_CLOCK*             estimatetime;  /**< time used for estimation */
   SCIP_CLOCK*             proptime;      /**< time used for propagation */
   SCIP_CLOCK*             intevaltime;   /**< time used for interval evaluation */
   SCIP_CLOCK*             simplifytime;  /**< time used for expression simplification */

   SCIP_DECL_CONSEXPR_EXPRCOPYHDLR((*copyhdlr));  /**< handler copy callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRFREEHDLR((*freehdlr));  /**< handler free callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA((*copydata));  /**< data copy callback, or NULL for expressions that have no data */
   SCIP_DECL_CONSEXPR_EXPRFREEDATA((*freedata));  /**< data free callback, or NULL for expressions that have no data or which data does not need to be freed */
   SCIP_DECL_CONSEXPR_EXPRSIMPLIFY((*simplify));  /**< simplify callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRCOMPARE((*compare));    /**< compare callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRPRINT((*print));        /**< print callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRPARSE((*parse));        /**< parse callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPREVAL((*eval));          /**< point evaluation callback (can never be NULL) */
   SCIP_DECL_CONSEXPR_EXPRBWDIFF((*bwdiff));      /**< derivative evaluation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRINTEVAL((*inteval));    /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRINITSEPA((*initsepa));  /**< separation initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPREXITSEPA((*exitsepa));  /**< separation deinitialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRESTIMATE((*estimate));  /**< estimation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRREVERSEPROP((*reverseprop)); /**< reverse propagation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRHASH((*hash));          /**< hash callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRCURVATURE((*curvature)); /**< curvature detection callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRMONOTONICITY((*monotonicity)); /**< monotonicity detection callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRINTEGRALITY((*integrality)); /**< integrality detection callback (can be NULL) */
};

/* expression iteration data */
struct SCIP_ConsExpr_Expr_IterData
{
   SCIP_CONSEXPR_EXPR*      parent;       /**< parent expression in DFS iteration */
   int                      currentchild; /**< child that is currently visited (or will be visited next) by DFS iteration */
   unsigned int             visitedtag;   /**< tag to identify whether an expression has been visited already */
   SCIP_CONSEXPRITERATOR_USERDATA userdata; /**< space for iterator user to store some (temporary) data */
};

/** a node in the expression graph that is handled by the expression constraint handler */
struct SCIP_ConsExpr_Expr
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;      /**< expression type (as pointer to its handler) */
   SCIP_CONSEXPR_EXPRDATA* exprdata;      /**< expression data */

   int                     nchildren;     /**< number of children */
   int                     childrensize;  /**< length of children array */
   SCIP_CONSEXPR_EXPR**    children;      /**< children expressions */

   int                     nuses;         /**< reference counter */
   int                     nlockspos;     /**< positive locks counter */
   int                     nlocksneg ;    /**< negative locks counter */

   /* enforcement of expr == auxvar (or expr <= auxvar, or expr >= auxvar) */
   SCIP_CONSEXPR_EXPRENFO** enfos;        /**< enforcements */
   int                     nenfos;        /**< number of enforcements */
   unsigned int            lastenforced;  /**< last enforcement round where expression was enforced successfully */
   SCIP_Bool               enfoinitialized;/**< whether enforcements have been initialized, i.e., expr ran through DETECT (may still have nenfos=0) */
   unsigned int            nactivityusesprop; /**< number of nonlinear handler whose activity computation (or domain propagation) depends on the activity of the expression */
   unsigned int            nactivityusessepa; /**< number of nonlinear handler whose separation (estimate or enfo) depends on the activity of the expression */
   unsigned int            nauxvaruses;   /**< number of nonlinear handlers whose separation uses an auxvar in the expression */

   /* separation */
   SCIP_VAR*               auxvar;        /**< auxiliary variable used for outer approximation cuts */
   int                     auxfilterpos;  /**< filter position of variable event data for auxiliary variable */

   /* branching */
   SCIP_Real               violscoresum;  /**< sum of violation scores for branching stored for this expression */
   SCIP_Real               violscoremax;  /**< max of violation scores for branching stored for this expression */
   int                     nviolscores;   /**< number of violation scores stored for this expression */
   unsigned int            violscoretag;  /**< tag to decide whether a violation score of an expression needs to be initialized */

   /* point-evaluation */
   unsigned int            evaltag;       /**< tag of point for which the expression has been evaluated last, or 0 */
   SCIP_Real               evalvalue;     /**< value of expression from last evaluation (corresponding to evaltag) */
   SCIP_Real               derivative;    /**< partial derivative of a "root path" w.r.t. this expression
                                            *  (see documentation of Differentiation methods in cons_expr.c) */
   unsigned int            difftag;       /**< when computing partial derivatives of an expression w.r.t. a variable,
                                            *  the tag allows us to decide whether the expression depends on the
                                            *  variable; the tag will be checked in SCIPgetConsExprExprPartialDiff() */

   /* activity */
   SCIP_INTERVAL           activity;      /**< activity of expression with respect to local variable bounds */
   unsigned int            activitytag;   /**< tag of local variable bounds for which activity is valid */
   SCIP_Bool               inqueue;       /**< flag to store whether an expression is in the queue of reverse propagation */

   /* expression iterators data */
   SCIP_CONSEXPR_EXPR_ITERDATA iterdata[SCIP_CONSEXPRITERATOR_MAXNACTIVE];  /**< data for expression iterators */

   /* curvature information */
   SCIP_EXPRCURV           curvature;     /**< curvature of the expression w.r.t. bounds that have been used in the last curvature detection */

   /* monotonicity information of each child */
   SCIP_MONOTONE*          monotonicity;  /**< array containing monotonicity of expression w.r.t. each children */
   int                     monotonicitysize; /**< length of monotonicity array */

   /* integrality information */
   SCIP_Bool               isintegral;     /**< flag to store whether an expression is integral */

   SCIP_CONSEXPR_QUADEXPR* quaddata;       /**< representation of expression as a quadratic, if checked and being quadratic */
   SCIP_Bool               quadchecked;    /**< whether we checked whether the expression is quadratic */
};

typedef struct SCIP_ConsExpr_QuadExprTerm  SCIP_CONSEXPR_QUADEXPRTERM;  /**< a single term associated to a quadratic variable */
typedef struct SCIP_ConsExpr_BilinExprTerm SCIP_CONSEXPR_BILINEXPRTERM; /**< a single bilinear term */

/** data for representation of an expression as quadratic */
struct SCIP_ConsExpr_QuadExpr
{
   SCIP_Real                    constant;        /**< a constant term */

   int                          nlinexprs;       /**< number of expressions that appear linearly */
   SCIP_CONSEXPR_EXPR**         linexprs;        /**< expressions that appear linearly */
   SCIP_Real*                   lincoefs;        /**< coefficients of expressions that appear linearly */

   int                          nquadexprs;      /**< number of expressions in quadratic terms */
   SCIP_CONSEXPR_QUADEXPRTERM*  quadexprterms;   /**< array with quadratic expression terms */

   int                          nbilinexprterms; /**< number of bilinear expressions terms */
   SCIP_CONSEXPR_BILINEXPRTERM* bilinexprterms;  /**< bilinear expression terms array */

   SCIP_Bool                    allexprsarevars; /**< whether all arguments (linexprs, quadexprterms[.].expr) are variable expressions */

   SCIP_EXPRCURV                curvature;       /**< curvature of the quadratic representation of the expression */
   SCIP_Bool                    curvaturechecked;/**< whether curvature has been checked */
};

/** data structure to store a single term associated to a quadratic variable */
struct SCIP_ConsExpr_QuadExprTerm
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< quadratic expression */
   SCIP_Real             lincoef;            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef;            /**< square coefficient of variable */

   int                   nadjbilin;          /**< number of bilinear terms this variable is involved in */
   int                   adjbilinsize;       /**< size of adjacent bilinear terms array */
   int*                  adjbilin;           /**< indices of associated bilinear terms */
};

/** data structure to store a single bilinear term coef * expr1 * expr2  (similar to SCIP_CONSEXPR_QUADELEM)
 * except for temporary reasons, we assume that the index of var1 is smaller than the index of var2
 */
struct SCIP_ConsExpr_BilinExprTerm
{
   SCIP_CONSEXPR_EXPR*   expr1;              /**< first factor of bilinear term */
   SCIP_CONSEXPR_EXPR*   expr2;              /**< second factor of bilinear term */
   SCIP_Real             coef;               /**< coefficient of bilinear term */
   int                   pos2;               /**< position of expr2's quadexprterm in quadexprterms */
};

/** generic data and callback methods of an nonlinear handler */
struct SCIP_ConsExpr_Nlhdlr
{
   char*                         name;             /**< nonlinearity handler name */
   char*                         desc;             /**< nonlinearity handler description (can be NULL) */
   SCIP_CONSEXPR_NLHDLRDATA*     data;             /**< data of handler */
   int                           detectpriority;   /**< detection priority of nonlinearity handler */
   int                           enfopriority;     /**< enforcement priority of nonlinearity handler */
   SCIP_Bool                     enabled;          /**< whether the nonlinear handler should be used */

   SCIP_Longint                  nenfocalls; /**< number of times, the enforcement or estimation callback was called */
   SCIP_Longint                  nintevalcalls; /**< number of times, the interval evaluation callback was called */
   SCIP_Longint                  npropcalls; /**< number of times, the propagation callback was called */
   SCIP_Longint                  nseparated; /**< number of times, the expression handler enforced by separation */
   SCIP_Longint                  ncutoffs;   /**< number of cutoffs found so far by this nonlinear handler */
   SCIP_Longint                  ndomreds;   /**< number of domain reductions found so far by this expression handler */
   SCIP_Longint                  ndetections;/**< number of detect calls in which structure was detected (success returned by detect call) (over all runs) */
   SCIP_Longint                  ndetectionslast;/**< number of detect calls in which structure was detected (success returned by detect call) (in last round) */
   SCIP_Longint                  nbranchscores; /**< number of times, branching scores were added by this nonlinear handler */
   SCIP_Longint                  nreformulates; /**< number of times, an expression has been successfully reformulated by a nonlinear handler */

   SCIP_CLOCK*                   detecttime; /**< time used for detection */
   SCIP_CLOCK*                   enfotime;   /**< time used for enforcement or estimation */
   SCIP_CLOCK*                   proptime;   /**< time used for reverse propagation */
   SCIP_CLOCK*                   intevaltime;/**< time used for interval evaluation */
   SCIP_CLOCK*                   reformulatetime;/**< time used for expression reformulation */

   SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA((*freehdlrdata));  /**< callback to free data of handler (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA((*freeexprdata));  /**< callback to free expression specific data (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR((*copyhdlr));          /**< callback to copy nonlinear handler (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRINIT((*init));                  /**< initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLREXIT((*exit));                  /**< deinitialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRDETECT((*detect));              /**< structure detection callback */
   SCIP_DECL_CONSEXPR_NLHDLREVALAUX((*evalaux));            /**< auxiliary evaluation callback */
   SCIP_DECL_CONSEXPR_NLHDLRINITSEPA((*initsepa));          /**< separation initialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRENFO((*enfo));                  /**< enforcement callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRESTIMATE((*estimate));          /**< estimator callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLREXITSEPA((*exitsepa));          /**< separation deinitialization callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRINTEVAL((*inteval));            /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP((*reverseprop));    /**< reverse propagation callback (can be NULL) */
   SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE((*reformulate));    /**< reformulation callback (can be NULL) */
};

/** enforcement data of an expression */
struct SCIP_ConsExpr_ExprEnfo
{
   SCIP_CONSEXPR_NLHDLR*         nlhdlr;          /**< nonlinear handler */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata;  /**< data of nonlinear handler */
   SCIP_CONSEXPR_EXPRENFO_METHOD nlhdlrparticipation; /**< methods where nonlinear handler participates */
   SCIP_Bool                     issepainit;      /**< was the initsepa callback of nlhdlr called */
   SCIP_Real                     auxvalue;        /**< auxiliary value of expression w.r.t. currently enforced solution */
   SCIP_Bool                     sepausesactivity;/**< whether separation uses activity of some expression */
};

/** expression tree iterator */
struct SCIP_ConsExpr_Iterator
{
   SCIP_CONSHDLR*              consexprhdlr; /**< expr constraint handler */
   BMS_BLKMEM*                 blkmem;       /**< block memory */

   SCIP_Bool                   initialized;  /**< whether the iterator has been initialized, that is, is in use */
   SCIP_CONSEXPRITERATOR_TYPE  itertype;     /**< type of expression iterator */
   SCIP_CONSEXPR_EXPR*         curr;         /**< current expression of the iterator */
   int                         iterindex;    /**< index of iterator data in expressions, or -1 if not using iterator data in expressions */
   unsigned int                visitedtag;   /**< tag to mark and recognize an expression as visited, or 0 if not avoiding multiple visits */

   /* data for rtopological mode */
   SCIP_CONSEXPR_EXPR**        dfsexprs;     /**< DFS stack */
   int*                        dfsnvisited;  /**< number of visited children for each expression in the stack */
   int                         dfsnexprs;    /**< total number of expression in stack */
   int                         dfssize;      /**< size DFS stack */

   /* data for BFS mode */
   SCIP_QUEUE*                 queue;        /**< BFS queue */

   /* data for DFS mode */
   SCIP_CONSEXPRITERATOR_STAGE dfsstage;     /**< current stage */
   unsigned int                stopstages;   /**< stages in which to interrupt iterator */
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_CONS_EXPR_H__ */
