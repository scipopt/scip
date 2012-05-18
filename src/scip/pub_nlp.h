/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_nlp.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for NLP management
 * @author Thorsten Gellermann
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_NLP_H__
#define __SCIP_PUB_NLP_H__

#include <stdio.h>

#include "scip/def.h"
#include "scip/type_message.h"
#include "blockmemshell/memory.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"
#include "nlpi/type_expr.h"
#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@name Expressions and Expression tree methods */
/**@{ */
   
/** returns variables of expression tree */
extern
SCIP_VAR** SCIPexprtreeGetVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
);

/** stores array of variables in expression tree */
extern
SCIP_RETCODE SCIPexprtreeSetVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
);

/** adds variables to the expression tree variables array */
extern
SCIP_RETCODE SCIPexprtreeAddVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
);

/** prints an expression tree using variable names from variables array */
extern
SCIP_RETCODE SCIPexprtreePrintWithNames(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file for printing, or NULL for stdout */
);

/** searches the variables array of an expression tree for a variable and returns its position, or -1 if not found
 * Note that this is an O(n) operation!
 */
extern
int SCIPexprtreeFindVar(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_VAR*             var                 /**< variable to search for */
);

/** removes fixed variables from an expression tree, so that at exit all variables are active */
extern
SCIP_RETCODE SCIPexprtreeRemoveFixedVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Bool*            changed,            /**< buffer to store whether the tree was changed, i.e., whether there was a fixed variable */
   int*                  varpos,             /**< array of length at least tree->nvars to store new indices of previously existing variables in expression tree, or -1 if variable was removed; set to NULL if not of interest */
   int*                  newvarsstart        /**< buffer to store index in tree->vars array where new variables begin, or NULL if not of interest */
);

/**@} */

/**@name Nonlinear row methods */
/**@{ */

/** gets constant */
extern
SCIP_Real SCIPnlrowGetConstant(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets number of variables of linear part */
extern
int SCIPnlrowGetNLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets array with variables of linear part */
extern
SCIP_VAR** SCIPnlrowGetLinearVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets array with coefficients in linear part */
extern
SCIP_Real* SCIPnlrowGetLinearCoefs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets number of quadratic variables in quadratic part */
extern
int SCIPnlrowGetNQuadVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets quadratic variables in quadratic part */
extern
SCIP_VAR** SCIPnlrowGetQuadVars(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gives position of variable in quadvars array of row, or -1 if not found */
extern
int SCIPnlrowSearchQuadVar(
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_VAR*             var                 /**< variable to search for */
   );

/** gets number of quadratic elements in quadratic part */
extern
int SCIPnlrowGetNQuadElems(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets quadratic elements in quadratic part */
extern
SCIP_QUADELEM* SCIPnlrowGetQuadElems(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets array with coefficients in linear part */
extern
void SCIPnlrowGetQuadData(
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int*                  nquadvars,          /**< buffer to store number of variables in quadratic term, or NULL if not of interest */
   SCIP_VAR***           quadvars,           /**< buffer to store pointer to array of variables in quadratic term, or NULL if not of interest */
   int*                  nquadelems,         /**< buffer to store number of entries in quadratic term, or NULL if not of interest */
   SCIP_QUADELEM**       quadelems           /**< buffer to store pointer to array of entries in quadratic term, or NULL if not of interest */
   );

/** gets expression tree */
extern
SCIP_EXPRTREE* SCIPnlrowGetExprtree(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns the left hand side of a nonlinear row */
extern
SCIP_Real SCIPnlrowGetLhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns the right hand side of a nonlinear row */
extern
SCIP_Real SCIPnlrowGetRhs(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns the name of a nonlinear row */
extern
const char* SCIPnlrowGetName(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets position of a nonlinear row in current NLP, or -1 if not in NLP */
extern
int SCIPnlrowGetNLPPos(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** returns TRUE iff row is member of current NLP */
extern
SCIP_Bool SCIPnlrowIsInNLP(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/** gets the dual NLP solution of a nlrow
 * for a ranged constraint, the dual value is positive if the right hand side is active and negative if the left hand side is active
 */
extern
SCIP_Real SCIPnlrowGetDualsol(
   SCIP_NLROW*           nlrow               /**< NLP row */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_PUB_NLP_H__ */
