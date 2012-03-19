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

/**@file   cons_polynomial.h
 * @brief  constraint handler for polynomial constraints
 * @author Lignfeng Niu, 
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_POLYNOMIAL_H__
#define __SCIP_CONS_POLYNOMIAL_H__


#include "scip/scip.h"

/** constraint data for polynomial constraints */
typedef struct MonomialTag
{
   /**	struct MonomialTag * MonomialName; */
	int nuses;
	int nvars;
	SCIP_VAR ** vars;
	SCIP_Real * power;
	SCIP_Real coefficient;
}Monomial;

typedef struct PolynomialTag
{
   /**	struct PolynomialTag * PolynomialName; */
	int nuses;
	int nMonomials;
	Monomial ** monomials;
}Polynomial;

/** creates the handler for polynomial constraints and includes it in SCIP */
extern
SCIP_DECL_INCLUDEPLUGIN(SCIPincludeConshdlrPolynomial);

/** creates and captures a polynomial constraint */
extern
SCIP_RETCODE SCIPcreateConsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   Polynomial*           polynomial,
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode      /**< should the node always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   );
#endif
