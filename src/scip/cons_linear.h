/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_linear.h
 * @brief  constraint handler for linear constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_LINEAR_H__
#define __CONS_LINEAR_H__


typedef struct LinCons LINCONS;         /**< linear constraint */
typedef struct LinConsUpgrade LINCONSUPGRADE; /**< linear constraint update method */



/** upgrading method for linear constraints into more specific constraints
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cons            : the linear constraint to upgrade
 *  - nvars           : number of variables in the constraint
 *  - vars            : array with constraint variables
 *  - vals            : array with constraint coefficients
 *  - lhs             : left hand side of linear constraint
 *  - rhs             : right hand side of linear constraint
 *  - nposbin         : number of binary variables with positive coefficient
 *  - nnegbin         : number of binary variables with negative coefficient
 *  - nposint         : number of integer variables with positive coefficient
 *  - nnegint         : number of integer variables with negative coefficient
 *  - nposimpl        : number of implicit integer variables with positive coefficient
 *  - nnegimpl        : number of implicit integer variables with negative coefficient
 *  - nposcont        : number of continous variables with positive coefficient
 *  - nnegcont        : number of continous variables with negative coefficient
 *  - ncoeffspone     : number of +1 coefficients
 *  - ncoeffsnone     : number of -1 coefficients
 *  - ncoeffspint     : number of positive integral coefficients other than +1
 *  - ncoeffsnint     : number of negative integral coefficients other than -1
 *  - ncoeffspfrac    : number of positive fractional coefficients
 *  - ncoeffsnfrac    : number of negative fractional coefficients
 *  - poscoeffsum     : sum of all positive coefficients
 *  - negcoeffsum     : sum of all negative coefficients
 *  - integral        : TRUE iff constraints activity value is always integral
 *  - upgdcons        : pointer to store the upgraded constraint
 *  - upgraded        : pointer to store TRUE iff the constraint was upgraded
 *
 *  possible return values for *result:
 *    SCIP_DIDNOTFIND : the linear constraint data was not upgraded to a more specific constraint
 *    SCIP_SUCCESS    : the linear constraint data was upgraded to the more specific constraint stored in *upgdcons
 */
#define DECL_LINCONSUPGD(x) RETCODE x (SCIP* scip, CONS* cons, int nvars, VAR** vars, Real* vals, Real lhs, Real rhs, \
            int nposbin, int nnegbin, int nposint, int nnegint, int nposimpl, int nnegimpl, int nposcont, int nnegcont, \
            int ncoeffspone, int ncoeffsnone, int ncoeffspint, int ncoeffsnint, int ncoeffspfrac, int ncoeffsnfrac, \
            Real poscoeffsum, Real negcoeffsum, Bool integral, CONS** upgdcons)


#include "scip.h"



/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConsHdlrLinear(
   SCIP*            scip                /**< SCIP data structure */
   );

/** includes a linear constraint update method into the linear constraint handler */
extern
RETCODE SCIPincludeLinconsUpgrade(
   SCIP*            scip,               /**< SCIP data structure */
   DECL_LINCONSUPGD((*linconsupgd)),    /**< method to call for upgrading linear constraint */
   int              priority            /**< priority of upgrading method */
   );

/** creates and captures a linear constraint */
extern
RETCODE SCIPcreateConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of constraint */
   Real             rhs,                /**< right hand side of constraint */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is linear constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable during node processing (subject to col generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

/** adds coefficient in linear constraint */
extern
RETCODE SCIPaddCoefConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   );

/** gets left hand side of linear constraint */
extern
RETCODE SCIPgetLhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real*            lhs                 /**< pointer to store left hand side */
   );

/** gets right hand side of linear constraint */
extern
RETCODE SCIPgetRhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real*            rhs                 /**< pointer to store right hand side */
   );

/** changes left hand side of linear constraint */
extern
RETCODE SCIPchgLhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             lhs                 /**< new left hand side */
   );

/** changes right hand side of linear constraint */
extern
RETCODE SCIPchgRhsConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   Real             rhs                 /**< new right hand side */
   );

/** tries to automatically convert a linear constraint into a more specific and more specialized constraint */
extern
RETCODE SCIPupgradeConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< source constraint to try to convert */
   CONS**           upgdcons            /**< pointer to store upgraded constraint, or NULL if not successful */
   );

#endif
