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



#include "scip.h"


/*
 * lincons methods
 */

/** creates a linear constraint object */
extern
RETCODE SCIPlinconsCreate(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS**        lincons,            /**< pointer to store the linear constraint */
   const char*      name,               /**< name of linear constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            var,                /**< array with variables of constraint entries */
   Real*            val,                /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             model,              /**< does linear constraint is a model constraint? */
   Bool             modifiable          /**< is constraint modifiable during node processing (sbj. to column generation)? */
   );

/** creates a linear constraint object from an LP row */
extern
RETCODE SCIPlinconsCreateLPRow(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS**        lincons,            /**< pointer to store the linear constraint */
   ROW*             row                 /**< LP row to create linear constraint from */
   );

/** frees a linear constraint object */
extern
RETCODE SCIPlinconsFree(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS**        lincons             /**< pointer to store the linear constraint */
   );

/** adds coefficient in linear constraint */
extern
RETCODE SCIPlinconsAddCoef(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   VAR*             var,                /**< variable of constraint entry */
   Real             val                 /**< coefficient of constraint entry */
   );

/** forbids roundings of variables in linear constraint that may violate the constraint */
extern
RETCODE SCIPlinconsForbidRounding(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons             /**< linear constraint */
   );

/** allows roundings of variables in linear constraint that may violate the constraint */
extern
RETCODE SCIPlinconsAllowRounding(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons             /**< linear constraint */
   );

/** gets activity bounds for constraint */
extern
RETCODE SCIPlinconsGetActivityBounds(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Real*            minactivity,        /**< pointer to store the minimal activity */
   Real*            maxactivity         /**< pointer to store the maximal activity */
   );

/** invalidates activity bounds, such that they are recalculated in next get */
extern
RETCODE SCIPlinconsInvalidActivityBounds(
   LINCONS*         lincons             /**< linear constraint */
   );

/** gets activity bounds for constraint after setting variable to zero */
extern
RETCODE SCIPlinconsGetActivityResiduals(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   VAR*             var,                /**< variable to calculate activity residual for */
   Real             val,                /**< coefficient value of variable in linear constraint */
   Real*            minresactivity,     /**< pointer to store the minimal residual activity */
   Real*            maxresactivity      /**< pointer to store the maximal residual activity */
   );

/** gets left hand side of linear constraint */
extern
RETCODE SCIPlinconsGetLhs(
   LINCONS*         lincons,            /**< linear constraint */
   Real*            lhs                 /**< pointer to store left hand side */
   );

/** gets right hand side of linear constraint */
extern
RETCODE SCIPlinconsGetRhs(
   LINCONS*         lincons,            /**< linear constraint */
   Real*            rhs                 /**< pointer to store right hand side */
   );

/** sets left hand side of linear constraint */
extern
RETCODE SCIPlinconsChgLhs(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Real             lhs                 /**< new left hand side */
   );

/** sets right hand side of linear constraint */
extern
RETCODE SCIPlinconsChgRhs(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Real             rhs                 /**< new right hand side */
   );

/** prints linear constraint to file stream */
extern
void SCIPlinconsPrint(
   LINCONS*         lincons,            /**< linear constraint */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** tightens variable's bounds due to activity bounds */
extern
RETCODE SCIPlinconsTightenBounds(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   RESULT*          result              /**< pointer to store SCIP_CUTOFF, if node is infeasible */
   );

/** separates linear constraint: adds linear constraint as cut, if violated by current LP solution */
extern
RETCODE SCIPlinconsSeparate(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   Bool*            violated            /**< pointer to store information, if linear constraint was violated */
   );

/** checks linear constraint for feasibility of given solution or pseudo solution */
extern
RETCODE SCIPlinconsCheck(
   SCIP*            scip,               /**< SCIP data structure */
   LINCONS*         lincons,            /**< linear constraint */
   SOL*             sol,                /**< solution to be checked, NULL to check pseudo solution */
   Bool             chcklprows,         /**< has linear constraint to be checked, if it is already in current LP? */
   Bool*            violated            /**< pointer to store information, if linear constraint is violated */
   );




/*
 * constraint specific interface methods
 */

/** creates the handler for linear constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConsHdlrLinear(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a linear constraint */
extern
RETCODE SCIPcreateConsLinear(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            var,                /**< array with variables of constraint entries */
   Real*            val,                /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of row */
   Real             rhs,                /**< right hand side of row */
   Bool             model,              /**< is constraint necessary for feasibility? */
   Bool             modifiable          /**< is row modifiable during node processing (subject to column generation)? */
   );

/** creates and captures a linear constraint from an LP row, captures the row */
extern
RETCODE SCIPcreateConsLinearLPRow(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   ROW*             row                 /**< LP row */
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


#endif
