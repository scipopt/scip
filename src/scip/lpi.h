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

/**@file   lpi.h
 * @brief  interface methods for specific LP solvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LPI_H__
#define __LPI_H__


typedef struct LPi LPI;                 /**< solver dependent LP interface */
typedef struct LPiState LPISTATE;       /**< complete LP state (i.e. basis information, dual norms) */

/** objective sense */
enum ObjSen
{
   SCIP_OBJSEN_MAXIMIZE = -1,           /**< maximize objective function */
   SCIP_OBJSEN_MINIMIZE = +1            /**< minimize objective function */
};
typedef enum ObjSen OBJSEN;

/** LP solver parameters */
enum LPParam
{
   SCIP_LPPAR_FROMSCRATCH =  0,         /**< solver should start from scratch at next call */
   SCIP_LPPAR_FASTMIP     =  1,         /**< fast mip setting of LP solver */
   SCIP_LPPAR_LPIT1       =  2,         /**< number of simplex iterations in phase 1 */
   SCIP_LPPAR_LPIT2       =  3,         /**< number of simplex iterations in phase 2 */
   SCIP_LPPAR_FEASTOL     =  4,         /**< feasibility tolerance */
   SCIP_LPPAR_LOBJLIM     =  5,         /**< lower objective limit */
   SCIP_LPPAR_UOBJLIM     =  6,         /**< upper objective limit */
   SCIP_LPPAR_LPITLIM     =  7,         /**< LP iteration limit */
   SCIP_LPPAR_LPTILIM     =  8,         /**< LP time limit */
   SCIP_LPPAR_PRICING     =  9,         /**< pricing strategy */
   SCIP_LPPAR_LPINFO      = 10          /**< should LP solver output information to the screen? */
};
typedef enum LPParam LPPARAM;

/** LP pricing strategy */
enum Pricing
{
   SCIP_PRICING_FULL        = 0,        /**< full pricing */
   SCIP_PRICING_STEEP       = 1,        /**< steepest edge pricing */
   SCIP_PRICING_STEEPQSTART = 2         /**< steepest edge pricing without initial dual norms */
};
typedef enum Pricing PRICING;



#include "def.h"
#include "mem.h"
#include "retcode.h"
#include "set.h"



/*
 * LP interface methods
 */

/** gets name and version of LP solver */
extern
const char* SCIPlpiGetSolverName(
   void
   );

/** creates an LP problem object */
extern 
RETCODE SCIPlpiCreate(
   LPI**            lpi,                /**< pointer to an LP interface structure */
   const char*      name                /**< problem name */
   );

/** deletes an LP problem object */
extern
RETCODE SCIPlpiFree(
   LPI**            lpi                 /**< pointer to an LP interface structure */
   );

/** copies data into LP problem object */
extern 
RETCODE SCIPlpiCopyData(
   LPI*             lpi,                /**< LP interface structure */
   int              ncol,               /**< number of columns */
   int              nrow,               /**< number of rows */
   OBJSEN           objsen,             /**< objective sense */
   const Real*      obj,                /**< objective function vector */
   const Real*      rhs,                /**< right hand side vector */
   const char*      sen,                /**< row sense vector */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       cnt,                /**< number of nonzeros for each column */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const Real*      val,                /**< values of constraint matrix entries */
   const Real*      lb,                 /**< lower bound vector */
   const Real*      ub,                 /**< upper bound vector */
   const char**     cname,              /**< column names */
   const char**     rname               /**< row names */
   );

/** adds columns to the LP */
extern 
RETCODE SCIPlpiAddCols(
   LPI*             lpi,                /**< LP interface structure */
   int              ncol,               /**< number of columns to be added */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const Real*      obj,                /**< objective function values of new columns */
   const Real*      lb,                 /**< lower bounds of new columns */
   const Real*      ub,                 /**< upper bounds of new columns */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const Real*      val,                /**< values of constraint matrix entries */
   char**           name,               /**< column names */
   Real             infinity            /**< value used as infinity */
   );

/** deletes columns from LP */
extern 
RETCODE SCIPlpiDelColset(
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of columns
                                         *   input:  1 if column should be deleted, 0 if not
                                         *   output: new position of column, -1 if column was deleted */
   );

/** deletes all columns in the given range from LP */
extern
RETCODE SCIPlpiDelCols(
   LPI*             lpi,                /**< LP interface structure */
   int              firstcol,           /**< first column to be deleted */
   int              lastcol             /**< last column to be deleted */
   );

/** adds rows to the LP */
extern 
RETCODE SCIPlpiAddRows(
   LPI*             lpi,                /**< LP interface structure */
   int              nrow,               /**< number of rows to be added */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const Real*      lhs,                /**< left hand sides of new rows */
   const Real*      rhs,                /**< right hand sides of new rows */
   const int*       beg,                /**< start index of each row in ind- and val-array */
   const int*       ind,                /**< column indices of constraint matrix entries */
   const Real*      val,                /**< values of constraint matrix entries */
   char**           name,               /**< row names */
   Real             infinity            /**< value used as infinity */
   );

/** deletes rows from LP */
extern 
RETCODE SCIPlpiDelRowset(
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of rows
                                         *   input:  1 if row should be deleted, 0 if not
                                         *   output: new position of row, -1 if row was deleted */
   );

/** deletes all rows in the given range from LP */
extern
RETCODE SCIPlpiDelRows(
   LPI*             lpi,                /**< LP interface structure */
   int              firstrow,           /**< first row to be deleted */
   int              lastrow             /**< last row to be deleted */
   );

/** get dense row of inverse basis matrix (A_B)^-1 */
extern 
RETCODE SCIPlpiGetBinvRow(
   LPI*             lpi,                /**< LP interface structure */
   int              i,                  /**< row number */
   Real*            val                 /**< vector to return coefficients */
   );

/** get dense row of inverse basis matrix times constraint matrix (A_B)^-1 * A */
extern 
RETCODE SCIPlpiGetBinvARow(
   LPI*             lpi,                /**< LP interface structure */
   int              i,                  /**< row number */
   const Real*      binv,               /**< dense row vector of row in (A_B)^-1 from prior call to SCIPgetrowBinv() */
   Real*            val                 /**< vector to return coefficients */
   );

/** changes lower and upper bounds of columns */
extern 
RETCODE SCIPlpiChgBounds(
   LPI*             lpi,                /**< LP interface structure */
   int              n,                  /**< number of columns to change bounds for */
   const int*       ind,                /**< column indices */
   const Real*      lb,                 /**< values for the new lower bounds */
   const Real*      ub,                 /**< values for the new upper bounds */
   Real             infinity            /**< value used as infinity */
   );

/** changes left and right hand sides of rows */
extern 
RETCODE SCIPlpiChgSides(
   LPI*             lpi,                /**< LP interface structure */
   int              n,                  /**< number of rows to change sides for */
   const int*       ind,                /**< row indices */
   const Real*      lhs,                /**< new values for left hand sides */
   const Real*      rhs,                /**< new values for right hand sides */
   Real             infinity            /**< value used as infinity */
   );

/** changes the objective sense */
extern 
RETCODE SCIPlpiChgObjsen(
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen              /**< new objective sense */
   );

/** returns the indices of the basic columns and rows */
extern 
RETCODE SCIPlpiGetBind(
   LPI*             lpi,                /**< LP interface structure */
   int*             bind                /**< basic column n gives value n, basic row m gives value -1-m */
   );

/** gets integer parameter of LP */
extern 
RETCODE SCIPlpiGetIntpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int*             ival                /**< buffer to store the parameter value */
   );

/** sets integer parameter of LP */
extern 
RETCODE SCIPlpiSetIntpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int              ival                /**< parameter value */
   );

/** gets floating point parameter of LP */
extern 
RETCODE SCIPlpiGetRealpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   Real*            dval                /**< buffer to store the parameter value */
   );

/** sets floating point parameter of LP */
extern 
RETCODE SCIPlpiSetRealpar(
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   Real             dval                /**< parameter value */
   );

/** gets objective value of solution */
extern
RETCODE SCIPlpiGetObjval(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval              /**< stores the objective value */
   );

/** gets primal and dual solution vectors */
extern 
RETCODE SCIPlpiGetSol(
   LPI*             lpi,                /**< LP interface structure */
   Real*            objval,             /**< stores the objective value */
   Real*            primsol,            /**< primal solution vector */
   Real*            dualsol,            /**< dual solution vector */
   Real*            activity,           /**< row activity vector */
   Real*            redcost             /**< reduced cost vector */
   );

/** gets primal ray for unbounded LPs */
extern 
RETCODE SCIPlpiGetPrimalRay(
   LPI*             lpi,                /**< LP interface structure */
   Real*            ray                 /**< primal ray */
   );

/** gets dual farkas proof for infeasibility */
extern
RETCODE SCIPlpiGetDualfarkas(
   LPI*             lpi,                /**< LP interface structure */
   Real*            dualfarkas          /**< dual farkas row multipliers */
   );

/** performs strong branching iterations on all candidates */
extern 
RETCODE SCIPlpiStrongbranch(
   LPI*             lpi,                /**< LP interface structure */
   const int*       cand,               /**< candidate list */
   int              ncand,              /**< size of candidate list */
   int              itlim,              /**< iteration limit for strong branchings */
   Real*            down,               /**< stores dual bound after branching candidate down */
   Real*            up                  /**< stores dual bound after branching candidate up */
   );

/** calls primal simplex to solve the LP */
extern 
RETCODE SCIPlpiSolvePrimal(
   LPI*             lpi                 /**< LP interface structure */
   );

/** calls dual simplex to solve the LP */
extern 
RETCODE SCIPlpiSolveDual(
   LPI*             lpi                 /**< LP interface structure */
   );

/** gets information about primal and dual feasibility of the LP basis */
extern 
RETCODE SCIPlpiGetBasisFeasibility(
   LPI*             lpi,                /**< LP interface structure */
   Bool*            primalfeasible,     /**< stores primal feasibility status */
   Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff LP is primal unbounded */
extern 
Bool SCIPlpiIsPrimalUnbounded(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP is primal infeasible */
extern 
Bool SCIPlpiIsPrimalInfeasible(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff LP was solved to optimality */
extern 
Bool SCIPlpiIsOptimal(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff actual LP basis is stable */
extern 
Bool SCIPlpiIsStable(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the objective limit was reached */
extern 
Bool SCIPlpiIsObjlimExc(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
extern 
Bool SCIPlpiIsIterlimExc(
   LPI*             lpi                 /**< LP interface structure */
   );

/** returns TRUE iff the time limit was reached */
extern 
Bool SCIPlpiIsTimelimExc(
   LPI*             lpi                 /**< LP interface structure */
   );

/** stores LP state (like basis information) into lpistate object */
extern
RETCODE SCIPlpiGetState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** loads LP state (like basis information) into solver */
extern
RETCODE SCIPlpiSetState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LPISTATE*        lpistate            /**< LP state information (like basis information) */
   );

/** frees LP state information */
extern
RETCODE SCIPlpiFreeState(
   LPI*             lpi,                /**< LP interface structure */
   MEMHDR*          memhdr,             /**< block memory */
   LPISTATE**       lpistate            /**< pointer to LP state information (like basis information) */
   );

/** writes LP state (like basis information) to a file */
extern 
RETCODE SCIPlpiWriteState(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );

/** writes LP to a file */
extern 
RETCODE SCIPlpiWriteLP(
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );


#endif
