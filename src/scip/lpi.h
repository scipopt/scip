/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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

#include "def.h"
#include "mem.h"
#include "retcode.h"

typedef struct LPi LPI;                 /**< solver dependent LP interface */
typedef struct LPState LPSTATE;         /**< complete LP state (i.e. basis information, dual norms) */


enum ObjSen                             /**< objective sense */
{
   SCIP_OBJSEN_MAXIMIZE = -1,           /**< maximize objective function */
   SCIP_OBJSEN_MINIMIZE = +1            /**< minimize objective function */
};
typedef enum ObjSen OBJSEN;

enum LPParam                            /**< LP solver parameters */
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

enum Pricing                            /**< LP pricing strategy */
{
   SCIP_PRICING_FULL        = 0,        /**< full pricing */
   SCIP_PRICING_STEEP       = 1,        /**< steepest edge pricing */
   SCIP_PRICING_STEEPQSTART = 2         /**< steepest edge pricing without initial dual norms */
};
typedef enum Pricing PRICING;


extern 
RETCODE SCIPopenLP(                     /**< creates an LP problem object */
   LPI**            lpi,                /**< pointer to an LP interface structure */
   const char*      name                /**< name of the LP */
   );

extern
RETCODE SCIPcloseLP(                    /**< deletes an LP problem object */
   LPI**            lpi                 /**< pointer to an LP interface structure */
   );

extern 
RETCODE SCIPcopyLP(                     /**< copies data into LP problem object */
   LPI*             lpi,                /**< LP interface structure */
   int              ncol,               /**< number of columns */
   int              nrow,               /**< number of rows */
   OBJSEN           objsen,             /**< objective sense */
   const double*    obj,                /**< objective function vector */
   const double*    rhs,                /**< right hand side vector */
   const char*      sen,                /**< row sense vector */
   const int*       beg,                /**< start index of each column in ind- and val-array */
   const int*       cnt,                /**< number of nonzeros for each column */
   const int*       ind,                /**< row indices of constraint matrix entries */
   const double*    val,                /**< values of constraint matrix entries */
   const double*    lb,                 /**< lower bound vector */
   const double*    ub,                 /**< upper bound vector */
   const char**     cname,              /**< column names */
   const char**     rname               /**< row names */
   );

extern 
RETCODE SCIPaddLPRows(                  /**< adds rows to the LP */
   LPI*             lpi,                /**< LP interface structure */
   int              nrow,               /**< number of rows to be added */
   int              nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const double*    rhs,                /**< right hand side vector of new rows */
   const char*      sen,                /**< row senses */
   const int*       beg,                /**< start index of each row in ind- and val-array */
   const int*       ind,                /**< column indices of constraint matrix entries */
   const double*    val,                /**< values of constraint matrix entries */
   const char**     name                /**< row names */
   );

extern 
RETCODE SCIPdelLPRows(                  /**< deletes rows from LP */
   LPI*             lpi,                /**< LP interface structure */
   int*             dstat               /**< deletion status of rows
                                         *   input:  neg. value if row should be deleted, non-neg. value if not
                                         *   output: new position of row, -1 if row was deleted */
   );

extern 
RETCODE SCIPgetrowBinv(                 /**< get dense row of inverse basis matrix (A_B)^-1 */
   LPI*             lpi,                /**< LP interface structure */
   int              i,                  /**< row number */
   double*          val                 /**< vector to return coefficients */
   );

extern 
RETCODE SCIPgetrowBinvA(                /**< get dense row of inverse basis matrix times constraint matrix (A_B)^-1 * A */
   LPI*             lpi,                /**< LP interface structure */
   int              i,                  /**< row number */
   const double*    binv,               /**< dense row vector of row in (A_B)^-1 from prior call to SCIPgetrowBinv() */
   double*          val                 /**< vector to return coefficients */
   );

extern 
RETCODE SCIPgetLb(                      /**< gets lower bounds of variables */
   LPI*             lpi,                /**< LP interface structure */
   int              beg,                /**< first variable to get bound for */
   int              end,                /**< last variable to get bound for */
   double*          lb                  /**< vector to store the bounds */
   );

extern 
RETCODE SCIPgetUb(                      /**< gets upper bounds of variables */
   LPI*             lpi,                /**< LP interface structure */
   int              beg,                /**< first variable to get bound for */
   int              end,                /**< last variable to get bound for */
   double*          ub                  /**< vector to store the bounds */
   );

extern 
RETCODE SCIPchgBd(                      /**< changes bounds of the variables in the LP */
   LPI*             lpi,                /**< LP interface structure */
   int              n,                  /**< number of bounds to be changed */
   const int*       ind,                /**< column indices */
   const char*      lu,                 /**< specifies, if 'L'ower or 'U'pper bound should be changed */
   const double*    bd                  /**< values for the new bounds */
   );

extern 
RETCODE SCIPchgRhs(                     /**< changes right hand sides of rows in the LP */
   LPI*             lpi,                /**< LP interface structure */
   int              n,                  /**< number of rows to change */
   const int*       ind,                /**< row indices */
   const double*    rhs                 /**< new values for right hand sides */
   );

extern 
RETCODE SCIPchgObjsen(                  /**< changes the objective sense */
   LPI*             lpi,                /**< LP interface structure */
   OBJSEN           objsen              /**< new objective sense */
   );

extern 
RETCODE SCIPgetBind(                    /**< returns the indices of the basic columns and rows */
   LPI*             lpi,                /**< LP interface structure */
   int*             bind                /**< basic column n gives value n, basic row m gives value -1-m */
   );

extern 
RETCODE SCIPgetLPIntpar(                /**< gets integer parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int*             ival                /**< buffer to store the parameter value */
   );

extern 
RETCODE SCIPsetLPIntpar(                /**< sets integer parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   int              ival                /**< parameter value */
   );

extern 
RETCODE SCIPgetLPDblpar(                /**< gets double parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   double*          dval                /**< buffer to store the parameter value */
   );

extern 
RETCODE SCIPsetLPDblpar(                /**< sets double parameter of LP */
   LPI*             lpi,                /**< LP interface structure */
   LPPARAM          type,               /**< parameter number */
   double           dval                /**< parameter value */
   );

extern 
RETCODE SCIPgetSol(                     /**< gets primal and dual solution vectors */
   LPI*             lpi,                /**< LP interface structure */
   double*          objval,             /**< stores the objective value */
   double*          psol,               /**< primal solution vector */
   double*          pi,                 /**< dual solution vector */
   double*          slck,               /**< slack vector */
   double*          redcost             /**< reduced cost vector */
   );

extern 
RETCODE SCIPstrongbranch(               /**< performs strong branching iterations on all candidates */
   LPI*             lpi,                /**< LP interface structure */
   const double*    psol,               /**< primal LP solution vector */
   int              ncand,              /**< size of candidate list */
   const int*       cand,               /**< candidate list */
   int              itlim,              /**< iteration limit for strong branchings */
   double*          down,               /**< stores dual bound after branching candidate down */
   double*          up                  /**< stores dual bound after branching candidate up */
   );

extern 
RETCODE SCIPoptLPPrimal(                /**< calls primal simplex to solve the LP */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
RETCODE SCIPoptLPDual(                  /**< calls dual simplex to solve the LP */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisPrimalUnbounded(             /**< returns TRUE iff LP is primal unbounded */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisPrimalInfeasible(            /**< returns TRUE iff LP is primal infeasible */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisOptimal(                     /**< returns TRUE iff LP was solved to optimality */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisDualValid(                   /**< returns TRUE iff actual LP solution is dual valid */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisStable(                      /**< returns TRUE iff actual LP basis is stable */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisLPError(                     /**< returns TRUE iff an error occured while solving the LP */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisObjlimExc(                   /**< returns TRUE iff the objective limit was reached */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisIterlimExc(                  /**< returns TRUE iff the iteration limit was reached */
   LPI*             lpi                 /**< LP interface structure */
   );

extern 
Bool SCIPisTimelimExc(                  /**< returns TRUE iff the time limit was reached */
   LPI*             lpi                 /**< LP interface structure */
   );

extern
RETCODE SCIPsaveLPState(                /**< stores LP state (like basis information) */
   MEM*             mem,                /**< block memory buffers */
   LPI*             lpi,                /**< LP interface structure */
   LPSTATE**        lpstate             /**< pointer to LP state information (like basis information) */
   );

extern
RETCODE SCIPloadLPState(                /**< loads LP state (like basis information) into solver */
   MEM*             mem,                /**< block memory buffers */
   LPI*             lpi,                /**< LP interface structure */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   );

extern
void SCIPcaptureLPState(                /**< increases usage counter of LP state */
   LPSTATE*         lpstate             /**< LP state information (like basis information) */
   );
   
extern
void SCIPreleaseLPState(                /**< decreases usage counter of LP state, and frees memory if necessary */
   MEM*             mem,                /**< block memory buffers */
   LPSTATE**        lpstate             /**< LP state information (like basis information) */
   );


extern 
RETCODE SCIPwriteLPState(               /**< writes LP state (like basis information) to a file */
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );

extern 
RETCODE SCIPwriteLP(                    /**< writes LP to a file */
   LPI*             lpi,                /**< LP interface structure */
   const char*      fname               /**< file name */
   );


#endif
