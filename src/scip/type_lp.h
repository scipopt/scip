/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_lp.h,v 1.5 2005/01/18 09:26:58 bzfpfend Exp $"

/**@file   type_lp.h
 * @brief  type definitions for LP management
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_LP_H__
#define __TYPE_LP_H__


/** solution status after solving LP */
enum LPSolStat
{
   SCIP_LPSOLSTAT_NOTSOLVED    = 0,     /**< LP was not solved, no solution exists */
   SCIP_LPSOLSTAT_OPTIMAL      = 1,     /**< LP was solved to optimality */
   SCIP_LPSOLSTAT_INFEASIBLE   = 2,     /**< LP is primal infeasible */
   SCIP_LPSOLSTAT_UNBOUNDEDRAY = 3,     /**< LP has a primal unbounded ray */
   SCIP_LPSOLSTAT_OBJLIMIT     = 4,     /**< objective limit was reached during optimization */
   SCIP_LPSOLSTAT_ITERLIMIT    = 5,     /**< iteration limit was reached during optimization */
   SCIP_LPSOLSTAT_TIMELIMIT    = 6,     /**< time limit was reached during optimization */
   SCIP_LPSOLSTAT_ERROR        = 7      /**< an error occured during optimization */
};
typedef enum LPSolStat LPSOLSTAT;

/** type of variable bound: lower or upper bound */
enum BoundType
{
   SCIP_BOUNDTYPE_LOWER = 0,            /**< lower bound */
   SCIP_BOUNDTYPE_UPPER = 1             /**< upper bound */
};
typedef enum BoundType BOUNDTYPE;

/** type of row side: left hand or right hand side */
enum SideType
{
   SCIP_SIDETYPE_LEFT  = 0,             /**< left hand side */
   SCIP_SIDETYPE_RIGHT = 1              /**< right hand side */
};
typedef enum SideType SIDETYPE;

typedef struct Col COL;                 /**< column of an LP */
typedef struct Row ROW;                 /**< row of an LP */
typedef struct Lp LP;                   /**< current LP data */


#endif
