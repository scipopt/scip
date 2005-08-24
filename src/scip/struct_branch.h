/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_branch.h,v 1.18 2005/08/24 17:26:59 bzfpfend Exp $"

/**@file   struct_branch.h
 * @brief  datastructures for branching rules and branching candidate storage
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BRANCH_H__
#define __SCIP_STRUCT_BRANCH_H__


#include "scip/def.h"
#include "scip/type_var.h"
#include "scip/type_branch.h"


/** branching candidate storage */
struct SCIP_BranchCand
{
   SCIP_VAR**            lpcands;            /**< candidates for branching on LP solution (fractional integer variables) */
   SCIP_Real*            lpcandssol;         /**< solution values of LP candidates */
   SCIP_Real*            lpcandsfrac;        /**< fractionalities of LP candidates */
   SCIP_VAR**            pseudocands;        /**< candidates for branching on pseudo solution (non-fixed integer variables) */
   int                   lpcandssize;        /**< number of available slots in lpcands array */
   int                   nlpcands;           /**< number of candidates for branching on LP solution */
   int                   npriolpcands;       /**< number of LP candidates with largest branch priority value */
   int                   npriolpbins;        /**< number of binary LP candidates with largest branch priority value */
   int                   lpmaxpriority;      /**< maximal branch priority of all LP candidates */
   int                   pseudocandssize;    /**< number of available slots in pseudocands array */
   int                   npseudocands;       /**< number of candidates for branching on pseudo solution */
   int                   npriopseudocands;   /**< number of pseudo candidates with largest branch priority value */
   int                   npriopseudobins;    /**< number of binary pseudo candidates with largest branch priority value */
   int                   npriopseudoints;    /**< number of integer pseudo candidates with largest branch priority value */
   int                   pseudomaxpriority;  /**< maximal branch priority of all pseudo candidates */
   int                   validlpcandslp;     /**< lp number for which lpcands are valid */
};

/** branching rule */
struct SCIP_Branchrule
{
   SCIP_Real             maxbounddist;       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_Longint          nlpcalls;           /**< number of times, this branching rule was called on an LP solution */
   SCIP_Longint          npseudocalls;       /**< number of times, this branching rule was called on a pseudo solution */
   SCIP_Longint          ncutoffs;           /**< number of cutoffs found so far by this branching rule */
   SCIP_Longint          ncutsfound;         /**< number of cutting planes found so far by this branching rule */
   SCIP_Longint          nconssfound;        /**< number of cutting constraints added so far by this branching rule (not
                                              *   counting constraint additions to child nodes used for branching) */
   SCIP_Longint          ndomredsfound;      /**< number of domain reductions found so far by this branching rule */
   SCIP_Longint          nchildren;          /**< number of children created so far by this branching rule */
   char*                 name;               /**< name of branching rule */
   char*                 desc;               /**< description of branching rule */
   SCIP_DECL_BRANCHFREE  ((*branchfree));    /**< destructor of branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit));    /**< initialize branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit));    /**< deinitialize branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol));/**< solving process initialization method of branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol));/**< solving process deinitialization method of branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp));  /**< branching execution method for fractional LP solutions */
   SCIP_DECL_BRANCHEXECPS((*branchexecps));  /**< branching execution method for not completely fixed pseudo solutions */
   SCIP_BRANCHRULEDATA*  branchruledata;     /**< branching rule data */
   SCIP_CLOCK*           branchclock;        /**< branching rule execution time */
   int                   priority;           /**< priority of the branching rule */
   int                   maxdepth;           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Bool             initialized;        /**< is branching rule initialized? */
};


#endif
