/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_branch.h,v 1.3 2004/02/04 17:27:43 bzfpfend Exp $"

/**@file   struct_branch.h
 * @brief  datastructures for branching rules and branching candidate storage
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_BRANCH_H__
#define __STRUCT_BRANCH_H__


#include "def.h"
#include "type_var.h"
#include "type_branch.h"


/** branching candidate storage */
struct BranchCand
{
   VAR**            lpcands;            /**< candidates for branching on LP solution (fractional integer variables) */
   Real*            lpcandssol;         /**< solution values of LP candidates */
   Real*            lpcandsfrac;        /**< fractionalities of LP candidates */
   VAR**            pseudocands;        /**< candidates for branching on pseudo solution (non-fixed integer variables) */
   int              lpcandssize;        /**< number of available slots in lpcands array */
   int              nlpcands;           /**< number of candidates for branching on LP solution */
   int              pseudocandssize;    /**< number of available slots in pseudocands array */
   int              npseudocands;       /**< number of candidates for branching on pseudo solution */
   int              npseudobins;        /**< number of binary candidates for branching on pseudo solution */
   int              npseudoints;        /**< number of integer candidates for branching on pseudo solution */
   int              npseudoimpls;       /**< number of implicit integer candidates for branching on pseudo solution */
   int              validlpcandslp;     /**< lp number for which lpcands are valid */
};

/** branching rule */
struct Branchrule
{
   char*            name;               /**< name of branching rule */
   char*            desc;               /**< description of branching rule */
   int              priority;           /**< priority of the branching rule */
   DECL_BRANCHFREE  ((*branchfree));    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit));    /**< initialize branching rule */
   DECL_BRANCHEXIT  ((*branchexit));    /**< deinitialize branching rule */
   DECL_BRANCHEXECLP((*branchexeclp));  /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXECPS((*branchexecps));  /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata;     /**< branching rule data */
   Bool             initialized;        /**< is branching rule initialized? */
};


#endif
