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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_conflict.h,v 1.16 2005/02/14 13:35:51 bzfpfend Exp $"

/**@file   struct_conflict.h
 * @brief  datastructures for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_CONFLICT_H__
#define __STRUCT_CONFLICT_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_lpi.h"
#include "scip/type_misc.h"
#include "scip/type_var.h"
#include "scip/type_conflict.h"



/** conflict handler */
struct Conflicthdlr
{
   char*            name;               /**< name of conflict handler */
   char*            desc;               /**< description of conflict handler */
   DECL_CONFLICTFREE((*conflictfree));  /**< destructor of conflict handler */
   DECL_CONFLICTINIT((*conflictinit));  /**< initialize conflict handler */
   DECL_CONFLICTEXIT((*conflictexit));  /**< deinitialize conflict handler */
   DECL_CONFLICTINITSOL((*conflictinitsol));/**< solving process initialization method of conflict handler */
   DECL_CONFLICTEXITSOL((*conflictexitsol));/**< solving process deinitialization method of conflict handler */
   DECL_CONFLICTEXEC((*conflictexec));  /**< conflict processing method of conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata;  /**< conflict handler data */
   int              priority;           /**< priority of the conflict handler */
   Bool             initialized;        /**< is conflict handler initialized? */
};

/** conflict analysis data structure */
struct Conflict
{
   Longint          npropcalls;         /**< number of calls to propagation conflict analysis */
   Longint          npropconfclauses;   /**< number of valid conflict clauses detected in propagation conflict analysis */
   Longint          npropconfliterals;  /**< total number of literals in valid propagation conflict clauses */
   Longint          npropreconvclauses; /**< number of reconvergence clauses detected in propagation conflict analysis */
   Longint          npropreconvliterals;/**< total number of literals in valid propagation reconvergence clauses */
   Longint          nlpcalls;           /**< number of calls to infeasible LP conflict analysis */
   Longint          nlpconfclauses;     /**< number of valid conflict clauses detected in infeas LP conflict analysis */
   Longint          nlpconfliterals;    /**< total number of literals in valid infeasible LP conflict clauses */
   Longint          nlpreconvclauses;   /**< number of reconvergence clauses detected in infeasible LP conflict analysis */
   Longint          nlpreconvliterals;  /**< total number of literals in valid infeasible LP reconvergence clauses */
   Longint          nlpiterations;      /**< total number of LP iterations used in LP conflict analysis */
   Longint          nsbcalls;           /**< number of calls to infeasible strong branching conflict analysis */
   Longint          nsbconfclauses;     /**< number of conflict clauses detected in strong branching conflict analysis */
   Longint          nsbconfliterals;    /**< total number of literals in valid strong branching conflict clauses */
   Longint          nsbreconvclauses;   /**< number of reconvergence clauses detected in strong branch conflict analysis */
   Longint          nsbreconvliterals;  /**< total number of literals in valid strong branching reconvergence clauses */
   Longint          nsbiterations;      /**< total number of LP iterations used in strong branching conflict analysis */
   Longint          npseudocalls;       /**< number of calls to pseudo solution conflict analysis */
   Longint          npseudoconfclauses; /**< number of valid conflict clauses detected in pseudo sol conflict analysis */
   Longint          npseudoconfliterals;/**< total number of literals in valid pseudo solution conflict clauses */
   Longint          npseudoreconvclauses; /**< number of reconvergence clauses detected in pseudo sol conflict analysis */
   Longint          npseudoreconvliterals;/**< total number of literals in valid pseudo solution reconvergence clauses */
   CLOCK*           propanalyzetime;    /**< time used for propagation conflict analysis */
   CLOCK*           lpanalyzetime;      /**< time used for infeasible LP conflict analysis */
   CLOCK*           sbanalyzetime;      /**< time used for infeasible LP conflict analysis */
   CLOCK*           pseudoanalyzetime;  /**< time used for pseudo solution conflict analysis */
   LPI*             lpi;                /**< LP problem object for the alternative polyhedron */
   PQUEUE*          binbdchgqueue;      /**< unprocessed conflict bound changes on binary variables */
   PQUEUE*          nonbinbdchgqueue;   /**< unprocessed conflict bound changes on non-binary variables */
   VAR**            conflictvars;       /**< variables resembling the conflict clause */
   int*             conflictvardepths;  /**< depth at which the conflict variables were fixed */
   int              conflictvarssize;   /**< size of conflictvars array */
   int              nconflictvars;      /**< number of variables in the conflict set (used slots of conflictvars array) */
   int              ntmpconflictvars;   /**< number of additional variables added temporarily to conflict set */
   int              count;              /**< conflict set counter to label conflict variables with */
};


#endif
