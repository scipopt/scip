/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_conflict.h,v 1.22 2006/01/03 12:22:56 bzfpfend Exp $"

/**@file   struct_conflict.h
 * @brief  datastructures for conflict analysis
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONFLICT_H__
#define __SCIP_STRUCT_CONFLICT_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_lpi.h"
#include "scip/type_misc.h"
#include "scip/type_var.h"
#include "scip/type_conflict.h"



/** conflict handler */
struct SCIP_Conflicthdlr
{
   char*                 name;               /**< name of conflict handler */
   char*                 desc;               /**< description of conflict handler */
   SCIP_DECL_CONFLICTFREE((*conflictfree));  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit));  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit));  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol));/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol));/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec));  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata;  /**< conflict handler data */
   int                   priority;           /**< priority of the conflict handler */
   SCIP_Bool             initialized;        /**< is conflict handler initialized? */
};

/** conflict clause */
struct SCIP_Clause
{
   SCIP_VAR**            vars;               /**< literals of the clause, sorted by variable index */
   int                   nvars;              /**< number of literals in the clause */
   int                   validdepth;         /**< depth in the tree where the clause is valid */
   int                   insertdepth;        /**< depth level where constraint should be added */
   int                   conflictdepth;      /**< depth in the tree where the clause yields a conflict */
   int                   repropdepth;        /**< depth at which the clause triggers a deduction */
   int                   score;              /**< score (usefulness) of the clause */
};

/** conflict analysis data structure */
struct SCIP_Conflict
{
   SCIP_Longint          nappliedglbclauses; /**< total number of conflict clauses added globally to the problem */
   SCIP_Longint          nappliedglbliterals;/**< total number of literals in globally applied conflict clauses */
   SCIP_Longint          nappliedlocclauses; /**< total number of conflict clauses added locally to the problem */
   SCIP_Longint          nappliedlocliterals;/**< total number of literals in locally applied conflict clauses */
   SCIP_Longint          npropcalls;         /**< number of calls to propagation conflict analysis */
   SCIP_Longint          npropsuccess;       /**< number of calls yielding at least one conflict clause */
   SCIP_Longint          npropconfclauses;   /**< number of valid conflict clauses detected in propagation conflict analysis */
   SCIP_Longint          npropconfliterals;  /**< total number of literals in valid propagation conflict clauses */
   SCIP_Longint          npropreconvclauses; /**< number of reconvergence clauses detected in propagation conflict analysis */
   SCIP_Longint          npropreconvliterals;/**< total number of literals in valid propagation reconvergence clauses */
   SCIP_Longint          nlpcalls;           /**< number of calls to infeasible LP conflict analysis */
   SCIP_Longint          nlpsuccess;         /**< number of calls yielding at least one conflict clause */
   SCIP_Longint          nlpconfclauses;     /**< number of valid conflict clauses detected in infeas LP conflict analysis */
   SCIP_Longint          nlpconfliterals;    /**< total number of literals in valid infeasible LP conflict clauses */
   SCIP_Longint          nlpreconvclauses;   /**< number of reconvergence clauses detected in infeasible LP conflict analysis */
   SCIP_Longint          nlpreconvliterals;  /**< total number of literals in valid infeasible LP reconvergence clauses */
   SCIP_Longint          nlpiterations;      /**< total number of LP iterations used in LP conflict analysis */
   SCIP_Longint          nsbcalls;           /**< number of calls to infeasible strong branching conflict analysis */
   SCIP_Longint          nsbsuccess;         /**< number of calls yielding at least one conflict clause */
   SCIP_Longint          nsbconfclauses;     /**< number of conflict clauses detected in strong branching conflict analysis */
   SCIP_Longint          nsbconfliterals;    /**< total number of literals in valid strong branching conflict clauses */
   SCIP_Longint          nsbreconvclauses;   /**< number of reconvergence clauses detected in strong branch conflict analysis */
   SCIP_Longint          nsbreconvliterals;  /**< total number of literals in valid strong branching reconvergence clauses */
   SCIP_Longint          nsbiterations;      /**< total number of LP iterations used in strong branching conflict analysis */
   SCIP_Longint          npseudocalls;       /**< number of calls to pseudo solution conflict analysis */
   SCIP_Longint          npseudosuccess;     /**< number of calls yielding at least one conflict clause */
   SCIP_Longint          npseudoconfclauses; /**< number of valid conflict clauses detected in pseudo sol conflict analysis */
   SCIP_Longint          npseudoconfliterals;/**< total number of literals in valid pseudo solution conflict clauses */
   SCIP_Longint          npseudoreconvclauses; /**< number of reconvergence clauses detected in pseudo sol conflict analysis */
   SCIP_Longint          npseudoreconvliterals;/**< total number of literals in valid pseudo solution reconvergence clauses */
   SCIP_CLOCK*           propanalyzetime;    /**< time used for propagation conflict analysis */
   SCIP_CLOCK*           lpanalyzetime;      /**< time used for infeasible LP conflict analysis */
   SCIP_CLOCK*           sbanalyzetime;      /**< time used for infeasible LP conflict analysis */
   SCIP_CLOCK*           pseudoanalyzetime;  /**< time used for pseudo solution conflict analysis */
   SCIP_CLAUSE**         clauses;            /**< conflict clauses found at the current node */
   SCIP_PQUEUE*          binbdchgqueue;      /**< unprocessed conflict bound changes on binary variables */
   SCIP_PQUEUE*          nonbinbdchgqueue;   /**< unprocessed conflict bound changes on non-binary variables */
   SCIP_VAR**            conflictvars;       /**< variables resembling the conflict clause */
   int                   conflictvarssize;   /**< size of conflictvars array */
   int                   nconflictvars;      /**< number of variables in the conflict set (used slots of conflictvars array) */
   int                   ntmpconflictvars;   /**< number of additional variables added temporarily to conflict set */
   int                   count;              /**< conflict set counter to label conflict variables with */
   int                   clausessize;        /**< size of clauses array */
   int                   nclauses;           /**< number of available conflict clauses (used slots in clauses array) */
};


#endif
