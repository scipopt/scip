/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_set.h,v 1.3 2004/01/07 13:14:15 bzfpfend Exp $"

/**@file   struct_set.h
 * @brief  datastructures for global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_SET_H__
#define __STRUCT_SET_H__


#include "def.h"
#include "message.h"
#include "type_set.h"
#include "type_buffer.h"
#include "type_clock.h"
#include "type_paramset.h"
#include "type_event.h"
#include "type_scip.h"
#include "type_branch.h"
#include "type_conflict.h"
#include "type_cons.h"
#include "type_disp.h"
#include "type_heur.h"
#include "type_nodesel.h"
#include "type_presol.h"
#include "type_pricer.h"
#include "type_reader.h"
#include "type_sepa.h"


/** global SCIP settings */
struct Set
{
   SCIP*            scip;               /**< very ugly: pointer to scip main data structure for callback methods */
   PARAMSET*        paramset;           /**< set of parameters */
   BUFFER*          buffer;             /**< memory buffers for short living temporary objects */
   READER**         readers;            /**< file readers */
   int              nreaders;           /**< number of file readers */
   int              readerssize;        /**< size of readers array */
   PRICER**         pricers;            /**< variable pricers */
   int              npricers;           /**< number of variable pricers */
   int              nactivepricers;     /**< number of variable pricers used in the current problem */
   int              pricerssize;        /**< size of pricers array */
   Bool             pricerssorted;      /**< are the pricers sorted by activity and priority? */
   CONSHDLR**       conshdlrs;          /**< constraint handlers */
   int              nconshdlrs;         /**< number of constraint handlers */
   int              conshdlrssize;      /**< size of conshdlrs array */
   CONFLICTHDLR**   conflicthdlrs;      /**< conflict handlers */
   int              nconflicthdlrs;     /**< number of conflict handlers */
   int              conflicthdlrssize;  /**< size of conflicthdlrs array */
   Bool             conflicthdlrssorted;/**< are the conflict handlers sorted by priority? */
   PRESOL**         presols;            /**< presolvers */
   int              npresols;           /**< number of presolvers */
   int              presolssize;        /**< size of presols array */
   Bool             presolssorted;      /**< are the presolvers sorted by priority? */
   SEPA**           sepas;              /**< separators */
   int              nsepas;             /**< number of separators */
   int              sepassize;          /**< size of sepas array */
   Bool             sepassorted;        /**< are the separators sorted by priority? */
   HEUR**           heurs;              /**< primal heuristics */
   int              nheurs;             /**< number of primal heuristics */
   int              heurssize;          /**< size of heurs array */
   Bool             heurssorted;        /**< are the heuristics sorted by priority? */
   EVENTHDLR**      eventhdlrs;         /**< event handlers */
   int              neventhdlrs;        /**< number of event handlers */
   int              eventhdlrssize;     /**< size of eventhdlrs array */
   NODESEL**        nodesels;           /**< node selectors */
   int              nnodesels;          /**< number of node selectors */
   int              nodeselssize;       /**< size of nodesels array */
   NODESEL*         actnodesel;         /**< currently used node selector, or NULL if invalid */
   BRANCHRULE**     branchrules;        /**< branching rules */
   int              nbranchrules;       /**< number of branching rules */
   int              branchrulessize;    /**< size of branchrules array */
   Bool             branchrulessorted;  /**< are the branching rules sorted by priority? */
   DISP**           disps;              /**< display columns */
   int              ndisps;             /**< number of display columns */
   int              dispssize;          /**< size of disps array */

   VERBLEVEL        verblevel;          /**< verbosity level of output */
   Bool             catchctrlc;         /**< should the CTRL-C interrupt be catched by SCIP? */
   Real             infinity;           /**< values larger than this are considered infinity */
   Real             epsilon;            /**< absolute values smaller than this are considered zero */
   Real             sumepsilon;         /**< absolute values of sums smaller than this are considered zero */
   Real             feastol;            /**< LP feasibility tolerance */
   Real             cutvioleps;         /**< epsilon for deciding if a cut is violated */
   Real             cutviolepsroot;     /**< epsilon for deciding if a cut is violated in the root node */
   Real             historyeps;         /**< minimal distance value to use for branching history updates */
   Real             memgrowfac;         /**< memory growing factor for dynamically allocated arrays */
   int              memgrowinit;        /**< initial size of dynamically allocated arrays */
   Real             treegrowfac;        /**< memory growing factor for tree array */
   int              treegrowinit;       /**< initial size of tree array */
   Real             pathgrowfac;        /**< memory growing factor for path array */
   int              pathgrowinit;       /**< initial size of path array */
   Real             branchscorefac;     /**< branching score factor to weigh downward and upward gain prediction */
   int              dispwidth;          /**< maximal number of characters in a node information line */
   int              dispfreq;           /**< frequency for displaying node information lines */
   int              dispheaderfreq;     /**< frequency for displaying header lines (every n'th node information line) */
   int              maxpresolrounds;    /**< maximal number of presolving rounds (-1: unlimited) */
   Real             presolabortfac;     /**< abort presolve, if l.t. this frac of the problem was changed in last round */
   int              maxpricevars;       /**< maximal number of variables priced in per pricing round */
   int              maxpricevarsroot;   /**< maximal number of priced variables at the root node */
   Real             abortpricevarsfac;  /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
   int              maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int              maxsepacutsroot;    /**< maximal number of separated cuts at the root node */
   Real             maxconfvarsfac;     /**< maximal fraction of binary variables involved in a conflict clause */
   int              minmaxconfvars;     /**< minimal absolute maximum of variables involved in a conflict clause */
   int              colagelimit;        /**< maximum age a column can reach before it is deleted from the LP */
   int              rowagelimit;        /**< maximum age a row can reach before it is deleted from the LP */
   int              cutagelimit;        /**< maximum age a cut can reach before it is deleted from the global cut pool */
   int              consagelimit;       /**< maximum age an unnecessary constraint can reach before it is deleted */
   int              maxsol;             /**< maximal number of solutions to store in the solution storage */
   Longint          nodelimit;          /**< maximal number of nodes to process (-1: no limit) */
   Real             timelimit;          /**< maximal time in seconds to run */
   Real             memlimit;           /**< maximal memory usage in MB */
   Real             gaplimit;           /**< solving stops, if the given gap is reached */
   int              sollimit;           /**< solving stops, if the given number of solutions were found (-1: no limit) */
   Real             memsavefac;         /**< fraction of maximal memory usage resulting in switch to memory saving mode */
   int              lpsolvefreq;        /**< frequency for solving LP at the nodes (-1: never; 0: only root LP) */
   int              lpsolvedepth;       /**< maximal depth for solving LP at the nodes (-1: no depth limit) */
   Bool             fastmip;            /**< should FASTMIP setting of LP solver be used? */
   Bool             scaling;            /**< should scaling of LP solver be used? */
   Bool             cleanupcols;        /**< should new non-basic columns be removed after LP solving? */
   Bool             cleanuprows;        /**< should new basic rows be removed after LP solving? */
   CLOCKTYPE        clocktype;          /**< default clock type to use */
   Bool             clocksenabled;      /**< is timing enabled? */
};


#endif
