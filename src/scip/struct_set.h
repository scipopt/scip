/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_set.h,v 1.23 2004/07/13 15:03:52 bzfpfend Exp $"

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
   Real             infinity;           /**< values larger than this are considered infinity */
   Real             epsilon;            /**< absolute values smaller than this are considered zero */
   Real             sumepsilon;         /**< absolute values of sums smaller than this are considered zero */
   Real             feastol;            /**< feasibility tolerance for constraints */
   Real             dualfeastol;        /**< feasibility tolerance for reduced costs */
   Real             boundstreps;        /**< minimal improve for strengthening bounds */
   Real             cutvioleps;         /**< epsilon for deciding if a cut is violated */
   Real             cutviolepsroot;     /**< epsilon for deciding if a cut is violated in the root node */
   Real             pseudocosteps;      /**< minimal variable distance value to use for pseudo cost updates */
   Real             pseudocostdelta;    /**< minimal objective distance value to use for pseudo cost updates */
   Real             memgrowfac;         /**< memory growing factor for dynamically allocated arrays */
   Real             treegrowfac;        /**< memory growing factor for tree array */
   Real             pathgrowfac;        /**< memory growing factor for path array */
   Real             branchscorefac;     /**< branching score factor to weigh downward and upward gain prediction */
   Real             presolabortfac;     /**< abort presolve, if l.t. this frac of the problem was changed in last round */
   Real             abortpricevarsfac;  /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
   Real             maxconfvarsfac;     /**< maximal fraction of binary variables involved in a conflict clause */
   Real             timelimit;          /**< maximal time in seconds to run */
   Real             memlimit;           /**< maximal memory usage in MB */
   Real             gaplimit;           /**< solving stops, if the given gap is reached */
   Real             memsavefac;         /**< fraction of maximal memory usage resulting in switch to memory saving mode */
   Longint          nodelimit;          /**< maximal number of nodes to process (-1: no limit) */
   SCIP*            scip;               /**< very ugly: pointer to scip main data structure for callback methods */
   PARAMSET*        paramset;           /**< set of parameters */
   BUFFER*          buffer;             /**< memory buffers for short living temporary objects */
   READER**         readers;            /**< file readers */
   PRICER**         pricers;            /**< variable pricers */
   CONSHDLR**       conshdlrs;          /**< constraint handlers */
   CONFLICTHDLR**   conflicthdlrs;      /**< conflict handlers */
   PRESOL**         presols;            /**< presolvers */
   SEPA**           sepas;              /**< separators */
   HEUR**           heurs;              /**< primal heuristics */
   EVENTHDLR**      eventhdlrs;         /**< event handlers */
   NODESEL**        nodesels;           /**< node selectors */
   NODESEL*         nodesel;            /**< currently used node selector, or NULL if invalid */
   BRANCHRULE**     branchrules;        /**< branching rules */
   DISP**           disps;              /**< display columns */
   char*            vbcfilename;        /**< name of the VBC Tool output file, or - if no output should be created */
   int              nreaders;           /**< number of file readers */
   int              readerssize;        /**< size of readers array */
   int              npricers;           /**< number of variable pricers */
   int              nactivepricers;     /**< number of variable pricers used in the current problem */
   int              pricerssize;        /**< size of pricers array */
   int              nconshdlrs;         /**< number of constraint handlers */
   int              conshdlrssize;      /**< size of conshdlrs array */
   int              nconflicthdlrs;     /**< number of conflict handlers */
   int              conflicthdlrssize;  /**< size of conflicthdlrs array */
   int              npresols;           /**< number of presolvers */
   int              presolssize;        /**< size of presols array */
   int              nsepas;             /**< number of separators */
   int              sepassize;          /**< size of sepas array */
   int              nheurs;             /**< number of primal heuristics */
   int              heurssize;          /**< size of heurs array */
   int              neventhdlrs;        /**< number of event handlers */
   int              eventhdlrssize;     /**< size of eventhdlrs array */
   int              nnodesels;          /**< number of node selectors */
   int              nodeselssize;       /**< size of nodesels array */
   int              nbranchrules;       /**< number of branching rules */
   int              branchrulessize;    /**< size of branchrules array */
   int              ndisps;             /**< number of display columns */
   int              dispssize;          /**< size of disps array */
   int              memgrowinit;        /**< initial size of dynamically allocated arrays */
   int              treegrowinit;       /**< initial size of tree array */
   int              pathgrowinit;       /**< initial size of path array */
   int              dispwidth;          /**< maximal number of characters in a node information line */
   int              dispfreq;           /**< frequency for displaying node information lines */
   int              dispheaderfreq;     /**< frequency for displaying header lines (every n'th node information line) */
   int              restartbdchgs;      /**< number of root node bound changes triggering a restart with preprocessing
                                         *   (-1: no restart, 0: restart only after complete root node evaluation) */
   int              maxpresolrounds;    /**< maximal number of presolving rounds (-1: unlimited) */
   int              maxproprounds;      /**< maximal number of propagation rounds per node (-1: unlimited) */
   int              maxproproundsroot;  /**< maximal number of propagation rounds in the root node (-1: unlimited) */
   int              maxpricevars;       /**< maximal number of variables priced in per pricing round */
   int              maxpricevarsroot;   /**< maximal number of priced variables at the root node */
   int              maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int              maxsepacutsroot;    /**< maximal number of separated cuts at the root node */
   int              minmaxconfvars;     /**< minimal absolute maximum of variables involved in a conflict clause */
   int              colagelimit;        /**< maximum age a column can reach before it is deleted from the LP */
   int              rowagelimit;        /**< maximum age a row can reach before it is deleted from the LP */
   int              cutagelimit;        /**< maximum age a cut can reach before it is deleted from the global cut pool */
   int              consagelimit;       /**< maximum age an unnecessary constraint can reach before it is deleted, or -1 */
   int              consobsoleteage;    /**< age of a constraint after which it is marked obsolete (not useful anymore) */
   int              maxsol;             /**< maximal number of solutions to store in the solution storage */
   int              sollimit;           /**< solving stops, if the given number of solutions were found (-1: no limit) */
   int              lpsolvefreq;        /**< frequency for solving LP at the nodes (-1: never; 0: only root LP) */
   int              lpsolvedepth;       /**< maximal depth for solving LP at the nodes (-1: no depth limit) */
   int              redcostfreq;        /**< frequency for applying reduced cost fixing (-1: never; 0: only root LP) */
   VERBLEVEL        verblevel;          /**< verbosity level of output */
   CLOCKTYPE        clocktype;          /**< default clock type to use */
   Bool             pricerssorted;      /**< are the pricers sorted by activity and priority? */
   Bool             conflicthdlrssorted;/**< are the conflict handlers sorted by priority? */
   Bool             presolssorted;      /**< are the presolvers sorted by priority? */
   Bool             sepassorted;        /**< are the separators sorted by priority? */
   Bool             heurssorted;        /**< are the heuristics sorted by priority? */
   Bool             branchrulessorted;  /**< are the branching rules sorted by priority? */
   Bool             catchctrlc;         /**< should the CTRL-C interrupt be caught by SCIP? */
   Bool             usepropconflict;    /**< should propagation conflict analysis be used? */
   Bool             uselpconflict;      /**< should infeasible LP conflict analysis be used? */
   Bool             usesbconflict;      /**< should infeasible strong branching conflict analysis be used? */
   Bool             usepseudoconflict;  /**< should pseudo solution conflict analysis be used? */
   Bool             checklpfeas;        /**< should LP solutions be checked, resolving LP when numerical troubles occur? */
   Bool             exactsolve;         /**< should the problem be solved exactly (with proven dual bounds)? */
   Bool             fastmip;            /**< should FASTMIP setting of LP solver be used? */
   Bool             scaling;            /**< should scaling of LP solver be used? */
   Bool             lpinfo;             /**< should the LP solver display status messages? */
   Bool             cleanupcols;        /**< should new non-basic columns be removed after LP solving? */
   Bool             cleanuprows;        /**< should new basic rows be removed after LP solving? */
   Bool             clocksenabled;      /**< is timing enabled? */
   Bool             preferbinbranch;    /**< should branching on binary variables be prefered? */
};


#endif
