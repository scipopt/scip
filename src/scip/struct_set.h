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
#pragma ident "@(#) $Id: struct_set.h,v 1.37 2004/11/17 13:09:48 bzfpfend Exp $"

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
#include "type_relax.h"
#include "type_sepa.h"
#include "type_prop.h"


/** global SCIP settings */
struct Set
{
   SCIP*            scip;               /**< very ugly: pointer to scip main data structure for callback methods */
   PARAMSET*        paramset;           /**< set of parameters */
   BUFFER*          buffer;             /**< memory buffers for short living temporary objects */
   READER**         readers;            /**< file readers */
   PRICER**         pricers;            /**< variable pricers */
   CONSHDLR**       conshdlrs;          /**< constraint handlers */
   CONFLICTHDLR**   conflicthdlrs;      /**< conflict handlers */
   PRESOL**         presols;            /**< presolvers */
   RELAX**          relaxs;             /**< relaxators */
   SEPA**           sepas;              /**< separators */
   PROP**           props;              /**< propagators */
   HEUR**           heurs;              /**< primal heuristics */
   EVENTHDLR**      eventhdlrs;         /**< event handlers */
   NODESEL**        nodesels;           /**< node selectors */
   NODESEL*         nodesel;            /**< currently used node selector, or NULL if invalid */
   BRANCHRULE**     branchrules;        /**< branching rules */
   DISP**           disps;              /**< display columns */
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
   int              nrelaxs;            /**< number of relaxators */
   int              relaxssize;         /**< size of relaxs array */
   int              nsepas;             /**< number of separators */
   int              sepassize;          /**< size of sepas array */
   int              nprops;             /**< number of propagators */
   int              propssize;          /**< size of props array */
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
   Bool             pricerssorted;      /**< are the pricers sorted by activity and priority? */
   Bool             conflicthdlrssorted;/**< are the conflict handlers sorted by priority? */
   Bool             presolssorted;      /**< are the presolvers sorted by priority? */
   Bool             relaxssorted;       /**< are the relaxators sorted by priority? */
   Bool             sepassorted;        /**< are the separators sorted by priority? */
   Bool             propssorted;        /**< are the propagators sorted by priority? */
   Bool             heurssorted;        /**< are the heuristics sorted by priority? */
   Bool             branchrulessorted;  /**< are the branching rules sorted by priority? */

   /* branching settings */
   Real             branch_scorefac;    /**< branching score factor to weigh downward and upward gain prediction */
   Bool             branch_preferbinary;/**< should branching on binary variables be prefered? */

   /* conflict analysis settings */
   Real             conf_maxvarsfac;    /**< maximal fraction of binary variables involved in a conflict clause */
   int              conf_minmaxvars;    /**< minimal absolute maximum of variables involved in a conflict clause */
   int              conf_maxlploops;    /**< maximal number of LP resolving loops during conflict analysis */
   int              conf_fuiplevels;    /**< number of depth levels up to which first UIP's are used in conflict
                                         *   analysis (-1: use All-FirstUIP rule) */
   int              conf_interclauses;  /**< maximal number of intermediate conflict clauses generated in conflict
                                         *   graph (-1: use every intermediate clause) */
   Bool             conf_reconvclauses; /**< should reconvergence clauses be created for UIPs of last depth level? */
   Bool             conf_useprop;       /**< should propagation conflict analysis be used? */
   Bool             conf_uselp;         /**< should infeasible LP conflict analysis be used? */
   Bool             conf_usesb;         /**< should infeasible strong branching conflict analysis be used? */
   Bool             conf_usepseudo;     /**< should pseudo solution conflict analysis be used? */
   Bool             conf_repropagate;   /**< should earlier nodes be repropagated in order to replace branching
                                         *   decisions by deductions */

   /* constraint settings */
   int              cons_agelimit;      /**< maximum age an unnecessary constraint can reach before it is deleted, or -1 */
   int              cons_obsoleteage;   /**< age of a constraint after which it is marked obsolete (not useful anymore) */

   /* display settings */
   VERBLEVEL        disp_verblevel;     /**< verbosity level of output */
   int              disp_width;         /**< maximal number of characters in a node information line */
   int              disp_freq;          /**< frequency for displaying node information lines */
   int              disp_headerfreq;    /**< frequency for displaying header lines (every n'th node information line) */
   Bool             disp_lpinfo;        /**< should the LP solver display status messages? */

   /* limit settings */
   Real             limit_time;         /**< maximal time in seconds to run */
   Real             limit_memory;       /**< maximal memory usage in MB */
   Real             limit_gap;          /**< solving stops, if the given gap is reached */
   Longint          limit_nodes;        /**< maximal number of nodes to process (-1: no limit) */
   int              limit_sol;          /**< solving stops, if the given number of solutions were found (-1: no limit) */
   int              limit_bestsol;      /**< solving stops, if the given number of solution improvements were found
                                         *   (-1: no limit) */
   int              limit_maxsol;       /**< maximal number of solutions to store in the solution storage */

   /* LP settings */
   int              lp_solvefreq;       /**< frequency for solving LP at the nodes (-1: never; 0: only root LP) */
   int              lp_solvedepth;      /**< maximal depth for solving LP at the nodes (-1: no depth limit) */
   int              lp_colagelimit;     /**< maximum age a column can reach before it is deleted from the LP
                                         *   (-1: don't delete columns due to aging) */
   int              lp_rowagelimit;     /**< maximum age a row can reach before it is deleted from the LP 
                                         *   (-1: don't delete rows due to aging) */
   Bool             lp_cleanupcols;     /**< should new non-basic columns be removed after LP solving? */
   Bool             lp_cleanuprows;     /**< should new basic rows be removed after LP solving? */
   Bool             lp_checkfeas;       /**< should LP solutions be checked, resolving LP when numerical troubles occur? */
   Bool             lp_fastmip;         /**< should FASTMIP setting of LP solver be used? */
   Bool             lp_scaling;         /**< should scaling of LP solver be used? */
   Bool             lp_presolving;      /**< should presolving of LP solver be used? */

   /* memory settings */
   Real             mem_savefac;        /**< fraction of maximal memory usage resulting in switch to memory saving mode */
   Real             mem_arraygrowfac;   /**< memory growing factor for dynamically allocated arrays */
   Real             mem_treegrowfac;    /**< memory growing factor for tree array */
   Real             mem_pathgrowfac;    /**< memory growing factor for path array */
   int              mem_arraygrowinit;  /**< initial size of dynamically allocated arrays */
   int              mem_treegrowinit;   /**< initial size of tree array */
   int              mem_pathgrowinit;   /**< initial size of path array */
   
   /* miscellaneous settings */
   Bool             misc_catchctrlc;    /**< should the CTRL-C interrupt be caught by SCIP? */
   Bool             misc_exactsolve;    /**< should the problem be solved exactly (with proven dual bounds)? */

   /* numerical settings */
   Real             num_infinity;       /**< values larger than this are considered infinity */
   Real             num_epsilon;        /**< absolute values smaller than this are considered zero */
   Real             num_sumepsilon;     /**< absolute values of sums smaller than this are considered zero */
   Real             num_feastol;        /**< feasibility tolerance for constraints */
   Real             num_dualfeastol;    /**< feasibility tolerance for reduced costs */
   Real             num_boundstreps;    /**< minimal improve for strengthening bounds */
   Real             num_pseudocosteps;  /**< minimal variable distance value to use for pseudo cost updates */
   Real             num_pseudocostdelta;/**< minimal objective distance value to use for pseudo cost updates */

   /* presolving settings */
   Real             presol_abortfac;    /**< abort presolve, if l.t. this frac of the problem was changed in last round */
   int              presol_maxrounds;   /**< maximal number of presolving rounds (-1: unlimited) */
   int              presol_restartbdchgs; /**< number of root node bound changes triggering a restart with preprocessing
                                           *   (-1: no restart, 0: restart only after complete root node evaluation) */

   /* pricing settings */
   Real             price_abortfac;     /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
   int              price_maxvars;      /**< maximal number of variables priced in per pricing round */
   int              price_maxvarsroot;  /**< maximal number of priced variables at the root node */

   /* propagation settings */
   int              prop_maxrounds;     /**< maximal number of propagation rounds per node (-1: unlimited) */
   int              prop_maxroundsroot; /**< maximal number of propagation rounds in the root node (-1: unlimited) */
   int              prop_redcostfreq;   /**< frequency for applying reduced cost fixing (-1: never; 0: only root LP) */

   /* separation settings */
   Real             sepa_maxbounddist;  /**< maximal relative distance from current node's dual bound to primal bound
                                         *   compared to best node's dual bound for applying separation
                                         *   (0.0: only on current best node, 1.0: on all nodes) */
   Real             sepa_minefficacy;   /**< minimal efficacy for a cut to enter the LP */
   Real             sepa_minefficacyroot; /**< minimal efficacy for a cut to enter the LP in the root node */
   Real             sepa_minortho;      /**< minimal orthogonality for a cut to enter the LP */
   Real             sepa_minorthoroot;  /**< minimal orthogonality for a cut to enter the LP in the root node */
   Real             sepa_orthofac;      /**< factor to scale orthogonality of cut in separation score calculation */
   char             sepa_efficacynorm;  /**< row norm to use for efficacy calculation ('e'uclidean, 'm'aximum, 's'um,
                                         *   'd'iscrete) */
   int              sepa_maxrounds;     /**< maximal number of separation rounds per node (-1: unlimited) */
   int              sepa_maxroundsroot; /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int              sepa_maxaddrounds;  /**< maximal additional number of separation rounds in subsequent price-and-cut
                                         *   loops (-1: no additional restriction) */
   int              sepa_maxstallrounds;/**< maximal number of consecutive separation rounds without objective
                                         *   improvement (-1: no additional restriction) */
   int              sepa_maxcuts;       /**< maximal number of cuts separated per separation round */
   int              sepa_maxcutsroot;   /**< maximal number of separated cuts at the root node */
   int              sepa_cutagelimit;   /**< maximum age a cut can reach before it is deleted from the global cut pool */
   int              sepa_poolfreq;      /**< separation frequency for the global cut pool */

   /* timing settings */
   CLOCKTYPE        time_clocktype;     /**< default clock type to use */
   Bool             time_enabled;       /**< is timing enabled? */

   /* VBC tool settings */
   char*            vbc_filename;       /**< name of the VBC Tool output file, or - if no output should be created */
   Bool             vbc_realtime;       /**< should the real solving time be used instead of time step counter in VBC output? */
};


#endif
