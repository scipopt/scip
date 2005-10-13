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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_set.h,v 1.59 2005/10/13 21:10:07 bzfpfend Exp $"

/**@file   struct_set.h
 * @brief  datastructures for global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SET_H__
#define __SCIP_STRUCT_SET_H__


#include "scip/def.h"
#include "scip/message.h"
#include "scip/type_set.h"
#include "scip/type_buffer.h"
#include "scip/type_clock.h"
#include "scip/type_paramset.h"
#include "scip/type_event.h"
#include "scip/type_scip.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_cons.h"
#include "scip/type_disp.h"
#include "scip/type_heur.h"
#include "scip/type_nodesel.h"
#include "scip/type_presol.h"
#include "scip/type_pricer.h"
#include "scip/type_reader.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_prop.h"


/** global SCIP settings */
struct SCIP_Set
{
   SCIP_STAGE            stage;              /**< SCIP operation stage */
   SCIP*                 scip;               /**< very ugly: pointer to scip main data structure for callback methods */
   SCIP_PARAMSET*        paramset;           /**< set of parameters */
   SCIP_BUFFER*          buffer;             /**< memory buffers for short living temporary objects */
   SCIP_READER**         readers;            /**< file readers */
   SCIP_PRICER**         pricers;            /**< variable pricers */
   SCIP_CONSHDLR**       conshdlrs;          /**< constraint handlers (sorted by check priority) */
   SCIP_CONSHDLR**       conshdlrs_sepa;     /**< constraint handlers (sorted by separation priority) */
   SCIP_CONSHDLR**       conshdlrs_enfo;     /**< constraint handlers (sorted by enforcement priority) */
   SCIP_CONFLICTHDLR**   conflicthdlrs;      /**< conflict handlers */
   SCIP_PRESOL**         presols;            /**< presolvers */
   SCIP_RELAX**          relaxs;             /**< relaxators */
   SCIP_SEPA**           sepas;              /**< separators */
   SCIP_PROP**           props;              /**< propagators */
   SCIP_HEUR**           heurs;              /**< primal heuristics */
   SCIP_EVENTHDLR**      eventhdlrs;         /**< event handlers */
   SCIP_NODESEL**        nodesels;           /**< node selectors */
   SCIP_NODESEL*         nodesel;            /**< currently used node selector, or NULL if invalid */
   SCIP_BRANCHRULE**     branchrules;        /**< branching rules */
   SCIP_DISP**           disps;              /**< display columns */
   int                   nreaders;           /**< number of file readers */
   int                   readerssize;        /**< size of readers array */
   int                   npricers;           /**< number of variable pricers */
   int                   nactivepricers;     /**< number of variable pricers used in the current problem */
   int                   pricerssize;        /**< size of pricers array */
   int                   nconshdlrs;         /**< number of constraint handlers */
   int                   conshdlrssize;      /**< size of conshdlrs array */
   int                   nconflicthdlrs;     /**< number of conflict handlers */
   int                   conflicthdlrssize;  /**< size of conflicthdlrs array */
   int                   npresols;           /**< number of presolvers */
   int                   presolssize;        /**< size of presols array */
   int                   nrelaxs;            /**< number of relaxators */
   int                   relaxssize;         /**< size of relaxs array */
   int                   nsepas;             /**< number of separators */
   int                   sepassize;          /**< size of sepas array */
   int                   nprops;             /**< number of propagators */
   int                   propssize;          /**< size of props array */
   int                   nheurs;             /**< number of primal heuristics */
   int                   heurssize;          /**< size of heurs array */
   int                   neventhdlrs;        /**< number of event handlers */
   int                   eventhdlrssize;     /**< size of eventhdlrs array */
   int                   nnodesels;          /**< number of node selectors */
   int                   nodeselssize;       /**< size of nodesels array */
   int                   nbranchrules;       /**< number of branching rules */
   int                   branchrulessize;    /**< size of branchrules array */
   int                   ndisps;             /**< number of display columns */
   int                   dispssize;          /**< size of disps array */
   SCIP_Bool             pricerssorted;      /**< are the pricers sorted by activity and priority? */
   SCIP_Bool             conflicthdlrssorted;/**< are the conflict handlers sorted by priority? */
   SCIP_Bool             presolssorted;      /**< are the presolvers sorted by priority? */
   SCIP_Bool             relaxssorted;       /**< are the relaxators sorted by priority? */
   SCIP_Bool             sepassorted;        /**< are the separators sorted by priority? */
   SCIP_Bool             propssorted;        /**< are the propagators sorted by priority? */
   SCIP_Bool             heurssorted;        /**< are the heuristics sorted by priority? */
   SCIP_Bool             branchrulessorted;  /**< are the branching rules sorted by priority? */

   /* branching settings */
   char                  branch_scorefunc;   /**< branching score function ('s'um, 'p'roduct) */
   SCIP_Real             branch_scorefac;    /**< branching score factor to weigh downward and upward gain prediction
                                              *   in sum score function */
   SCIP_Bool             branch_preferbinary;/**< should branching on binary variables be preferred? */

   /* conflict analysis settings */
   SCIP_Real             conf_maxvarsfac;    /**< maximal fraction of binary variables involved in a conflict clause */
   int                   conf_minmaxvars;    /**< minimal absolute maximum of variables involved in a conflict clause */
   int                   conf_maxunfixed;    /**< maximal number of unfixed variables at insertion depth of conflict clause
                                              *   (-1: no limit) */
   int                   conf_maxlploops;    /**< maximal number of LP resolving loops during conflict analysis */
   int                   conf_fuiplevels;    /**< number of depth levels up to which first UIP's are used in conflict
                                              *   analysis (-1: use All-FirstUIP rule) */
   int                   conf_interclauses;  /**< maximal number of intermediate conflict clauses generated in conflict
                                              *   graph (-1: use every intermediate clause) */
   int                   conf_maxclauses;    /**< maximal number of conflict clauses accepted at an infeasible node
                                              *   (-1: use all generated conflict clauses) */
   SCIP_Bool             conf_reconvclauses; /**< should reconvergence clauses be created for UIPs of last depth level? */
   SCIP_Bool             conf_useprop;       /**< should propagation conflict analysis be used? */
   SCIP_Bool             conf_uselp;         /**< should infeasible LP conflict analysis be used? */
   SCIP_Bool             conf_usesb;         /**< should infeasible strong branching conflict analysis be used? */
   SCIP_Bool             conf_usepseudo;     /**< should pseudo solution conflict analysis be used? */
   SCIP_Bool             conf_allowlocal;    /**< should conflict clauses be generated that are only valid locally? */
   SCIP_Bool             conf_repropagate;   /**< should earlier nodes be repropagated in order to replace branching
                                              *   decisions by deductions */
   SCIP_Bool             conf_dynamic;       /**< should the conflict constraints be subject to aging? */
   SCIP_Bool             conf_removeable;    /**< should the conflict's relaxations be subject to LP aging and cleanup? */
   SCIP_Real             conf_scorefac;      /**< factor to decrease importance of variables' earlier conflict scores */

   /* constraint settings */
   int                   cons_agelimit;      /**< maximum age an unnecessary constraint can reach before it is deleted
                                              *   (0: dynamic, -1: disable aging) */
   int                   cons_obsoleteage;   /**< age of a constraint after which it is marked obsolete
                                              *   (0: dynamic, -1: disable obsoletion) */

   /* display settings */
   SCIP_VERBLEVEL        disp_verblevel;     /**< verbosity level of output */
   int                   disp_width;         /**< maximal number of characters in a node information line */
   int                   disp_freq;          /**< frequency for displaying node information lines */
   int                   disp_headerfreq;    /**< frequency for displaying header lines (every n'th node information line) */
   SCIP_Bool             disp_lpinfo;        /**< should the LP solver display status messages? */

   /* limit settings */
   SCIP_Real             limit_time;         /**< maximal time in seconds to run */
   SCIP_Real             limit_memory;       /**< maximal memory usage in MB */
   SCIP_Real             limit_gap;          /**< solving stops, if the given gap is reached */
   SCIP_Longint          limit_nodes;        /**< maximal number of nodes to process (-1: no limit) */
   int                   limit_solutions;    /**< solving stops, if the given number of solutions were found (-1: no limit) */
   int                   limit_bestsol;      /**< solving stops, if the given number of solution improvements were found
                                              *   (-1: no limit) */
   int                   limit_maxsol;       /**< maximal number of solutions to store in the solution storage */

   /* LP settings */
   int                   lp_solvefreq;       /**< frequency for solving LP at the nodes (-1: never; 0: only root LP) */
   int                   lp_solvedepth;      /**< maximal depth for solving LP at the nodes (-1: no depth limit) */
   char                  lp_initalgorithm;   /**< LP algorithm for solving initial LP relaxations ('s'implex, 'b'arrier,
                                              *   barrier with 'c'rossover) */
   char                  lp_resolvealgorithm;/**< LP algorithm for resolving LP relaxations if a starting basis exists
                                              *   ('s'implex, 'b'arrier, barrier with 'c'rossover) */
   char                  lp_pricing;         /**< LP pricing strategy ('a'uto, 'f'ull pricing, 's'teepest edge pricing,
                                              *   'q'uickstart steepest edge pricing, 'd'evex pricing) */
   int                   lp_colagelimit;     /**< maximum age a column can reach before it is deleted from the SCIP_LP
                                              *   (-1: don't delete columns due to aging) */
   int                   lp_rowagelimit;     /**< maximum age a row can reach before it is deleted from the LP 
                                              *   (-1: don't delete rows due to aging) */
   SCIP_Bool             lp_cleanupcols;     /**< should new non-basic columns be removed after LP solving? */
   SCIP_Bool             lp_cleanupcolsroot; /**< should new non-basic columns be removed after root LP solving? */
   SCIP_Bool             lp_cleanuprows;     /**< should new basic rows be removed after LP solving? */
   SCIP_Bool             lp_cleanuprowsroot; /**< should new basic rows be removed after root LP solving? */
   SCIP_Bool             lp_checkstability;  /**< should LP solver's return status be checked for stability? */
   SCIP_Bool             lp_checkfeas;       /**< should LP solutions be checked, resolving LP when numerical troubles occur? */
   SCIP_Bool             lp_fastmip;         /**< should FASTMIP setting of LP solver be used? */
   SCIP_Bool             lp_scaling;         /**< should scaling of LP solver be used? */
   SCIP_Bool             lp_presolving;      /**< should presolving of LP solver be used? */

   /* memory settings */
   SCIP_Real             mem_savefac;        /**< fraction of maximal memory usage resulting in switch to memory saving mode */
   SCIP_Real             mem_arraygrowfac;   /**< memory growing factor for dynamically allocated arrays */
   SCIP_Real             mem_treegrowfac;    /**< memory growing factor for tree array */
   SCIP_Real             mem_pathgrowfac;    /**< memory growing factor for path array */
   int                   mem_arraygrowinit;  /**< initial size of dynamically allocated arrays */
   int                   mem_treegrowinit;   /**< initial size of tree array */
   int                   mem_pathgrowinit;   /**< initial size of path array */
   
   /* miscellaneous settings */
   SCIP_Bool             misc_catchctrlc;    /**< should the CTRL-C interrupt be caught by SCIP? */
   SCIP_Bool             misc_exactsolve;    /**< should the problem be solved exactly (with proven dual bounds)? */

   /* numerical settings */
   SCIP_Real             num_infinity;       /**< values larger than this are considered infinity */
   SCIP_Real             num_epsilon;        /**< absolute values smaller than this are considered zero */
   SCIP_Real             num_sumepsilon;     /**< absolute values of sums smaller than this are considered zero */
   SCIP_Real             num_feastol;        /**< feasibility tolerance for constraints */
   SCIP_Real             num_dualfeastol;    /**< feasibility tolerance for reduced costs */
   SCIP_Real             num_barrierconvtol; /**< convergence tolerance used in barrier algorithm */
   SCIP_Real             num_boundstreps;    /**< minimal improve for strengthening bounds */
   SCIP_Real             num_pseudocosteps;  /**< minimal variable distance value to use for pseudo cost updates */
   SCIP_Real             num_pseudocostdelta;/**< minimal objective distance value to use for pseudo cost updates */

   /* presolving settings */
   SCIP_Real             presol_abortfac;    /**< abort presolve, if l.t. this frac of the problem was changed in last round */
   int                   presol_maxrounds;   /**< maximal number of presolving rounds (-1: unlimited) */
   int                   presol_maxrestarts; /**< maximal number of restarts (-1: unlimited) */
   SCIP_Real             presol_restartfac;  /**< fraction of bounds that changed in the root node triggering a restart with
                                              *   preprocessing (0.0: restart only after complete root node evaluation) */

   /* pricing settings */
   SCIP_Real             price_abortfac;     /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
   int                   price_maxvars;      /**< maximal number of variables priced in per pricing round */
   int                   price_maxvarsroot;  /**< maximal number of priced variables at the root node */

   /* propagation settings */
   int                   prop_maxrounds;     /**< maximal number of propagation rounds per node (-1: unlimited) */
   int                   prop_maxroundsroot; /**< maximal number of propagation rounds in the root node (-1: unlimited) */
   int                   prop_redcostfreq;   /**< frequency for applying reduced cost fixing (-1: never; 0: only root LP) */

   /* separation settings */
   SCIP_Real             sepa_maxbounddist;  /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying separation
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_Real             sepa_minefficacy;   /**< minimal efficacy for a cut to enter the LP */
   SCIP_Real             sepa_minefficacyroot; /**< minimal efficacy for a cut to enter the LP in the root node */
   SCIP_Real             sepa_minortho;      /**< minimal orthogonality for a cut to enter the LP */
   SCIP_Real             sepa_minorthoroot;  /**< minimal orthogonality for a cut to enter the LP in the root node */
   SCIP_Real             sepa_objparalfac;   /**< factor to scale objective parallelism of cut in separation score calc. */
   SCIP_Real             sepa_orthofac;      /**< factor to scale orthogonality of cut in separation score calculation */
   char                  sepa_efficacynorm;  /**< row norm to use for efficacy calculation ('e'uclidean, 'm'aximum, 's'um,
                                              *   'd'iscrete) */
   int                   sepa_maxrounds;     /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   sepa_maxroundsroot; /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   sepa_maxaddrounds;  /**< maximal additional number of separation rounds in subsequent price-and-cut
                                              *   loops (-1: no additional restriction) */
   int                   sepa_maxstallrounds;/**< maximal number of consecutive separation rounds without objective
                                              *   or integrality improvement (-1: no additional restriction) */
   int                   sepa_maxcuts;       /**< maximal number of cuts separated per separation round */
   int                   sepa_maxcutsroot;   /**< maximal number of separated cuts at the root node */
   int                   sepa_cutagelimit;   /**< maximum age a cut can reach before it is deleted from the global cut pool */
   int                   sepa_poolfreq;      /**< separation frequency for the global cut pool */

   /* timing settings */
   SCIP_CLOCKTYPE        time_clocktype;     /**< default clock type to use */
   SCIP_Bool             time_enabled;       /**< is timing enabled? */

   /* VBC tool settings */
   char*                 vbc_filename;       /**< name of the VBC Tool output file, or - if no output should be created */
   SCIP_Bool             vbc_realtime;       /**< should the real solving time be used instead of time step counter in VBC output? */
};


#endif
