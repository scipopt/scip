/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   set.h
 * @brief  global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SET_H__
#define __SET_H__


/** possible settings for enabling/disabling algorithms and other features */
enum Setting
{
   SCIP_UNDEFINED = 0,                  /**< undefined setting */
   SCIP_DISABLED  = 1,                  /**< feature is disabled */
   SCIP_AUTO      = 2,                  /**< feature is set to automatic mode */
   SCIP_ENABLED   = 3                   /**< feature is enabled */
};
typedef enum Setting SETTING;

typedef struct Set SET;                 /**< global SCIP settings */



#include <math.h>

#include "def.h"
#include "stat.h"
#include "misc.h"
#include "scip.h"
#include "paramset.h"
#include "buffer.h"
#include "reader.h"
#include "pricer.h"
#include "cons.h"
#include "conflict.h"
#include "presol.h"
#include "sepa.h"
#include "heur.h"
#include "event.h"
#include "nodesel.h"
#include "branch.h"
#include "disp.h"
#include "lp.h"
#include "message.h"


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
   int              pricerssize;        /**< size of pricers array */
   CONSHDLR**       conshdlrs;          /**< constraint handlers */
   int              nconshdlrs;         /**< number of constraint handlers */
   int              conshdlrssize;      /**< size of conshdlrs array */
   CONFLICTHDLR**   conflicthdlrs;      /**< conflict handlers */
   int              nconflicthdlrs;     /**< number of conflict handlers */
   int              conflicthdlrssize;  /**< size of conflicthdlrs array */
   PRESOL**         presols;            /**< presolvers */
   int              npresols;           /**< number of presolvers */
   int              presolssize;        /**< size of presols array */
   SEPA**           sepas;              /**< separators */
   int              nsepas;             /**< number of separators */
   int              sepassize;          /**< size of sepas array */
   HEUR**           heurs;              /**< primal heuristics */
   int              nheurs;             /**< number of primal heuristics */
   int              heurssize;          /**< size of heurs array */
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
   Bool             cleanupcols;        /**< should new non-basic columns be removed after LP solving? */
   Bool             cleanuprows;        /**< should new basic rows be removed after LP solving? */
   CLOCKTYPE        clocktype;          /**< default clock type to use */
   Bool             clocksenabled;      /**< is timing enabled? */
};


/** creates global SCIP settings */
extern
RETCODE SCIPsetCreate(
   SET**            set,                /**< pointer to SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** frees global SCIP settings */
extern
RETCODE SCIPsetFree(
   SET**            set,                /**< pointer to SCIP settings */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a Bool parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPsetAddBoolParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Bool             defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPsetAddIntParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   int*             valueptr,           /**< pointer to store the current parameter value, or NULL */
   int              defaultvalue,       /**< default value of the parameter */
   int              minvalue,           /**< minimum value for parameter */
   int              maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a Longint parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPsetAddLongintParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   Longint          defaultvalue,       /**< default value of the parameter */
   Longint          minvalue,           /**< minimum value for parameter */
   Longint          maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a Real parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPsetAddRealParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   Real             defaultvalue,       /**< default value of the parameter */
   Real             minvalue,           /**< minimum value for parameter */
   Real             maxvalue,           /**< maximum value for parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPsetAddCharParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   char             defaultvalue,       /**< default value of the parameter */
   const char*      allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPsetAddStringParam(
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of the parameter */
   const char*      desc,               /**< description of the parameter */
   char**           valueptr,           /**< pointer to store the current parameter value, or NULL */
   const char*      defaultvalue,       /**< default value of the parameter */
   DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** gets the value of an existing Bool parameter */
RETCODE SCIPsetGetBoolParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Bool*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Int parameter */
extern
RETCODE SCIPsetGetIntParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   int*             value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Longint parameter */
extern
RETCODE SCIPsetGetLongintParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Longint*         value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Real parameter */
extern
RETCODE SCIPsetGetRealParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Real*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing Char parameter */
extern
RETCODE SCIPsetGetCharParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   char*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing String parameter */
extern
RETCODE SCIPsetGetStringParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   char**           value               /**< pointer to store the parameter */
   );

/** changes the value of an existing Bool parameter */
extern
RETCODE SCIPsetSetBoolParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing Int parameter */
extern
RETCODE SCIPsetSetIntParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   int              value               /**< new value of the parameter */
   );

/** changes the value of an existing Longint parameter */
extern
RETCODE SCIPsetSetLongintParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing Real parameter */
extern
RETCODE SCIPsetSetRealParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing Char parameter */
extern
RETCODE SCIPsetSetCharParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   char             value               /**< new value of the parameter */
   );

/** changes the value of an existing String parameter */
extern
RETCODE SCIPsetSetStringParam(
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of the parameter */
   const char*      value               /**< new value of the parameter */
   );

/** reads parameters from a file */
extern
RETCODE SCIPsetReadParams(
   SET*             set,                /**< global SCIP settings */
   const char*      filename            /**< file name */
   );

/** writes all parameters in the parameter set to a file */
extern
RETCODE SCIPsetWriteParams(
   SET*             set,                /**< global SCIP settings */
   const char*      filename,           /**< file name, or NULL for stdout */
   Bool             comments            /**< should parameter descriptions be written as comments? */
   );

/** returns the array of all available SCIP parameters */
extern
PARAM** SCIPsetGetParams(
   SET*             set                 /**< global SCIP settings */
   );

/** returns the total number of all available SCIP parameters */
extern
int SCIPsetGetNParams(
   SET*             set                 /**< global SCIP settings */
   );

/** inserts file reader in file reader list */
extern
RETCODE SCIPsetIncludeReader(
   SET*             set,                /**< global SCIP settings */
   READER*          reader              /**< file reader */
   );

/** returns the file reader of the given name, or NULL if not existing */
extern
READER* SCIPsetFindReader(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of file reader */
   );

/** inserts variable pricer in variable pricer list */
extern
RETCODE SCIPsetIncludePricer(
   SET*             set,                /**< global SCIP settings */
   PRICER*          pricer              /**< variable pricer */
   );

/** returns the variable pricer of the given name, or NULL if not existing */
extern
PRICER* SCIPsetFindPricer(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of variable pricer */
   );

/** inserts constraint handler in constraint handler list */
extern
RETCODE SCIPsetIncludeConshdlr(
   SET*             set,                /**< global SCIP settings */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** returns the constraint handler of the given name, or NULL if not existing */
extern
CONSHDLR* SCIPsetFindConshdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of constraint handler */
   );

/** inserts conflict handler in conflict handler list */
extern
RETCODE SCIPsetIncludeConflicthdlr(
   SET*             set,                /**< global SCIP settings */
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   );

/** returns the conflict handler of the given name, or NULL if not existing */
extern
CONFLICTHDLR* SCIPsetFindConflicthdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of conflict handler */
   );

/** inserts presolver in presolver list */
extern
RETCODE SCIPsetIncludePresol(
   SET*             set,                /**< global SCIP settings */
   PRESOL*          presol              /**< presolver */
   );

/** returns the presolver of the given name, or NULL if not existing */
extern
PRESOL* SCIPsetFindPresol(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of presolver */
   );

/** inserts separator in separator list */
extern
RETCODE SCIPsetIncludeSepa(
   SET*             set,                /**< global SCIP settings */
   SEPA*            sepa                /**< separator */
   );

/** returns the separator of the given name, or NULL if not existing */
extern
SEPA* SCIPsetFindSepa(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of separator */
   );

/** inserts primal heuristic in primal heuristic list */
extern
RETCODE SCIPsetIncludeHeur(
   SET*             set,                /**< global SCIP settings */
   HEUR*            heur                /**< primal heuristic */
   );

/** returns the primal heuristic of the given name, or NULL if not existing */
extern
HEUR* SCIPsetFindHeur(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of primal heuristic */
   );

/** inserts event handler in event handler list */
extern
RETCODE SCIPsetIncludeEventhdlr(
   SET*             set,                /**< global SCIP settings */
   EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** returns the event handler of the given name, or NULL if not existing */
extern
EVENTHDLR* SCIPsetFindEventhdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   );

/** inserts node selector in node selector list */
extern
RETCODE SCIPsetIncludeNodesel(
   SET*             set,                /**< global SCIP settings */
   NODESEL*         nodesel             /**< node selector */
   );

/** returns the node selector of the given name, or NULL if not existing */
extern
NODESEL* SCIPsetFindNodesel(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   );

/** returns node selector with highest priority in the current mode */
extern
NODESEL* SCIPsetGetActNodesel(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** inserts branching rule in branching rule list */
extern
RETCODE SCIPsetIncludeBranchrule(
   SET*             set,                /**< global SCIP settings */
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** returns the branching rule of the given name, or NULL if not existing */
extern
BRANCHRULE* SCIPsetFindBranchrule(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   );

/** sorts branching rules by priorities */
extern
void SCIPsetSortBranchrules(
   SET*             set                 /**< global SCIP settings */
   );

/** inserts display column in display column list */
extern
RETCODE SCIPsetIncludeDisp(
   SET*             set,                /**< global SCIP settings */
   DISP*            disp                /**< display column */
   );

/** returns the display column of the given name, or NULL if not existing */
extern
DISP* SCIPsetFindDisp(
   const SET*       set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   );

/** initializes all user callback functions */
extern
RETCODE SCIPsetInitCallbacks(
   const SET*       set                 /**< global SCIP settings */
   );

/** calls exit methods of all user callback functions */
extern
RETCODE SCIPsetExitCallbacks(
   const SET*       set                 /**< global SCIP settings */
   );

/** calculate memory size for dynamically allocated arrays */
extern
int SCIPsetCalcMemGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** calculate memory size for tree array */
extern
int SCIPsetCalcTreeGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** calculate memory size for path array */
extern
int SCIPsetCalcPathGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** sets verbosity level for message output */
extern
RETCODE SCIPsetSetVerbLevel(
   SET*             set,                /**< global SCIP settings */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   );

/** sets LP feasibility tolerance */
extern
RETCODE SCIPsetSetFeastol(
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data (or NULL) */
   Real             feastol             /**< new feasibility tolerance */
   );

/** returns the maximal number of cuts separated per round */
extern
int SCIPsetGetMaxsepacuts(
   const SET*       set,                /**< global SCIP settings */
   Bool             root                /**< are we at the root node? */
   );

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
extern
Real SCIPsetRelDiff(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );



#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** checks, if values are in range of epsilon */
extern
Bool SCIPsetIsEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) lower than val2 */
extern
Bool SCIPsetIsLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) greater than val2 */
extern
Bool SCIPsetIsLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) greater than val2 */
extern
Bool SCIPsetIsGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) lower than val2 */
extern
Bool SCIPsetIsGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range epsilon of 0.0 */
extern
Bool SCIPsetIsZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than epsilon */
extern
Bool SCIPsetIsPositive(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -epsilon */
extern
Bool SCIPsetIsNegative(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if values are in range of sumepsilon */
extern
Bool SCIPsetIsSumEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) lower than val2 */
extern
Bool SCIPsetIsSumLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
extern
Bool SCIPsetIsSumLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) greater than val2 */
extern
Bool SCIPsetIsSumGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
extern
Bool SCIPsetIsSumGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range sumepsilon of 0.0 */
extern
Bool SCIPsetIsSumZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than sumepsilon */
extern
Bool SCIPsetIsSumPositive(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -sumepsilon */
extern
Bool SCIPsetIsSumNegative(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if values are in range of feasibility tolerance */
extern
Bool SCIPsetIsFeasEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPsetIsFeasLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPsetIsFeasLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPsetIsFeasGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPsetIsFeasGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range feasibility tolerance of 0.0 */
extern
Bool SCIPsetIsFeasZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than feasibility tolerance */
extern
Bool SCIPsetIsFeasPositive(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -feasibility tolerance */
extern
Bool SCIPsetIsFeasNegative(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if the cut's activity is more then cutvioleps larger than the given right hand side;
 *  both, the activity and the rhs, should be normed
 */
extern
Bool SCIPsetIsCutViolated(
   const SET*       set,                /**< global SCIP settings */
   Bool             root,               /**< should the root's cutvioleps be used? */
   Real             cutactivity,        /**< activity of the cut */
   Real             cutrhs              /**< right hand side value of the cut */
   );

/** checks, if relative difference of values is in range of epsilon */
extern
Bool SCIPsetIsRelEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than epsilon */
extern
Bool SCIPsetIsRelLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
extern
Bool SCIPsetIsRelLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than epsilon */
extern
Bool SCIPsetIsRelGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
extern
Bool SCIPsetIsRelGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of values is in range of sumepsilon */
extern
Bool SCIPsetIsSumRelEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
extern
Bool SCIPsetIsSumRelLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
extern
Bool SCIPsetIsSumRelLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
extern
Bool SCIPsetIsSumRelGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
extern
Bool SCIPsetIsSumRelGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is (positive) infinite */
extern
Bool SCIPsetIsInfinity(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   );

/** checks, if value is non-negative within the LP feasibility bounds */
extern
Bool SCIPsetIsFeasible(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is integral within the LP feasibility bounds */
extern
Bool SCIPsetIsIntegral(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if given fractional part is smaller than feastol */
extern
Bool SCIPsetIsFracIntegral(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value + feasibility tolerance down to the next integer */
extern
Real SCIPsetFloor(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value - feasibility tolerance up to the next integer */
extern
Real SCIPsetCeil(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. ceil(x) - x */
extern
Real SCIPsetFrac(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to return fractional part for */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsetIsEQ(set, val1, val2)       ( EPSEQ(val1, val2, (set)->epsilon) )
#define SCIPsetIsLT(set, val1, val2)       ( EPSLT(val1, val2, (set)->epsilon) )
#define SCIPsetIsLE(set, val1, val2)       ( EPSLE(val1, val2, (set)->epsilon) )
#define SCIPsetIsGT(set, val1, val2)       ( EPSGT(val1, val2, (set)->epsilon) )
#define SCIPsetIsGE(set, val1, val2)       ( EPSGE(val1, val2, (set)->epsilon) )
#define SCIPsetIsZero(set, val)            ( EPSZ(val, (set)->epsilon) )
#define SCIPsetIsPositive(set, val)        ( EPSP(val, (set)->epsilon) )
#define SCIPsetIsNegative(set, val)        ( EPSN(val, (set)->epsilon) )

#define SCIPsetIsSumEQ(set, val1, val2)    ( EPSEQ(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumLT(set, val1, val2)    ( EPSLT(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumLE(set, val1, val2)    ( EPSLE(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumGT(set, val1, val2)    ( EPSGT(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumGE(set, val1, val2)    ( EPSGE(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumZero(set, val)         ( EPSZ(val, (set)->sumepsilon) )
#define SCIPsetIsSumPositive(set, val)     ( EPSP(val, (set)->sumepsilon) )
#define SCIPsetIsSumNegative(set, val)     ( EPSN(val, (set)->sumepsilon) )

#define SCIPsetIsFeasEQ(set, val1, val2)   ( EPSZ(SCIPsetRelDiff(set, val1, val2), (set)->feastol) )
#define SCIPsetIsFeasLT(set, val1, val2)   ( EPSN(SCIPsetRelDiff(set, val1, val2), (set)->feastol) )
#define SCIPsetIsFeasLE(set, val1, val2)   ( !EPSP(SCIPsetRelDiff(set, val1, val2), (set)->feastol) )
#define SCIPsetIsFeasGT(set, val1, val2)   ( EPSP(SCIPsetRelDiff(set, val1, val2), (set)->feastol) )
#define SCIPsetIsFeasGE(set, val1, val2)   ( !EPSN(SCIPsetRelDiff(set, val1, val2), (set)->feastol) )
#define SCIPsetIsFeasZero(set, val)        ( EPSZ(val, (set)->feastol) )
#define SCIPsetIsFeasPositive(set, val)    ( EPSP(val, (set)->feastol) )
#define SCIPsetIsFeasNegative(set, val)    ( EPSN(val, (set)->feastol) )

#define SCIPsetIsCutViolated(set, root, act,rhs) ( root ? EPSGT(act, rhs, (set)->cutviolepsroot) \
                                                        : EPSGT(act, rhs, (set)->cutvioleps) )

#define SCIPsetIsRelEQ(set, val1, val2)    ( EPSZ(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelLT(set, val1, val2)    ( EPSN(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelLE(set, val1, val2)    ( !EPSP(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelGT(set, val1, val2)    ( EPSP(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelGE(set, val1, val2)    ( !EPSN(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )

#define SCIPsetIsSumRelEQ(set, val1, val2) ( EPSZ(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelLT(set, val1, val2) ( EPSN(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelLE(set, val1, val2) ( !EPSP(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelGT(set, val1, val2) ( EPSP(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelGE(set, val1, val2) ( !EPSN(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )

#define SCIPsetIsInfinity(set, val)        ( (val) >= (set)->infinity )
#define SCIPsetIsFeasible(set, val)        ( (val) >= -(set)->feastol )
#define SCIPsetIsIntegral(set, val)        ( EPSISINT((val), (set)->feastol) )
#define SCIPsetIsFracIntegral(set, val)    ( !EPSP(val, (set)->feastol) )

#define SCIPsetFloor(set, val)             ( EPSFLOOR((val), (set)->feastol) )
#define SCIPsetCeil(set, val)              ( EPSCEIL((val), (set)->feastol) )
#define SCIPsetFrac(set, val)              ( EPSFRAC((val), (set)->feastol) )

#endif


#define SCIPsetCaptureBufferArray(set,ptr,num)   ( SCIPbufferCapture((set)->buffer, set, (void**)(ptr), \
                                                   (int)((num)*sizeof(**(ptr)))) )
#define SCIPsetReleaseBufferArray(set,ptr)       ( SCIPbufferRelease((set)->buffer, (void**)(ptr), 0*sizeof(**(ptr))) )
#define SCIPsetCaptureBufferSize(set,ptr,size)   ( SCIPbufferCapture((set)->buffer, set, (void**)(ptr), size) )
#define SCIPsetReleaseBufferSize(set,ptr)        ( SCIPbufferRelease((set)->buffer, (void**)(ptr), 0) )


#endif
