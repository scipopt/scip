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
#pragma ident "@(#) $Id: set.h,v 1.61 2004/01/22 14:42:30 bzfpfend Exp $"

/**@file   set.h
 * @brief  internal methods for global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SET_H__
#define __SET_H__


#include "def.h"
#include "message.h"
#include "memory.h"
#include "type_set.h"
#include "type_stat.h"
#include "buffer.h"
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

#include "struct_set.h"



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
   Bool             comments,           /**< should parameter descriptions be written as comments? */
   Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
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

/** sorts pricers by priorities */
extern
void SCIPsetSortPricers(
   SET*             set                 /**< global SCIP settings */
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

/** sorts conflict handlers by priorities */
extern
void SCIPsetSortConflicthdlrs(
   SET*             set                 /**< global SCIP settings */
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

/** sorts presolvers by priorities */
extern
void SCIPsetSortPresols(
   SET*             set                 /**< global SCIP settings */
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

/** sorts separators by priorities */
extern
void SCIPsetSortSepas(
   SET*             set                 /**< global SCIP settings */
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

/** sorts heuristics by priorities */
extern
void SCIPsetSortHeurs(
   SET*             set                 /**< global SCIP settings */
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
   SET*             set                 /**< global SCIP settings */
   );

/** calls exit methods of all user callback functions */
extern
RETCODE SCIPsetExitCallbacks(
   SET*             set                 /**< global SCIP settings */
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
   Real             feastol             /**< new feasibility tolerance */
   );

/** sets LP feasibility tolerance for reduced costs */
extern
RETCODE SCIPsetSetDualfeastol(
   SET*             set,                /**< global SCIP settings */
   Real             dualfeastol         /**< new reduced costs feasibility tolerance */
   );

/** returns the maximal number of variables priced into the LP per round */
extern
int SCIPsetGetMaxpricevars(
   const SET*       set,                /**< global SCIP settings */
   Bool             root                /**< are we at the root node? */
   );

/** returns the maximal number of cuts separated per round */
extern
int SCIPsetGetMaxsepacuts(
   const SET*       set,                /**< global SCIP settings */
   Bool             root                /**< are we at the root node? */
   );



#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
extern
Real SCIPsetRelDiff(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

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

/** checks, if the first given lower bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
extern
Bool SCIPsetIsLbBetter(
   const SET*       set,                /**< global SCIP settings */
   Real             lb1,                /**< first lower bound to compare */
   Real             lb2                 /**< second lower bound to compare */
   );

/** checks, if the first given upper bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
extern
Bool SCIPsetIsUbBetter(
   const SET*       set,                /**< global SCIP settings */
   Real             ub1,                /**< first upper bound to compare */
   Real             ub2                 /**< second upper bound to compare */
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

#define SCIPsetRelDiff(set, val1, val2)    ( (val1-val2)/(MAX3(1.0,ABS(val1),ABS(val2))) )
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

#define SCIPsetIsLbBetter(set, lb1, lb2)   ( EPSGT(lb1, lb2, (set)->boundstreps) )
#define SCIPsetIsUbBetter(set, ub1, ub2)   ( EPSLT(ub1, ub2, (set)->boundstreps) )
#define SCIPsetIsCutViolated(set, root, act, rhs) ( root ? EPSGT(act, rhs, (set)->cutviolepsroot) \
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


#define SCIPsetAllocBufferArray(set,ptr,num)    ( SCIPbufferAllocMem((set)->buffer, set, (void**)(ptr), \
                                                    (int)((num)*sizeof(**(ptr)))) )
#define SCIPsetDuplicateBufferArray(set,ptr,source,num) \
                                                ( SCIPbufferDuplicateMem((set)->buffer, set, (void**)(ptr), source, \
                                                    (int)((num)*sizeof(**(ptr)))) )
#define SCIPsetReallocBufferArray(set,ptr,num)  ( SCIPbufferReallocMem((set)->buffer, set, (void**)(ptr), \
                                                    (int)((num)*sizeof(**(ptr)))) )
#define SCIPsetFreeBufferArray(set,ptr)         ( SCIPbufferFreeMem((set)->buffer, (void**)(ptr), 0*sizeof(**(ptr))) )
#define SCIPsetAllocBufferSize(set,ptr,size)    ( SCIPbufferAllocMem((set)->buffer, set, (void**)(ptr), size) )
#define SCIPsetDuplicateBufferSize(set,ptr,source,size) \
                                                ( SCIPbufferDuplicateMem((set)->buffer, set, (void**)(ptr), source, size) )
#define SCIPsetReallocBufferSize(set,ptr,size)  ( SCIPbufferReallocMem((set)->buffer, set, (void**)(ptr), size) )
#define SCIPsetFreeBufferSize(set,ptr)          ( SCIPbufferFreeMem((set)->buffer, (void**)(ptr), 0) )


#endif
