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
#pragma ident "@(#) $Id: set.h,v 1.84 2005/02/28 13:26:23 bzfpfend Exp $"

/**@file   set.h
 * @brief  internal methods for global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SET_H__
#define __SET_H__


#include "scip/def.h"
#include "scip/message.h"
#include "scip/memory.h"
#include "scip/buffer.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
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

#include "scip/struct_set.h"

#ifdef NDEBUG
#include "scip/pub_misc.h"
#endif




/** creates global SCIP settings */
extern
RETCODE SCIPsetCreate(
   SET**            set,                /**< pointer to SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** frees global SCIP settings */
extern
RETCODE SCIPsetFree(
   SET**            set,                /**< pointer to SCIP settings */
   BLKMEM*          blkmem              /**< block memory */
   );

/** creates a Bool parameter, sets it to its default value, and adds it to the parameter set */
extern
RETCODE SCIPsetAddBoolParam(
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
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
   BLKMEM*          blkmem,             /**< block memory */
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< name of presolver */
   );

/** sorts presolvers by priorities */
extern
void SCIPsetSortPresols(
   SET*             set                 /**< global SCIP settings */
   );

/** inserts relaxator in relaxator list */
extern
RETCODE SCIPsetIncludeRelax(
   SET*             set,                /**< global SCIP settings */
   RELAX*           relax               /**< relaxator */
   );

/** returns the relaxator of the given name, or NULL if not existing */
extern
RELAX* SCIPsetFindRelax(
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< name of relaxator */
   );

/** sorts relaxators by priorities */
extern
void SCIPsetSortRelaxs(
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
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< name of separator */
   );

/** sorts separators by priorities */
extern
void SCIPsetSortSepas(
   SET*             set                 /**< global SCIP settings */
   );

/** inserts propagator in propagator list */
extern
RETCODE SCIPsetIncludeProp(
   SET*             set,                /**< global SCIP settings */
   PROP*            prop                /**< propagator */
   );

/** returns the propagator of the given name, or NULL if not existing */
extern
PROP* SCIPsetFindProp(
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< name of propagator */
   );

/** sorts propagators by priorities */
extern
void SCIPsetSortProps(
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   );

/** returns node selector with highest priority in the current mode */
extern
NODESEL* SCIPsetGetNodesel(
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
   SET*             set,                /**< global SCIP settings */
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
   SET*             set,                /**< global SCIP settings */
   const char*      name                /**< name of event handler */
   );

/** calls init methods of all plugins */
extern
RETCODE SCIPsetInitPlugins(
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** calls exit methods of all plugins */
extern
RETCODE SCIPsetExitPlugins(
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** calls initpre methods of all plugins */
extern
RETCODE SCIPsetInitprePlugins(
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
   STAT*            stat,               /**< dynamic problem statistics */
   Bool*            unbounded,          /**< pointer to store TRUE, if presolving detected unboundness */
   Bool*            infeasible          /**< pointer to store TRUE, if presolving detected infeasibility */
   );

/** calls exitpre methods of all plugins */
extern
RETCODE SCIPsetExitprePlugins(
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
   STAT*            stat,               /**< dynamic problem statistics */
   Bool*            unbounded,          /**< pointer to store TRUE, if presolving detected unboundness */
   Bool*            infeasible          /**< pointer to store TRUE, if presolving detected infeasibility */
   );

/** calls initsol methods of all plugins */
extern
RETCODE SCIPsetInitsolPlugins(
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** calls exitsol methods of all plugins */
extern
RETCODE SCIPsetExitsolPlugins(
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** calculate memory size for dynamically allocated arrays */
extern
int SCIPsetCalcMemGrowSize(
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** calculate memory size for tree array */
extern
int SCIPsetCalcTreeGrowSize(
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** calculate memory size for path array */
extern
int SCIPsetCalcPathGrowSize(
   SET*             set,                /**< global SCIP settings */
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

/** sets LP convergence tolerance used in barrier algorithm */
extern
RETCODE SCIPsetSetBarrierconvtol(
   SET*             set,                /**< global SCIP settings */
   Real             barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   );

/** returns the maximal number of variables priced into the LP per round */
extern
int SCIPsetGetPriceMaxvars(
   SET*             set,                /**< global SCIP settings */
   Bool             root                /**< are we at the root node? */
   );

/** returns the maximal number of cuts separated per round */
extern
int SCIPsetGetSepaMaxcuts(
   SET*             set,                /**< global SCIP settings */
   Bool             root                /**< are we at the root node? */
   );



#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns value treated as infinity */
extern
Real SCIPsetInfinity(
   SET*             set                 /**< global SCIP settings */
   );

/** returns value treated as zero */
extern
Real SCIPsetEpsilon(
   SET*             set                 /**< global SCIP settings */
   );

/** returns value treated as zero for sums of floating point values */
extern
Real SCIPsetSumepsilon(
   SET*             set                 /**< global SCIP settings */
   );

/** returns feasibility tolerance for constraints */
extern
Real SCIPsetFeastol(
   SET*             set                 /**< global SCIP settings */
   );

/** returns feasibility tolerance for reduced costs */
extern
Real SCIPsetDualfeastol(
   SET*             set                 /**< global SCIP settings */
   );

/** returns convergence tolerance used in barrier algorithm */
extern
Real SCIPsetBarrierconvtol(
   SET*             set                 /**< global SCIP settings */
   );

/** returns minimal variable distance value to use for pseudo cost updates */
extern
Real SCIPsetPseudocosteps(
   SET*             set                 /**< global SCIP settings */
   );

/** returns minimal minimal objective distance value to use for pseudo cost updates */
extern
Real SCIPsetPseudocostdelta(
   SET*             set                 /**< global SCIP settings */
   );

/** checks, if values are in range of epsilon */
extern
Bool SCIPsetIsEQ(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) lower than val2 */
extern
Bool SCIPsetIsLT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) greater than val2 */
extern
Bool SCIPsetIsLE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) greater than val2 */
extern
Bool SCIPsetIsGT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) lower than val2 */
extern
Bool SCIPsetIsGE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is (positive) infinite */
extern
Bool SCIPsetIsInfinity(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   );

/** checks, if value is in range epsilon of 0.0 */
extern
Bool SCIPsetIsZero(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than epsilon */
extern
Bool SCIPsetIsPositive(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -epsilon */
extern
Bool SCIPsetIsNegative(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is integral within epsilon */
extern
Bool SCIPsetIsIntegral(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
extern
Bool SCIPsetIsScalingIntegral(
   SET*             set,                /**< global SCIP settings */
   Real             val,                /**< unscaled value to check for scaled integrality */
   Real             scalar              /**< value to scale val with for checking for integrality */
   );

/** checks, if given fractional part is smaller than epsilon */
extern
Bool SCIPsetIsFracIntegral(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value + feasibility tolerance down to the next integer in epsilon tolerance */
extern
Real SCIPsetFloor(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value - feasibility tolerance up to the next integer in epsilon tolerance */
extern
Real SCIPsetCeil(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
extern
Real SCIPsetFrac(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to return fractional part for */
   );

/** checks, if values are in range of sumepsilon */
extern
Bool SCIPsetIsSumEQ(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) lower than val2 */
extern
Bool SCIPsetIsSumLT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
extern
Bool SCIPsetIsSumLE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) greater than val2 */
extern
Bool SCIPsetIsSumGT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
extern
Bool SCIPsetIsSumGE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range sumepsilon of 0.0 */
extern
Bool SCIPsetIsSumZero(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than sumepsilon */
extern
Bool SCIPsetIsSumPositive(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -sumepsilon */
extern
Bool SCIPsetIsSumNegative(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if values are in range of feasibility tolerance */
extern
Bool SCIPsetIsFeasEQ(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPsetIsFeasLT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPsetIsFeasLE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPsetIsFeasGT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPsetIsFeasGE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range feasibility tolerance of 0.0 */
extern
Bool SCIPsetIsFeasZero(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than feasibility tolerance */
extern
Bool SCIPsetIsFeasPositive(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -feasibility tolerance */
extern
Bool SCIPsetIsFeasNegative(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is integral within the LP feasibility bounds */
extern
Bool SCIPsetIsFeasIntegral(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if given fractional part is smaller than feastol */
extern
Bool SCIPsetIsFeasFracIntegral(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value + feasibility tolerance down to the next integer in feasibility tolerance */
extern
Real SCIPsetFeasFloor(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value - feasibility tolerance up to the next integer in feasibility tolerance */
extern
Real SCIPsetFeasCeil(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. x - floor(x) in feasibility tolerance */
extern
Real SCIPsetFeasFrac(
   SET*             set,                /**< global SCIP settings */
   Real             val                 /**< value to return fractional part for */
   );

/** checks, if the first given lower bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
extern
Bool SCIPsetIsLbBetter(
   SET*             set,                /**< global SCIP settings */
   Real             lb1,                /**< first lower bound to compare */
   Real             lb2                 /**< second lower bound to compare */
   );

/** checks, if the first given upper bound is tighter (w.r.t. bound strengthening epsilon) than the second one */
extern
Bool SCIPsetIsUbBetter(
   SET*             set,                /**< global SCIP settings */
   Real             ub1,                /**< first upper bound to compare */
   Real             ub2                 /**< second upper bound to compare */
   );

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy */
extern
Bool SCIPsetIsEfficacious(
   SET*             set,                /**< global SCIP settings */
   Bool             root,               /**< should the root's minimal cut efficacy be used? */
   Real             efficacy            /**< efficacy of the cut */
   );

/** checks, if relative difference of values is in range of epsilon */
extern
Bool SCIPsetIsRelEQ(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than epsilon */
extern
Bool SCIPsetIsRelLT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
extern
Bool SCIPsetIsRelLE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than epsilon */
extern
Bool SCIPsetIsRelGT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
extern
Bool SCIPsetIsRelGE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of values is in range of sumepsilon */
extern
Bool SCIPsetIsSumRelEQ(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
extern
Bool SCIPsetIsSumRelLT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
extern
Bool SCIPsetIsSumRelLE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
extern
Bool SCIPsetIsSumRelGT(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
extern
Bool SCIPsetIsSumRelGE(
   SET*             set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsetInfinity(set)               ( (set)->num_infinity )
#define SCIPsetEpsilon(set)                ( (set)->num_epsilon )
#define SCIPsetSumepsilon(set)             ( (set)->num_sumepsilon )
#define SCIPsetFeastol(set)                ( (set)->num_feastol )
#define SCIPsetDualfeastol(set)            ( (set)->num_dualfeastol )
#define SCIPsetBarrierconvtol(set)         ( (set)->num_barrierconvtol )
#define SCIPsetPseudocosteps(set)          ( (set)->num_pseudocosteps )
#define SCIPsetPseudocostdelta(set)        ( (set)->num_pseudocostdelta )
#define SCIPsetIsEQ(set, val1, val2)       ( EPSEQ(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsLT(set, val1, val2)       ( EPSLT(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsLE(set, val1, val2)       ( EPSLE(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsGT(set, val1, val2)       ( EPSGT(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsGE(set, val1, val2)       ( EPSGE(val1, val2, (set)->num_epsilon) )
#define SCIPsetIsInfinity(set, val)        ( (val) >= (set)->num_infinity )
#define SCIPsetIsZero(set, val)            ( EPSZ(val, (set)->num_epsilon) )
#define SCIPsetIsPositive(set, val)        ( EPSP(val, (set)->num_epsilon) )
#define SCIPsetIsNegative(set, val)        ( EPSN(val, (set)->num_epsilon) )
#define SCIPsetIsIntegral(set, val)        ( EPSISINT(val, (set)->num_epsilon) )
#define SCIPsetIsScalingIntegral(set, val, scalar)                      \
   ( EPSISINT((scalar)*(val), MAX(REALABS(scalar), 1.0)*(set)->num_epsilon) )
#define SCIPsetIsFracIntegral(set, val)    ( !EPSP(val, (set)->num_epsilon) )
#define SCIPsetFloor(set, val)             ( EPSFLOOR(val, (set)->num_epsilon) )
#define SCIPsetCeil(set, val)              ( EPSCEIL(val, (set)->num_epsilon) )
#define SCIPsetFrac(set, val)              ( EPSFRAC(val, (set)->num_epsilon) )

#define SCIPsetIsSumEQ(set, val1, val2)    ( EPSEQ(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumLT(set, val1, val2)    ( EPSLT(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumLE(set, val1, val2)    ( EPSLE(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumGT(set, val1, val2)    ( EPSGT(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumGE(set, val1, val2)    ( EPSGE(val1, val2, (set)->num_sumepsilon) )
#define SCIPsetIsSumZero(set, val)         ( EPSZ(val, (set)->num_sumepsilon) )
#define SCIPsetIsSumPositive(set, val)     ( EPSP(val, (set)->num_sumepsilon) )
#define SCIPsetIsSumNegative(set, val)     ( EPSN(val, (set)->num_sumepsilon) )

#define SCIPsetIsFeasEQ(set, val1, val2)   ( EPSZ(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasLT(set, val1, val2)   ( EPSN(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasLE(set, val1, val2)   ( !EPSP(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasGT(set, val1, val2)   ( EPSP(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasGE(set, val1, val2)   ( !EPSN(SCIPrelDiff(val1, val2), (set)->num_feastol) )
#define SCIPsetIsFeasZero(set, val)        ( EPSZ(val, (set)->num_feastol) )
#define SCIPsetIsFeasPositive(set, val)    ( EPSP(val, (set)->num_feastol) )
#define SCIPsetIsFeasNegative(set, val)    ( EPSN(val, (set)->num_feastol) )
#define SCIPsetIsFeasIntegral(set, val)    ( EPSISINT(val, (set)->num_feastol) )
#define SCIPsetIsFeasFracIntegral(set, val) ( !EPSP(val, (set)->num_feastol) )
#define SCIPsetFeasFloor(set, val)         ( EPSFLOOR(val, (set)->num_feastol) )
#define SCIPsetFeasCeil(set, val)          ( EPSCEIL(val, (set)->num_feastol) )
#define SCIPsetFeasFrac(set, val)          ( EPSFRAC(val, (set)->num_feastol) )

#define SCIPsetIsLbBetter(set, lb1, lb2)   ( EPSGT(lb1, lb2, (set)->num_boundstreps) )
#define SCIPsetIsUbBetter(set, ub1, ub2)   ( EPSLT(ub1, ub2, (set)->num_boundstreps) )
#define SCIPsetIsEfficacious(set, root, efficacy) \
   ( root ? EPSP(efficacy, (set)->sepa_minefficacyroot) : EPSP(efficacy, (set)->sepa_minefficacy) )

#define SCIPsetIsRelEQ(set, val1, val2)    ( EPSZ(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelLT(set, val1, val2)    ( EPSN(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelLE(set, val1, val2)    ( !EPSP(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelGT(set, val1, val2)    ( EPSP(SCIPrelDiff(val1, val2), (set)->num_epsilon) )
#define SCIPsetIsRelGE(set, val1, val2)    ( !EPSN(SCIPrelDiff(val1, val2), (set)->num_epsilon) )

#define SCIPsetIsSumRelEQ(set, val1, val2) ( EPSZ(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelLT(set, val1, val2) ( EPSN(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelLE(set, val1, val2) ( !EPSP(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelGT(set, val1, val2) ( EPSP(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )
#define SCIPsetIsSumRelGE(set, val1, val2) ( !EPSN(SCIPrelDiff(val1, val2), (set)->num_sumepsilon) )

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
