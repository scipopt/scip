/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_cumulative.h
 * @brief  constraint handler for cumulative constraints
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Jens Schulz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_CUMULATIVE_H__
#define __SCIP_CONS_CUMULATIVE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/* cumulative profile */
struct CumulativeProfile
{
   int*                  timepoints;         /**< time point array */
   int*                  freecapacities;     /**< array holding corresponding available capacity */
   int                   ntimepoints;        /**< current number of entries */
   int                   arraysize;          /**< current array size */
};
typedef struct CumulativeProfile CUMULATIVEPROFILE;

/** creates the handler for cumulative constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrCumulative(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a cumulative constraint */
extern
SCIP_RETCODE SCIPcreateConsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );
   
/** returns the activities of the cumulative constraint */
extern
SCIP_VAR** SCIPgetVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the number of activities of the cumulative constraint */
extern
int SCIPgetNVarsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the capacity of the cumulative constraint */
extern
int SCIPgetCapacityCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the durations of the cumulative constraint */
extern
int* SCIPgetDurationsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the demands of the cumulative constraint */
extern
int* SCIPgetDemandsCumulative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** check for the given starting time variables with their demands and durations if the cumulative conditions for the
 *  given solution is satisfied 
 */
extern
SCIP_RETCODE SCIPcheckCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_Bool*            violated,           /**< pointer to store if the cumulative condition is violated */
   SCIP_CONS*            cons,               /**< constraint which is checked */
   SCIP_Bool             printreason         /**< should the reason for the violation be printed? */
   );
   
/** propagate the given cumulative condition */
extern
SCIP_RETCODE SCIPpropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables (jobs) */
   SCIP_VAR**            vars,               /**< array of integer variable which corresponds to starting times for a job */
   int*                  durations,          /**< array containing corresponding durations */
   int*                  demands,            /**< array containing corresponding demands */
   int                   capacity,           /**< available cumulative capacity */
   SCIP_CONS*            cons,               /**< constraint which gets propagated */
   int*                  nchgbds,            /**< pointer to store the number of variable bound changes */
   SCIP_Bool*            initialized,        /**< was conflict analysis initialized */
   SCIP_Bool*            cutoff              /**< pointer to store if the cumulative condition is violated */
   );

/** resolve propagation w.r.t. the cumulative condition */
SCIP_RETCODE SCIPrespropCumulativeCondition(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of start time variables (activities) */
   SCIP_VAR**            vars,               /**< array of start time variables */
   int*                  durations,          /**< array of durations */
   int*                  demands,            /**< array of demands */
   int                   capacity,           /**< cumulative capacity */
   SCIP_VAR*             infervar,           /**< the conflict variable whose bound change has to be resolved */
   int                   inferinfo,          /**< the user information */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< the index of the bound change, representing the point of time where the change took place */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   );

/** create a new cumulative profile for the given capacity */
extern
SCIP_RETCODE SCIPprofileCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE**   profile,            /**< pointer to store the create profile */
   int                   capacity,           /**< Capacity for this profile */
   int                   maxtimepoints       /**< maximum number of time points */
   );

/** frees given profile */
extern
void SCIPprofileFree(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE**   profile             /**< pointer to the profile */
   );

/** resizes the cumulative profile array */
extern
SCIP_RETCODE SCIPprofileResize(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< cumulative profile to resize */
   int                   newminsize          /**< minimum size to ensure */
   );

/** from the given job, the core time is computed. If core is non-empty the cumulative profile will be updated otherwise
 *  nothing happens
 */
extern
void SCIPprofileInsertCore(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   SCIP_VAR*             var,                /**< integer variable which corresponds to the starting point of the job */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            core,               /**< pointer to store if the corresponds job has a core */       
   SCIP_Bool*	         fixed,              /**< pointer to store if the job is fixed due to its bounds */ 
   SCIP_Bool*            infeasible          /**< pointer to store if the job does not fit due to capacity */
   );

/** subtracts the demand from the profile during core time of the job */
extern
void SCIPprofileDeleteCore(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   SCIP_VAR*             var,                /**< integer variable which corresponds to the starting point of the job */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            core                /**< pointer to store if the corresponds job has a core, or NULL */       
   );

/** print profile to the given file stream (for debugging) */
extern
void SCIPprofilePrint(
   SCIP*                 scip,               /**< SCIP data structure */
   CUMULATIVEPROFILE*    profile,            /**< profile to output */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** return if the given time point exists in the profile and stores the position of the given time point if it exists;
 *  otherwise the position of the next smaller existing time point */
extern
SCIP_Bool SCIPprofileFindLowerBound(
   CUMULATIVEPROFILE*    profile,              /**< profile to search in */
   int                   timepoint,            /**< time point to search for */
   int*                  pos                   /**< pointer to store the position */
   );

/** inserts the given time point into the profile if it this time point does not exists yet; returns its position in the
 *  time point array */
extern
int SCIPprofileInsertTimepoint(
   CUMULATIVEPROFILE*    profile,            /**< profile to insert the time point */
   int                   timepoint           /**< time point to insert */
   );

/** updates the profile due to inserting and removing a new job */
extern
void SCIPprofileUpdate(
   CUMULATIVEPROFILE*    profile,            /**< profile to update */
   int                   starttime,          /**< time point to start */
   int                   endtime,            /**< time point to end */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            infeasible          /**< pointer to store if the update is infeasible */
   );

/** returns TRUE if the job (given by its  demand and duration) can be inserted at the given time point; otherwise FALSE */
extern
SCIP_Bool SCIPprofileIsFeasibleStart(
   CUMULATIVEPROFILE*    profile,            /**< Cumulative profile to use */
   int                   timepoint,          /**< time point to start */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< the demand of the job */
   int*                  pos                 /**< pointer to store the earliest position where the job does not fit */
   );

/** return the earliest possible starting point within the time interval [lb,ub] for a given job (given by its duration
 *  and demand) */
extern
int SCIPprofileGetEarliestFeasibleStart(
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            infeasible          /**< pointer store if the job can not be scheduled */
   );

/** return the latest possible starting point within the time interval [lb,ub] for a given job (given by its duration
 *  and demand) */
extern
int SCIPprofileGetLatestFeasibleStart(
   CUMULATIVEPROFILE*    profile,            /**< profile to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the job */
   int                   demand,             /**< demand of the job */
   SCIP_Bool*            infeasible          /**< pointer store if the job can not be scheduled */
   );

#ifdef __cplusplus
}
#endif

#endif
