/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   common_branch_Lookahead.c
 * @brief  Common functions used by the lookahead branching rules
 * @author Christoph Schubert
 */
#ifndef __SCIP_COMMON_BRANCH_LOOKAHEADABBREVIATED_H__
#define __SCIP_COMMON_BRANCH_LOOKAHEADABBREVIATED_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * This enum is used to represent whether an upper bound, lower bound or both are set for a variable in the ValidDomRedData
 * container.
 */
enum BoundStatus
{
   BOUNDSTATUS_NONE        = 0,              /**< The corresponding variable has no new bound set. Has to have the value 0! */
   BOUNDSTATUS_UPPERBOUND  = 1,              /**< The corresponding variable only has a new upper bound set. */
   BOUNDSTATUS_LOWERBOUND  = 2,              /**< The corresponding variable only has a new lower bound set. */
   BOUNDSTATUS_BOTH        = 3               /**< The corresponding variable has both, an upper and a lower bound set. */
};
typedef enum BoundStatus BOUNDSTATUS;

/**
 * This struct collect the bounds, which can be used in the root problem. Therefore it contains two types of bounds:
 * - The bounds that occur when both up and down branches of a variable after a first level branch are cutoff. In this case
 *   the whole first level branch can be added as a restriction.
 * - The bounds that occur implicitly while branching on the second level. See SupposedBoundData for more information.
 */
struct ValidDomRedData
{
   SCIP_Real*            upperbounds;        /**< The current upper bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   SCIP_Real*            lowerbounds;        /**< The current lower bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_LOWERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   BOUNDSTATUS*          boundstatus;        /**< The current boundstatus for each active variable. Depending on this value
                                              *   the corresponding upperbound and lowerbound values are meaningful.*/
   SCIP_Bool*            violatedbybaselp;   /**< Indicates for each variable, whether the bound change would be violated by
                                              *   solution of the base lp. */
   int*                  boundedvars;        /**< Contains the var indices that have entries in the other arrays. This array
                                              *   may be used to only iterate over the relevant variables. */
   int                   nboundedvars;       /**< The length of the boundedvars array. */
   int                   nviolatedbybaselp;  /**< The number of variables, that have a bound change in this struct, that
                                              *   would be violated by the solution of the base lp. This is the number of
                                              *   TRUEs in the 'violatedbybaselp' array. */
};
typedef struct ValidDomRedData VALIDDOMREDDATA;

/**
 * Allocates buffer memory for the given ValidDomRedData and the contained arrays.
 */
EXTERN
SCIP_RETCODE allocValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   VALIDDOMREDDATA**     validbounddata      /**< The struct to be allocated. */
);

/**
 * Frees the buffer memory of the given ValidDomRedData and the contained arrays.
 */
void freeValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   VALIDDOMREDDATA**     validbounddata      /**< The struct that should be freed. */
);

/**
 * Executes the branching on a given variable with a given value.
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE branchOnVar(
   SCIP*                 scip                /**< SCIP data structure */,
   SCIP_VAR*             var,                /**< the variable to branch on */
   SCIP_Real             val,                /**< the value to branch on */
   SCIP_Real             bestdown,
   SCIP_Bool             bestdownvalid,
   SCIP_Real             bestup,
   SCIP_Real             bestupvalid,
   SCIP_Real             provedbound
);

/**
 * Handles the assignment of ne bounds (valid and supposed). Therefore the new bound is written directly over the old
 * bound. Analog the new bound status is written directly over the old one.
 *
 * @return TRUE, if no upper and lower bound for the given var were yet set; FALSE, else.
 */
EXTERN
SCIP_Bool addBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             branchvar,          /**< The variable for which a bound should be added. */
   SCIP_Real             newbound,           /**< The value of the new bound. */
   SCIP_Bool             keepminbound,       /**< In case there is already a bound for the variable, this Indicates
                                              *   whether the min or the max value of the new and the old bound should
                                              *   be kept. */
   BOUNDSTATUS           boundtype,          /**< The type of the new bound. Must be BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_LOWERBOUND. */
   SCIP_Real*            oldbound,           /**< Pointer to the old bound. Depending on the oldboundstatus this may contain
                                              *   no meaningful data. Also gets the new bound set */
   BOUNDSTATUS*          oldboundstatus      /**< Pointer to the old boundstatus. Also gets the new status set*/
);

/**
 * Adds the given upper bound as a valid bound to the ValidDomRedData container.
 * A valid bound comes from a cutoff on the first level.
 */
EXTERN
void addValidUpperBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             baselpsolval,       /**< the lp solution of the base node */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newupperbound,      /**< the new upper bound */
   VALIDDOMREDDATA*      validbounds         /**< the container to a add the bound to */
);

/**
 * Adds the given lower bound as a valid bound to the ValidDomRedData container.
 * A valid bound comes from a cutoff on the first level.
 */
EXTERN
void addValidLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             baselpsolval,       /**< the lp solution of the base node */
   SCIP_VAR*             branchvar,          /**< the var to assign the new bound to */
   SCIP_Real             newlowerbound,      /**< the new lower bound */
   VALIDDOMREDDATA*      validbounds         /**< the container to a add the bound to */
);

/**
 * Adds the domain reductions found throughout the execution of the branching rule.
 * Domain reductions of a variable occur if:
 * - one branch on the first level is cutoff (directly or because both branches of a second level variable were cutoff)
 * - both second level branches in the same direction for the same first level variable are cutoff
 *
 * @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 * @see VALIDDOMREDDATA
 */
EXTERN
SCIP_RETCODE addDomainReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   VALIDDOMREDDATA*      validbounds,        /**< The struct containing all bounds that should be added. */
   SCIP_Bool*            domredcutoff,
   SCIP_Bool*            domred
   );

EXTERN
const char* getStatusString(
   SCIP_RESULT           result
);

#ifdef __cplusplus
}
#endif

#endif
