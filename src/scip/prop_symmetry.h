/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_symmetry.h
 * @ingroup PROPAGATORS
 * @brief  propagator for symmetry handling
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_SYMMETRY_H_
#define __SCIP_PROP_SYMMETRY_H_

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <symmetry/type_symmetry.h>

/** selection rules for leaders in Schreier Sims cuts */
enum SCIP_LeaderRule
{
   SCIP_LEADERRULE_FIRSTINORBIT        = 0,       /**< first var in orbit */
   SCIP_LEADERRULE_LASTINORBIT         = 1,       /**< last var in orbit */
   SCIP_LEADERRULE_MAXCONFLICTSINORBIT = 2,       /**< var with most conflicting vars in its orbit */
   SCIP_LEADERRULE_MAXCONFLICTS        = 3        /**< var with most conflicting vars in problem */
};
typedef enum SCIP_LeaderRule SCIP_LEADERRULE;

/** tie breaks for leader rule based on the leader's orbit */
enum SCIP_LeaderTiebreakRule
{
   SCIP_LEADERTIEBREAKRULE_MINORBIT            = 0,    /**< orbit of minimum size */
   SCIP_LEADERTIEBREAKRULE_MAXORBIT            = 1,    /**< orbit of maximum size */
   SCIP_LEADERTIEBREAKRULE_MAXCONFLICTSINORBIT = 2     /**< orbit with maximum number of vars in conflict with leader */
};

/** variable types for leader in Schreier Sims cuts */
enum SCIP_SchreierSimsType
{
   SCIP_SCHREIERSIMSTYPE_BINARY                 = 1,    /**< binary variables */
   SCIP_SCHREIERSIMSTYPE_INTEGER                = 2,    /**< integer variables */
   SCIP_SCHREIERSIMSTYPE_IMPLINT                = 4,    /**< implicitly integer variables */
   SCIP_SCHREIERSIMSTYPE_CONTINUOUS             = 8     /**< continuous variables */
};

typedef enum SCIP_SchreierSimsType SCIP_SCHREIERSIMSTYPE;

/** include symmetry propagator */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return currently available symmetry group information */
SCIP_EXPORT
SCIP_RETCODE SCIPgetSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   SCIP_HASHMAP**        permvarmap,         /**< pointer to store hash map of permvars (or NULL) */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix (or NULL)*/
   int***                permstrans,         /**< pointer to store permutation generators as (npermvars x nperms) matrix (or NULL)*/
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of group size (or NULL) */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether binary variables are affected */
   int**                 components,         /**< pointer to store components of symmetry group (or NULL) */
   int**                 componentbegins,    /**< pointer to store begin positions of components in components array (or NULL) */
   int**                 vartocomponent,     /**< pointer to store assignment from variable to its component (or NULL) */
   int*                  ncomponents         /**< pointer to store number of components (or NULL) */
   );

/** return whether orbital fixing is enabled */
SCIP_EXPORT
SCIP_Bool SCIPisOrbitalfixingEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return number of the symmetry group's generators */
SCIP_EXPORT
int SCIPgetSymmetryNGenerators(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
