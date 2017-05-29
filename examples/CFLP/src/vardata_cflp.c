/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   vardata_cflp.c
 * @brief  Variable data containing the ids of constraints in which the variable appears
 * @author Stephen J. Maher
 *
 * This file implements the handling of the variable data which is attached to each file. See SCIP_VarData and \ref CFLP_PRICER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "probdata_cflp.h"
#include "vardata_cflp.h"

/** Variable data which is attached to all variables.
 *
 *  This variable data is used to store in which constraints this variable appears. Therefore, the variable data
 *  contains the ids of constraints in which the variable is part of. Hence, that data give us a column view.
 */
struct SCIP_VarData
{
   VARPROB               prob;
   SCIP_VAR**            vars;
   int                   nvars;
};

/**@name Local methods
 *
 * @{
 */

/** create a vardata */
static
SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   VARPROB               prob,               /**< the problem where this variable exists */
   int                   nvars               /**< the number of corresponding variables,
                                                  this is the number of subproblems or 1 if the variable is a subproblem var */
   )
{
   int i;

   SCIP_CALL( SCIPallocBlockMemory(scip, vardata) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*vardata)->vars, nvars) );

   for( i = 0; i < nvars; i++ )
      (*vardata)->vars[i] = NULL;

   (*vardata)->prob = prob;
   (*vardata)->nvars = nvars;

   return SCIP_OKAY;
}

/** frees user data of variable */
static
SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
   )
{
   SCIPfreeBlockMemoryArray(scip, &(*vardata)->vars, (*vardata)->nvars);
   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** frees user data of transformed variable (called when the transformed variable is freed) */
static
SCIP_DECL_VARDELORIG(vardataDelOrig)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
}/*lint !e715*/

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** create variable data */
SCIP_RETCODE SCIPvardataCreateCFLP(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   VARPROB               prob,               /**< the problem where this variable exists */
   int                   nvars               /**< the number of corresponding variables,
                                                  this is the number of subproblems or 1 if the variable is a subproblem var */
   )
{
   SCIP_CALL( vardataCreate(scip, vardata, prob, nvars) );



   return SCIP_OKAY;
}


/** creates variable */
SCIP_RETCODE SCIPcreateVarCFLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< the lower bound */
   SCIP_Real             ub,                 /**< the upper bound */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          type                /**< the variable type */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   /* create a basic variable object */
   SCIP_CALL( SCIPcreateVarBasic(scip, var, name, lb, ub, obj, type) );
   assert(*var != NULL);

   /* set callback functions */
   SCIPvarSetDelorigData(*var, vardataDelOrig);

   return SCIP_OKAY;
}


/** adds a variable to the variable mapping array */
void SCIPvardataAddVarMapping(
   SCIP_VARDATA*         vardata,            /**< variable data */
   SCIP_VAR*             mappedvar,          /**< the variable that is mapped to this variable */
   int                   probnum             /**< the problem number that the mapped var belongs to */
   )
{
   assert(vardata != NULL);

   if( vardata->prob == SUBPROB )
      vardata->vars[0] = mappedvar;
   else
      vardata->vars[probnum] = mappedvar;
}

/** returns the mapped variable for a given problem */
SCIP_VAR* SCIPvardataGetMappedVar(
   SCIP_VARDATA*         vardata,            /**< variable data */
   int                   probnum             /**< the problem number for the required mapped var, -1 if the master prob */
   )
{
   SCIP_VAR* var;
   assert(vardata != NULL);

   if( vardata->prob == SUBPROB )
      var = vardata->vars[0];
   else
      var = vardata->vars[probnum];

   return var;
}

/** returns the mapped variable for a given problem */
VARPROB SCIPvardataGetProb(
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   return vardata->prob;
}


/**@} */
