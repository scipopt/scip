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

/**@file   benders_default.c
 * @brief  default Benders' decomposition algorithm
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/benders_default.h"
#include "scip/pub_benders.h"
#include "scip/bendersdefcuts.h"


#define BENDERS_NAME            "default"
#define BENDERS_DESC            "default implementation of Benders' decomposition"
#define BENDERS_PRIORITY        0
#define BENDERS_CUTLP        TRUE   /**< should Benders' cut be generated for LP solutions */
#define BENDERS_CUTPSEUDO    TRUE   /**< should Benders' cut be generated for pseudo solutions */
#define BENDERS_CUTRELAX     TRUE   /**< should Benders' cut be generated for relaxation solutions */


/*
 * Data structures
 */

/* TODO: fill in the necessary Benders' decomposition data */

/** Benders' decomposition data */
struct SCIP_BendersData
{
   SCIP**                subproblems;        /**< the Benders' decomposition subproblems */
   SCIP_HASHMAP*         mastervartosubindex;/**< hash map from the master variable to an index for the subproblemn variables */
   SCIP_HASHMAP*         subvartomastervar;  /**< hashmap from the subproblem variable to the master variable */
   SCIP_VAR***           subproblemvars;     /**< the subproblem variables corresponding to master problem variables */
   int                   nmastervars;        /**< the number of variables in the master problem */
   int                   nsubproblems;       /**< the number of subproblems */
};




/*
 * Local methods
 */

/* creates the Benders' decomposition data */
static
SCIP_RETCODE createBendersData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   SCIP_BENDERSDATA**    bendersdata,        /**< the Benders' decomposition data */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   int i;

   assert(scip != NULL);
   assert(subproblems != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, bendersdata) );

   (*bendersdata)->nsubproblems = nsubproblems;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*bendersdata)->subproblems, nsubproblems) );

   /* Copying the subproblem to the benders data. */
   for( i = 0; i < nsubproblems; i++ )
      (*bendersdata)->subproblems[i] = subproblems[i];


   return SCIP_OKAY;
}


/* Creates the variable mappings between the  */
static
SCIP_RETCODE createVariableMappings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< the Benders' decomposition structure */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   int nsubproblems;
   int nvars;
   int nsubvars;
   char* varname;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   nsubproblems = bendersdata->nsubproblems;

   /* getting the master problem variable data */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   bendersdata->nmastervars = nvars;

   /* creating the hashmaps for the mapping between the master variables and the sub variables */
   SCIP_CALL( SCIPhashmapCreate(&bendersdata->mastervartosubindex, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&bendersdata->subvartomastervar, SCIPblkmem(scip), nvars*nsubproblems) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bendersdata->subproblemvars, nsubproblems) );
   for( i = 0; i < nsubproblems; i++ )
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bendersdata->subproblemvars, nvars) );

   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* subvar;
      SCIP_Real scalar;
      SCIP_Real constant;
      char* endptr;

      /* getting the original variable for the master variable */
      /* NOTE: This retreived the original variable. It may be a bug in regards to other parts of the code.
       * The process maps the subproblem variable to the original master variable. It was original supposed to be a
       * mapping between the subproblem variables and the transformed master variable. */
      SCIP_CALL( SCIPvarGetOrigvarSum(&vars[i], &scalar, &constant) );

      /* retreiving the var name */
      varname = SCIPvarGetName(vars[i]);

      /* retreiving the subproblem variable for the given master variable */
      for( j = 0; j < nsubproblems; j++ )
      {
         /* parsing the master variable name to find the corresponding subproblem variable */
         /* NOTE: I am using SCIPparseVarName here. This is necessary if SCIP changes the variables during the reading.
          * I expect that SCIPfindVar could be a better option. */
         SCIP_CALL( SCIPparseVarName(bendersdata->subproblems[j], varname, &subvar, &endptr) );

         /* adding the subvariable to master variable mapping into the hash map */
         if( subvar != NULL )
            SCIP_CALL( SCIPhashmapInsert(bendersdata->subvartomastervar, subvar, vars[i]) );

         /* storing the subproblem variable */
         bendersdata->subproblemvars[j][i] = subvar;
      }

      /* storing the mapping of the master variable to the variable index */
      SCIP_CALL( SCIPhashmapInsert(bendersdata->mastervartosubindex, vars[i], i) );
   }

   return SCIP_OKAY;
}



/*
 * Callback methods for Benders' decomposition
 */

/* TODO: Implement all necessary Benders' decomposition methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for benders plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BENDERSCOPY(bendersCopyDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decompostion not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersCopyDefault NULL
#endif

/** destructor of Benders' decomposition to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BENDERSFREE(bendersFreeDefault)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(*benders != NULL);

   bendersdata = SCIPbendersGetData(*benders);

   assert(bendersdata != NULL);

   for( i = bendersdata->nsubproblems - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblemvars[i], bendersdata->nmastervars);
   SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblemvars, bendersdata->nsubproblems);

   /* free hash map */
   SCIPhashmapFree(&bendersdata->subvartomastervar);
   SCIPhashmapFree(&bendersdata->mastervartosubindex);

   SCIPfreeBlockMemoryArray(scip, &bendersdata->subproblems, bendersdata->nsubproblems);

   SCIPfreeBlockMemory(scip, &bendersdata);

   return SCIP_OKAY;
}


/** initialization method of Benders' decomposition (called after problem was transformed) */
#if 0
static
SCIP_DECL_BENDERSINIT(bendersInitDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitDefault NULL
#endif


/** deinitialization method of Benders' decomposition (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BENDERSEXIT(bendersExitDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitDefault NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin)
 *
 *  This function is called immediately after the auxiliary variables are created in the master problem. The callback
 *  provides the user an opportunity to add variable data to the auxiliary variables.
 */
#if 0
static
SCIP_DECL_BENDERSINITPRE(bendersInitpreDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitpreDefault NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_BENDERSEXITPRE(bendersExitpreDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitpreDefault NULL
#endif


/** solving process initialization method of Benders' decomposition (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BENDERSINITSOL(bendersInitsolDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersInitsolDefault NULL
#endif


/** solving process deinitialization method of Benders' decomposition (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BENDERSEXITSOL(bendersExitsolDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersExitsolDefault NULL
#endif


/** mapping method between the master problem variables and the subproblem variables of Benders' decomposition */
static
SCIP_DECL_BENDERSGETVAR(bendersGetvarDefault)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;
   SCIP_VAR* mappedvar;
   SCIP_VAR* origvar;
   SCIP_Real scalar;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(benders != NULL);
   assert(var != NULL);

   bendersdata = SCIPbendersGetData(benders);

   if( probnumber == -1 )
   {
      origvar = var;
      /* The variable needs to be transformed back into an original variable. If the variable is already original, then
       * this function just returns the same variable */
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* using the original variable, the master variable can be retrieved from the hash map */
      mappedvar = (SCIP_VAR*) SCIPhashmapGetImage(bendersdata->subvartomastervar, origvar);
   }
   else
   {
      int masterindex;
      origvar = var;
      /* The variable needs to be transformed back into an original variable. If the variable is already original, then
       * this function just returns the same variable */
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* we are requesting the subproblem variable for a master problem variable */
      /* The master problem variable is a transformed variable. The original variable is not required. */
      /* NOTE: Currently the original variable is being used. This may not be correct and should be the transformed variable. */
      masterindex = (int) SCIPhashmapGetImage(bendersdata->mastervartosubindex, origvar);
      mappedvar = bendersdata->subproblemvars[probnumber][masterindex];
   }

   return mappedvar;
}

/** the execution method for Benders' decomposition */
static
SCIP_DECL_BENDERSEXEC(bendersExecDefault)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** the method for creating the Benders' decomposition subproblem. This method is called during the initialisation stage
 *  (after the master problem was transformed)
 *
 *  This method must create the SCIP instance for the subproblem and add the required variables and constraints. In
 *  addition, the settings required for the solving the problem must be set here. However, some settings will be
 *  overridden by the standard solving method included in the Benders' decomposition framework. If a special solving
 *  method is desired, the user can implement the bendersSolvesubDefault callback.
 */
static
SCIP_DECL_BENDERSCREATESUB(bendersCreatesubDefault)
{  /*lint --e{715}*/
   SCIP_BENDERSDATA* bendersdata;

   assert(scip != NULL);
   assert(benders != NULL);

   bendersdata = SCIPbendersGetData(benders);
   assert(bendersdata != NULL);

   /* adding the subproblem to the Benders' decomposition structure */
   SCIP_CALL( SCIPaddBendersSubproblem(scip, benders, bendersdata->subproblems[probnumber]) );

   /* creating the variable mappings */
   SCIP_CALL( createVariableMappings(scip, benders) );

   return SCIP_OKAY;
}

/** the subproblem solving method for Benders' decomposition. In this method the subproblem is setup with the given
 *  solution and then solved.
 *  NOTE: if the bendersSolvesubDefault callback is implemented then the bendersFreesubDefault callback must be implemented */
#if 0
static
SCIP_DECL_BENDERSSOLVESUB(bendersSolvesubDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersSolvesubDefault NULL
#endif

/** the post-solve method for Benders' decomposition */
#if 0
static
SCIP_DECL_BENDERSPOSTSOLVE(bendersPostsolveDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersPostsolveDefault NULL
#endif


/** the subproblem freeing method for Benders' decomposition. This is called between subproblem solves to clear the
 *  solving data. Generally this will only require a call to SCIPfreeTransform. However, depending on the problem it
 *  could additional freeing methods.
 *  NOTE: the bendersFreesubDefault callback must be implemented if the bendersSolvesubDefault is implemented */
#if 0
static
SCIP_DECL_BENDERSFREESUB(bendersFreesubDefault)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of default Benders' decomposition not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define bendersFreesubDefault NULL
#endif




/*
 * Benders' decomposition specific interface methods
 */

/** creates the default Benders' decomposition and includes it in SCIP */
SCIP_RETCODE SCIPincludeBendersDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   )
{
   SCIP_BENDERSDATA* bendersdata;
   SCIP_BENDERS* benders;

   /* create default Benders' decomposition data */
   bendersdata = NULL;
   SCIP_CALL( createBendersData(scip, subproblems, &bendersdata, nsubproblems) );

   benders = NULL;

   /* include Benders' decomposition */
#if 0
   /* use SCIPincludeBenders() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBenders(scip, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         bendersCopyDefault, bendersFreeDefault, bendersInitDefault, bendersExitDefault, bendersInitpreDefault, bendersExitpreDefault,
         bendersInitsolDefault, bendersExitsolDefault, bendersGetvarDefault, bendersExecDefault, bendersCreatesubDefault,
         bendersSolvesubDefault, bendersPostsolveDefault, bendersFreesubDefault, bendersdata) );
#else
   /* use SCIPincludeBendersBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBendersBasic(scip, &benders, BENDERS_NAME, BENDERS_DESC, BENDERS_PRIORITY, nsubproblems,
         BENDERS_CUTLP, BENDERS_CUTPSEUDO, BENDERS_CUTRELAX, bendersGetvarDefault, bendersExecDefault, bendersCreatesubDefault,
         bendersdata) );
   assert(benders != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetBendersCopy(scip, benders, bendersCopyDefault) );
   SCIP_CALL( SCIPsetBendersFree(scip, benders, bendersFreeDefault) );
   SCIP_CALL( SCIPsetBendersInit(scip, benders, bendersInitDefault) );
   SCIP_CALL( SCIPsetBendersExit(scip, benders, bendersExitDefault) );
   SCIP_CALL( SCIPsetBendersInitpre(scip, benders, bendersInitpreDefault) );
   SCIP_CALL( SCIPsetBendersExitpre(scip, benders, bendersExitpreDefault) );
   SCIP_CALL( SCIPsetBendersInitsol(scip, benders, bendersInitsolDefault) );
   SCIP_CALL( SCIPsetBendersExitsol(scip, benders, bendersExitsolDefault) );
   SCIP_CALL( SCIPsetBendersSolvesub(scip, benders, bendersSolvesubDefault) );
   SCIP_CALL( SCIPsetBendersPostsolve(scip, benders, bendersPostsolveDefault) );
   SCIP_CALL( SCIPsetBendersFreesub(scip, benders, bendersFreesubDefault) );
#endif

   /* OPTIONAL: including the default cuts for Benders' decomposition */
#if 0
   SCIP_CALL( SCIPincludeBendersDefaultCuts(scip, benders) );
#endif

   /* add default Benders' decomposition parameters */
   /* TODO: (optional) add Benders' decomposition specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
