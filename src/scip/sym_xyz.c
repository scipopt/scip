/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sym_xyz.c
 * @ingroup DEFPLUGINS_SYM
 * @brief  symmetry handler for xyz constraints
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sym_xyz.h"


/* fundamental symmetry handler properties */
#define SYM_NAME                 "xyz"
#define SYM_DESC                 "symmetry handler template"
#define SYM_PRIORITY                 0       /**< priority of try-add function*/

/* optional symmetry handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define SYM_SEPAPRIORITY         0 /**< priority of the symmetry handler for separation */
#define SYM_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define SYM_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define SYM_MAXBOUNDDIST       1.0 /**< maximal relative distance from current node's dual bound to primal bound compared
                                    *   to best node's dual bound for applying separation */


#define SYM_PROPPRIORITY         0 /**< priority of the symmetry handler for propagation */
#define SYM_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define SYM_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define SYM_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define SYM_PRESOLPRIORITY       0 /**< priority of the symmetry handler for presolving */
#define SYM_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the symmetry handler (fast, medium, or exhaustive) */
#define SYM_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */

/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** symmetry component data for components handled by xyz symmetry handler */
struct SCIP_SymCompData
{
};

/** symmetry handler data */
struct SCIP_SymhdlrData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of symmetry handler
 */

/* TODO: Implement all necessary symmetry handler methods. The methods with #if 0 ... #else #define ... are optional */

/** addition method for symmetry method handler plugins (tries to add symmetry handling method for given symmetries) */
static
SCIP_DECL_SYMHDLRTRYADD(symTryaddXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** copy method for symmetry handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_SYMHDLRCOPY(symCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symCopyXyz NULL
#endif

/** destructor of symmetry handler to free symmetry handler data (called when SCIP is exiting) */
#if 0
static
#define SCIP_DECL_SYMHDLRFREE(symFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symFreeXyz NULL
#endif

/** initialization method of symmetry handler (called after problem was transformed) */
#if 0
static
#define SCIP_DECL_SYMHDLRINIT(symInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symInitXyz NULL
#endif

/** deinitialization method of symmetry handler (called before transformed problem is freed) */
#if 0
static
#define SCIP_DECL_SYMHDLREXIT(symExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symExitXyz NULL
#endif

/** solving process initialization method of symmetry handler (called when branch and bound process is about to begin) */
#if 0
static
#define SCIP_DECL_SYMHDLRINITSOL(symInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symInitsolXyz NULL
#endif

/** solving process deinitialization method of symmetry handler (called before branch and bound process data is freed) */
#if 0
static
#define SCIP_DECL_SYMHDLREXITSOL(symExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symExitsolXyz NULL
#endif

/** LP solution separation method of symmetry handler */
#if 0
static
#define SCIP_DECL_SYMHDLRSEPALP(symSepalpXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symSepalpXyz NULL
#endif

/** arbitrary primal solution separation method of symmetry handler */
#if 0
static
#define SCIP_DECL_SYMHDLRSEPASOL(symSepasolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symSepasolXyz NULL
#endif

/** domain propagation method of symmetry handler */
#if 0
static
#define SCIP_DECL_SYMHDLRPROP(symPropXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symPropXyz NULL
#endif

/** propagation conflict resolving method of symmetry handler */
#if 0
static
#define SCIP_DECL_SYMHDLRRESPROP(symRespropXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symRespropXyz NULL
#endif

/** presolving method of symmetry handler */
#if 0
static
#define SCIP_DECL_SYMHDLRPRESOL(symPresolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz symmetry handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define symPresolXyz NULL
#endif

/*
 * constraint specific interface methods
 */

/** creates the handler for xyz symmetry handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeSymhdlrXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SYMHDLRDATA* symhdlrdata;
   SCIP_SYMHDLR* symhdlr;

   /* create xyz symmetry handler data */
   symhdlrdata = NULL;
   /* TODO: (optional) create symmetry handler specific data here */

   symhdlr = NULL;

   /* include symmetry handler */
#if 0
   /* use SCIPincludeSymhdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeSymhdlr(scip, SYM_NAME, SYM_DESC, SYM_PRIORITY,
         SYM_PROPPRIORITY, SYM_SEPAPRIORITY, SYM_PRESOLPRIORITY,
         SYM_PROPFREQ, SYM_SEPAFREQ, SYM_DELAYPROP, SYM_DELAYSEPA, SYM_MAXBOUNDDIST,
         SYM_MAXPREROUNDS, SYM_PROPTIMING, SYM_PRESOLTIMING,
         symTryaddXyz, symCopyXyz, symFreeXyz, symInitXyz, symExitXyz, symInitsolXyz,
         symExitsolXyz, symSepalpXyz, symSepasolXyz, symPropXyz, symRespropXyz,
         symPresolXyz, symhdlrdata) );
#else
   /* use SCIPincludeSymhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeSymhdlrBasic(scip, &symhdlr, SYM_NAME, SYM_DESC, SYM_PRIORITY, symTryaddXyz, symhdlrdata) );
   assert(symhdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetSymhdlrCopy(scip, symhdlr, symCopyXyz) );
   SCIP_CALL( SCIPsetSymhdlrExit(scip, symhdlr, symExitXyz) );
   SCIP_CALL( SCIPsetSymhdlrExitsol(scip, symhdlr, symExitsolXyz) );
   SCIP_CALL( SCIPsetSymhdlrFree(scip, symhdlr, symFreeXyz) );
   SCIP_CALL( SCIPsetSymhdlrInit(scip, symhdlr, symInitXyz) );
   SCIP_CALL( SCIPsetSymhdlrInitsol(scip, symhdlr, symInitsolXyz) );
   SCIP_CALL( SCIPsetSymhdlrPresol(scip, symhdlr, symPresolXyz, SYM_MAXPREROUNDS,
         SYM_PRESOLPRIORITY, SYM_PRESOLTIMING) );
   SCIP_CALL( SCIPsetSymhdlrProp(scip, symhdlr, symPropXyz, SYM_PROPFREQ, SYM_DELAYPROP,
         SYM_PROPPRIORITY, SYM_PROP_TIMING) );
   SCIP_CALL( SCIPsetSymhdlrResprop(scip, symhdlr, symRespropXyz) );
   SCIP_CALL( SCIPsetSymhdlrSepa(scip, symhdlr, symSepalpXyz, symSepasolXyz, SYM_SEPAFREQ,
         SYM_SEPAPRIORITY, SYM_DELAYSEPA) );

#endif

   /* add xyz symmetry handler parameters */
   /* TODO: (optional) add symmetry handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
