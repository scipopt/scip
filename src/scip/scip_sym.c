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

/**@file   scip_sym.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for symmetry handler plugins
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/set.h"
#include "scip/scip_sym.h"
#include "scip/struct_mem.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_sym.h"
#include "scip/sym.h"

/** creates a symmetry handler and includes it in SCIP.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPincludeSymhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of symmetry handler */
   const char*           desc,               /**< description of symmetry handler */
   int                   priority,           /**< priority of the symmetry handler */
   int                   proppriority,       /**< priority of the symmetry handler for propagation */
   int                   sepapriority,       /**< priority of the symmetry handler for separation */
   int                   presolpriority,     /**< priority of the symmetry handler for presolving */
   int                   propfreq,           /**< frequency for calling propagator of symmetry handler */
   int                   sepafreq,           /**< frequency for calling separator of symmetry handler */
   SCIP_Bool             delayprop,          /**< should propagation be delayed, if other sym-propagators found reductions? */
   SCIP_Bool             delaysepa,          /**< should separation be delayed, if other sym-separators found reductions? */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   int                   maxprerounds,       /**< maximal number of presolving rounds the symmetry handler participates in (-1: no limit) */
   SCIP_PROPTIMING       proptiming,         /**< positions in the node solving loop where propagation method of symmetry handlers should be executed */
   SCIP_PRESOLTIMING     presoltiming,       /**< timing mask of the symmetry handler's presolving method */
   SCIP_DECL_SYMHDLRTRYADD((*symtryadd)),    /**< addition method for symmetry method handler plugins */
   SCIP_DECL_SYMHDLRCOPY ((*symcopy)),       /**< copy method of symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symfree)),       /**< destructor method of symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*syminit)),       /**< initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXIT ((*symexit)),       /**< deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRINITSOL((*syminitsol)),  /**< solving process initialization method of symmetry handler */
   SCIP_DECL_SYMHDLREXITSOL((*symexitsol)),  /**< solving process deinitialization method of symmetry handler */
   SCIP_DECL_SYMHDLRTRANS((*symtrans)),      /**< transformation method of symmetry hanlder */
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separator for LP solutions */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separator for arbitrary primal solutions */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagation method of symmetry handler */
   SCIP_DECL_SYMHDLRRESPROP((*symresprop)),  /**< propagation conflict resolving method */
   SCIP_DECL_SYMHDLRPRESOL((*sympresol)),    /**< presolving method of symmetry handler */
   SCIP_DECL_SYMHDLRPRINT((*symprint)),      /**< print method of symmetry handler */
   SCIP_SYMHDLRDATA*     symhdlrdata         /**< symmetry handler data */
   )
{
   SCIP_SYMHDLR* symhdlr;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeSymhdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether symmetry handler is already present */
   if( SCIPfindSymhdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("symmetry handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPsymhdlrCreate(&symhdlr, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, proppriority, sepapriority, presolpriority, propfreq, sepafreq,
         delayprop, delaysepa, maxbounddist, maxprerounds, proptiming, presoltiming,
         symtryadd, symcopy, symfree, syminit, symexit, syminitsol, symexitsol, symtrans,
         symsepalp, symsepasol, symprop, symresprop, sympresol, symprint, symhdlrdata) );
   SCIP_CALL( SCIPsetIncludeSymhdlr(scip->set, symhdlr) );

   return SCIP_OKAY;
}

/** creates a symmetry handler and includes it in SCIP. All non-fundamental (or optional) callbacks will be set to NULL.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeProp() instead
 */
SCIP_RETCODE SCIPincludeSymhdlrBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR**        symhdlr,            /**< pointer to stoe symmetry handler */
   const char*           name,               /**< name of symmetry handler */
   const char*           desc,               /**< description of symmetry handler */
   int                   priority,           /**< priority of the symmetry handler */
   SCIP_DECL_SYMHDLRTRYADD((*symtryadd)),    /**< addition method for symmetry method handler plugins */
   SCIP_SYMHDLRDATA*     symhdlrdata         /**< symmetry handler data */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeSymhdlrBasic", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether symmetry handler is already present */
   if( SCIPfindSymhdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("symhdlr handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPsymhdlrCreate(symhdlr, scip->set, scip->messagehdlr, scip->mem->setmem,
         name, desc, priority, 0, 0, 0, -1, -1, TRUE, TRUE, 0.0, 0, SCIP_PROPTIMING_BEFORELP,
         SCIP_PRESOLTIMING_ALWAYS, symtryadd, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL, symhdlrdata) );
   SCIP_CALL( SCIPsetIncludeSymhdlr(scip->set, *symhdlr) );

   return SCIP_OKAY;

}

/** sets all separation related callbacks/parameters of the symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRSEPALP((*symsepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_SYMHDLRSEPASOL((*symsepasol)),  /**< separate cutting planes for arbitrary primal solution */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   sepapriority,       /**< priority of the symmetry handler for separation */
   SCIP_Bool             delaysepa           /**< should separation method be delayed, if other separators found cuts? */
   )
{
   int oldsepapriority;
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(symhdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrSepa", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   oldsepapriority = SCIPsymhdlrGetSepaPriority(symhdlr);
   SCIPsymhdlrSetSepa(symhdlr, symsepalp, symsepasol, sepafreq, sepapriority, delaysepa);

   if( oldsepapriority != sepapriority )
      scip->set->symhdlrssepasorted = FALSE;

   name = SCIPsymhdlrGetName(symhdlr);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/sepafreq", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, sepafreq) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/delaysepa", name);
   SCIP_CALL( SCIPsetSetDefaultBoolParam(scip->set, paramname, delaysepa) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/sepapriority", name);
   SCIP_CALL( SCIPsetSetDefaultBoolParam(scip->set, paramname, delaysepa) );

   return SCIP_OKAY;
}

/** sets both the propagation callback and the propagation frequency of the symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr ,           /**< symmetry handler */
   SCIP_DECL_SYMHDLRPROP ((*symprop)),       /**< propagate variable domains */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   int                   proppriority,       /**< priority of the symmetry handler for propagation */
   SCIP_PROPTIMING       proptiming          /**< positions in the node solving loop where propagation should be executed */
   )
{
   int oldpriority;
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(symhdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrProp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   oldpriority = SCIPsymhdlrGetPropPriority(symhdlr);
   SCIPsymhdlrSetProp(symhdlr, symprop, propfreq, delayprop, proppriority, proptiming);

   if( oldpriority != proppriority )
      scip->set->symhdlrspropsorted = FALSE;

   name = SCIPsymhdlrGetName(symhdlr);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/propfreq", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, propfreq) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/proptiming", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, (int) proptiming) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/delayprop", name);
   SCIP_CALL( SCIPsetSetDefaultBoolParam(scip->set, paramname, delayprop) );

   return SCIP_OKAY;
}

/** sets destructor method of symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRFREE ((*symfree))        /**< destructor of symmetry handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrFree", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(symhdlr != NULL);

   SCIPsymhdlrSetFree(symhdlr, symfree);

   return SCIP_OKAY;
}

/** sets initialization method of symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRINIT ((*syminit))        /**< initialize symmetry handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrInit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(symhdlr != NULL);

   SCIPsymhdlrSetInit(symhdlr, syminit);

   return SCIP_OKAY;
}

/** sets deinitialization method of symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrExit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLREXIT ((*symexit))        /**< deinitialize symmetry handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrExit", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(symhdlr != NULL);

   SCIPsymhdlrSetExit(symhdlr, symexit);

   return SCIP_OKAY;
}

/** sets solving process initialization method of symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrInitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRINITSOL((*syminitsol))   /**< solving process initialization method of symmetry handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrInitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(symhdlr != NULL);

   SCIPsymhdlrSetInitsol(symhdlr, syminitsol);

   return SCIP_OKAY;
}

/** sets solving process deinitialization method of symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrExitsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLREXITSOL((*symexitsol))   /**< solving process deinitialization method of symmetry handler */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrExitsol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert(symhdlr != NULL);

   SCIPsymhdlrSetExitsol(symhdlr, symexitsol);

   return SCIP_OKAY;
}

/** sets presolving method of symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRPRESOL((*sympresol)),    /**< presolving method of symmetry handler */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   int                   presolpriority,     /**< priority of the symmetry handler for presolving */
   SCIP_PRESOLTIMING     presoltiming        /**< timing mask of the constraint handler's presolving method */
   )
{
   int oldpriority;
   const char* name;
   char paramname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(symhdlr != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrPresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   oldpriority = SCIPsymhdlrGetPresolPriority(symhdlr);
   SCIP_CALL( SCIPsymhdlrSetPresol(symhdlr, sympresol, maxprerounds, presolpriority, presoltiming) );

   if( oldpriority != presolpriority )
      scip->set->symhdlrspresolsorted = FALSE;

   name = SCIPsymhdlrGetName(symhdlr);

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/maxprerounds", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, maxprerounds) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "symmetries/%s/presoltiming", name);
   SCIP_CALL( SCIPsetSetDefaultIntParam(scip->set, paramname, (int) presoltiming) );

   return SCIP_OKAY;
}

/** sets propagation conflict resolving method of symmetry handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 */
SCIP_RETCODE SCIPsetSymhdlrResprop(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMHDLR*         symhdlr,            /**< symmetry handler */
   SCIP_DECL_SYMHDLRRESPROP((*symresprop))   /**< propagation conflict resolving method */
   )
{
   assert(scip != NULL);
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetSymhdlrResprop", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIPsymhdlrSetResprop(symhdlr, symresprop);

   return SCIP_OKAY;
}

/** returns the symmetry handler of the given name, or NULL if not existing */
SCIP_SYMHDLR* SCIPfindSymhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of symmetry handler */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindSymhdlr(scip->set, name);
}

/** returns the array of currently available symmetry handlers */
SCIP_SYMHDLR** SCIPgetSymhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortSymhdlrs(scip->set);

   return scip->set->symhdlrs;
}

/** returns the number of currently available symmetry handlers */
int SCIPgetNSymhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nsymhdlrs;
}

/** outputs symmetry component information to file stream via the message handler system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note If the message handler is set to a NULL pointer nothing will be printed.
 *  @note The file stream will not be flushed directly, this can be achieved by calling SCIPinfoMessage() printing a
 *        newline character.
 */
SCIP_RETCODE SCIPprintSymcomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SYMCOMP*         symcomp,            /**< symmetry component */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(scip != NULL);
   assert(symcomp != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintSymcomp", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPsymcompPrint(symcomp, scip->set, scip->messagehdlr, file) );

   return SCIP_OKAY;
}

/** returns the symmetry components */
SCIP_SYMCOMP** SCIPgetSymcomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->syminfo != NULL);

   return scip->syminfo->symcomps;
}

/** returns number of symmetry components */
int SCIPgetNSymcomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->syminfo != NULL);

   return scip->syminfo->nsymcomps;
}
