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
#pragma ident "@(#) $Id: objconshdlr.h,v 1.2 2003/12/02 14:17:21 bzfpfend Exp $"

/**@file   objconshdlr.h
 * @brief  C++ wrapper for constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJCONSHDLR_H__
#define __OBJCONSHDLR_H__


#include <cassert>

extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for constraint handlers */
class ObjConshdlr
{
public:
   /** name of the constraint handler */
   const char* const scip_name_;
   
   /** description of the constraint handler */
   const char* const scip_desc_;
   
   /** default separation priority of the constraint handler */
   const int scip_sepapriority_;

   /** default enforcing priority of the constraint handler */
   const int scip_enfopriority_;

   /** default checking priority of the constraint handler */
   const int scip_checkpriority_;

   /** default separation frequency of the constraint handler */
   const int scip_sepafreq_;

   /** default propagation frequency of the constraint handler */
   const int scip_propfreq_;

   /** should the constraint handler be skipped, if no constraints are available? */
   const Bool scip_needscons_;

   /** default constructor */
   ObjConshdlr(
      const char*   name,               /**< name of constraint handler */
      const char*   desc,               /**< description of constraint handler */
      int           sepapriority,       /**< priority of the constraint handler for separation */
      int           enfopriority,       /**< priority of the constraint handler for constraint enforcing */
      int           checkpriority,      /**< priority of the constraint handler for checking infeasibility */
      int           sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
      int           propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
      Bool          needscons           /**< should the constraint handler be skipped, if no constraints are available? */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_sepapriority_(sepapriority),
        scip_enfopriority_(enfopriority),
        scip_checkpriority_(checkpriority),
        scip_sepafreq_(sepafreq),
        scip_propfreq_(propfreq),
        scip_needscons_(needscons)
   {
   }

   /** destructor of constraint handler to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr            /**< the constraint handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of constraint handler (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr            /**< the constraint handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of constraint handler (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr            /**< the constraint handler itself */
      )
   {
      return SCIP_OKAY;
   }

   /** solving start notification method of constraint handler (called when presolving was finished) */
   virtual RETCODE scip_solstart(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< final array of constraints in transformed problem */
      int           nconss              /**< final number of constraints in transformed problem */
      )
   {
      return SCIP_OKAY;
   }
   
   /** frees specific constraint data */
   virtual RETCODE scip_delete(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONSDATA**    consdata            /**< pointer to the constraint data to free */
      )
   {
      return SCIP_OKAY;
   }

   /** transforms constraint data into data belonging to the transformed problem */
   virtual RETCODE scip_trans(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         sourcecons,         /**< source constraint to transform */
      CONS**        targetcons          /**< pointer to store created target constraint */
      ) = 0;

   /** LP initialization method of constraint handler */
   virtual RETCODE scip_initlp(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss              /**< number of constraints to process */
      )
   {
      return SCIP_OKAY;
   }

   /** separation method of constraint handler */
   virtual RETCODE scip_sepa(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      RESULT*       result              /**< pointer to store the result of the separation call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** constraint enforcing method of constraint handler for LP solutions */
   virtual RETCODE scip_enfolp(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      RESULT*       result              /**< pointer to store the result of the enforcing call */
      ) = 0;

   /** constraint enforcing method of constraint handler for pseudo solutions */
   virtual RETCODE scip_enfops(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
      RESULT*       result              /**< pointer to store the result of the enforcing call */
      ) = 0;

   /** feasibility check method of constraint handler for integral solutions */
   virtual RETCODE scip_check(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      SOL*          sol,                /**< the solution to check feasibility for */
      Bool          checkintegrality,   /**< has integrality to be checked? */
      Bool          checklprows,        /**< have current LP rows to be checked? */
      RESULT*       result              /**< pointer to store the result of the feasibility checking call */
      ) = 0;

   /** domain propagation method of constraint handler */
   virtual RETCODE scip_prop(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      RESULT*       result              /**< pointer to store the result of the propagation call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** presolving method of constraint handler */
   virtual RETCODE scip_presol(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< no. of constraints to process */
      int           nrounds,            /**< no. of presolving rounds already done */
      int           nnewfixedvars,      /**< no. of variables fixed since last call to presolving method */
      int           nnewaggrvars,       /**< no. of variables aggregated since last call to presolving method */
      int           nnewchgvartypes,    /**< no. of variable type changes since last call to presolving method */
      int           nnewchgbds,         /**< no. of variable bounds tightend since last call to presolving method */
      int           nnewholes,          /**< no. of domain holes added since last call to presolving method */
      int           nnewdelconss,       /**< no. of deleted constraints since last call to presolving method */
      int           nnewupgdconss,      /**< no. of upgraded constraints since last call to presolving method */
      int           nnewchgcoefs,       /**< no. of changed coefficients since last call to presolving method */
      int           nnewchgsides,       /**< no. of changed left or right hand sides since last call to presolving method */
      int*          nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
      int*          naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
      int*          nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
      int*          nchgbds,            /**< pointer to count total number of variable bounds tightend of all presolvers */
      int*          naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
      int*          ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
      int*          nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
      int*          nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
      int*          nchgsides,          /**< pointer to count total number of changed sides of all presolvers */
      RESULT*       result              /**< pointer to store the result of the presolving call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** conflict variable resolving method of constraint handler */
   virtual RETCODE scip_rescvar(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons,               /**< the constraint that deduced the assignment of the conflict variable */
      VAR*          infervar            /**< the binary conflict variable that has to be resolved */
      )
   {
      errorMessage("conflict variable resolving method of constraint handler <%s> not implemented\n",
         SCIPconshdlrGetName(conshdlr));
      return SCIP_INVALIDCALL;
   }

   /** variable rounding lock method of constraint handler */
   virtual RETCODE scip_lock(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons,               /**< the constraint that should lock rounding of its variables */
      int           nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
      int           nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
      ) = 0;

   /** variable rounding unlock method of constraint handler */
   virtual RETCODE scip_unlock(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons,               /**< the constraint that should unlock rounding of its variables */
      int           nlockspos,          /**< no. of times, the roundings should be unlocked for the constraint */
      int           nlocksneg           /**< no. of times, the roundings should be unlocked for the constraint's negation */
      ) = 0;

   /** constraint activation notification method of constraint handler */
   virtual RETCODE scip_active(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been activated */
      )
   {
      return SCIP_OKAY;
   }

   /** constraint deactivation notification method of constraint handler */
   virtual RETCODE scip_deactive(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been deactivated */
      )
   {
      return SCIP_OKAY;
   }

   /** constraint enabling notification method of constraint handler */
   virtual RETCODE scip_enable(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been enabled */
      )
   {
      return SCIP_OKAY;
   }

   /** constraint disabling notification method of constraint handler */
   virtual RETCODE scip_disable(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been disabled */
      )
   {
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the constraint handler for the given constraint handler object and includes it in SCIP */
RETCODE SCIPincludeObjConshdlr(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjConshdlr* objconshdlr           /**< constraint handler object */
   );

#endif
