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
#pragma ident "@(#) $Id: objheur.h,v 1.2 2003/12/08 11:51:03 bzfpfend Exp $"

/**@file   objheur.h
 * @brief  C++ wrapper for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJHEUR_H__
#define __OBJHEUR_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for primal heuristics */
class ObjHeur
{
public:
   /** name of the primal heuristic */
   const char* const scip_name_;
   
   /** description of the primal heuristic */
   const char* const scip_desc_;
   
   /** display character of primal heuristic */
   const char scip_dispchar_;

   /** default priority of the primal heuristic */
   const int scip_priority_;

   /** frequency for calling primal heuristic */
   const int scip_freq_;

   /** call heuristic at nodes where only a pseudo solution exist? */
   const Bool scip_pseudonodes_;

   /** default constructor */
   ObjHeur(
      const char*   name,               /**< name of primal heuristic */
      const char*   desc,               /**< description of primal heuristic */
      char          dispchar,           /**< display character of primal heuristic */
      int           priority,           /**< priority of the primal heuristic */
      int           freq,               /**< frequency for calling primal heuristic */
      Bool          pseudonodes         /**< call heuristic at nodes where only a pseudo solution exist? */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_dispchar_(dispchar),
        scip_priority_(priority),
        scip_freq_(freq),
        scip_pseudonodes_(pseudonodes)
   {
   }

   /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      HEUR*         heur                /**< the primal heuristic itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of primal heuristic (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      HEUR*         heur                /**< the primal heuristic itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of primal heuristic (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      HEUR*         heur                /**< the primal heuristic itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** execution method of primal heuristic */
   virtual RETCODE scip_exec(
      SCIP*         scip,               /**< SCIP data structure */
      HEUR*         heur,               /**< the primal heuristic itself */
      RESULT*       result              /**< pointer to store the result of the heuristic call */
      ) = 0;
};

} /* namespace scip */


   
/** creates the primal heuristic for the given primal heuristic object and includes it in SCIP */
extern
RETCODE SCIPincludeObjHeur(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjHeur*   objheur,            /**< primal heuristic object */
   Bool             deleteobject        /**< should the primal heuristic object be deleted when heuristic is freed? */
   );

#endif
