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
#pragma ident "@(#) $Id: objnodesel.h,v 1.2 2003/12/08 11:51:03 bzfpfend Exp $"

/**@file   objnodesel.h
 * @brief  C++ wrapper for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJNODESEL_H__
#define __OBJNODESEL_H__


extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for node selectors */
class ObjNodesel
{
public:
   /** name of the node selector */
   const char* const scip_name_;
   
   /** description of the node selector */
   const char* const scip_desc_;
   
   /** priority of the node selector in standard mode */
   const int scip_stdpriority_;

   /** priority of the node selector in memory saving mode */
   const int scip_memsavepriority_;

   /** does node comparison sorts w.r.t. lower bound as primal criterion? */
   const Bool scip_lowestboundfirst_;

   /** default constructor */
   ObjNodesel(
      const char*   name,               /**< name of node selector */
      const char*   desc,               /**< description of node selector */
      int           stdpriority,        /**< priority of the node selector in standard mode */
      int           memsavepriority,    /**< priority of the node selector in memory saving mode */
      Bool          lowestboundfirst    /**< does node comparison sorts w.r.t. lower bound as primal criterion? */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_stdpriority_(stdpriority),
        scip_memsavepriority_(memsavepriority),
        scip_lowestboundfirst_(lowestboundfirst)
   {
   }

   /** destructor of node selector to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel             /**< the node selector itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of node selector (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel             /**< the node selector itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of node selector (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel             /**< the node selector itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** node selection method of node selector */
   virtual RETCODE scip_select(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel,            /**< the node selector itself */
      NODE**        selnode             /**< pointer to store the selected node */
      ) = 0;
   
   /** node comparison method of node selector */
   virtual int scip_comp(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel,            /**< the node selector itself */
      NODE*         node1,              /**< first node to compare */
      NODE*         node2               /**< second node to compare */
      ) = 0;
};

} /* namespace scip */


   
/** creates the node selector for the given node selector object and includes it in SCIP */
extern
RETCODE SCIPincludeObjNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjNodesel* objnodesel,        /**< node selector object */
   Bool             deleteobject        /**< should the node selector object be deleted when node selector is freed? */
   );

#endif
