/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objnodesel.h,v 1.7 2004/12/14 12:08:00 bzfpfend Exp $"

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

   /** destructor */
   virtual ~ObjNodesel()
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
   
   /** initialization method of node selector (called after problem was transformed) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel             /**< the node selector itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of node selector (called before transformed problem is freed) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel             /**< the node selector itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** node selection method of node selector
    *
    *  This method is called to select the next leaf of the branch and bound tree to be processed.
    *
    *  possible return values for *selnode:
    *  - NULL    : problem is solved, because tree is empty
    *  - non-NULL: node to be solved next
    */
   virtual RETCODE scip_select(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel,            /**< the node selector itself */
      NODE**        selnode             /**< pointer to store the selected node */
      ) = 0;
   
   /** node comparison method of node selector
    *
    *  This method is called to compare two nodes regarding their order in the node priority queue.
    *
    *  possible return values:
    *  - value < 0: node1 comes before (is better than) node2
    *  - value = 0: both nodes are equally good
    *  - value > 0: node2 comes after (is worse than) node2
    */
   virtual int scip_comp(
      SCIP*         scip,               /**< SCIP data structure */
      NODESEL*      nodesel,            /**< the node selector itself */
      NODE*         node1,              /**< first node to compare */
      NODE*         node2               /**< second node to compare */
      ) = 0;
};

} /* namespace scip */


   
/** creates the node selector for the given node selector object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyNodesel* mynodesel = new MyNodesel(...);
 *       CHECK_OKAY( SCIPincludeObjNodesel(scip, &mynodesel, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete mynodesel;    // delete nodesel AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjNodesel(scip, new MyNodesel(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyNodesel is called here
 */
extern
RETCODE SCIPincludeObjNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjNodesel* objnodesel,        /**< node selector object */
   Bool             deleteobject        /**< should the node selector object be deleted when node selector is freed? */
   );

/** returns the nodesel object of the given name, or NULL if not existing */
extern
scip::ObjNodesel* SCIPfindObjNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of node selector */
   );

/** returns the nodesel object for the given node selector */
extern
scip::ObjNodesel* SCIPgetObjNodesel(
   SCIP*            scip,               /**< SCIP data structure */
   NODESEL*         nodesel             /**< node selector */
   );

#endif
