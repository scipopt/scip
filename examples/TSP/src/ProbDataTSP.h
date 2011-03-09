/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   ProbDataTSP.h
 * @brief  C++ problem data for TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TSPPROBDATA_H__
#define __TSPPROBDATA_H__

#include "objscip/objscip.h"
#include "GomoryHuTree.h"

namespace tsp
{

/** SCIP user problem data for TSP */
class ProbDataTSP : public scip::ObjProbData
{
   GRAPH*                graph_;             /**< graph data */

public:

   /** default constructor */
   ProbDataTSP(
      GRAPH*             g                   /**< graph data */
      )
        : graph_(g)
   {
      capture_graph(graph_);
   }

   /** destructor */
   virtual ~ProbDataTSP()
   {
   }

   /** destructor of user problem data to free original user data (called when original problem is freed)
    *
    *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_delorig(
      SCIP*              scip                /**< SCIP data structure */
      );

   GRAPH* getGraph()
   {
      return graph_;
   }

};


} /* namespace tsp */

#endif
