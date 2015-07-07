/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Objectives.h
 * @brief  Objective data structure
 * @author Timo Strunk
 *
 * Data structure storing objective data
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef OBJECTIVES
#define OBJECTIVES

#include "scip/scip.h"
#include "objscip/objprobdata.h"
#include "scip/def.h"
#include <vector>
#include <map>
#include <string>

/** contains information about objective functions */
class Objectives
{
 public:
   /** default constructor */
   Objectives();

   /** default destructor */
   ~Objectives();

   /** add new objective name */
   void addObjective(
      const char*        name                /**< identifier of objective in mps file */
      );

   /** set objective coefficient corresponding to given variable and objective name */
   void addCost(
      SCIP_VAR*          var                 /**< pointer to SCIP variable */
      const char*        objname             /**< identifier of objective in mps file */
      SCIP_Real          val                 /**< cost coefficient */
      );

   /** change objective function of scip instance to new weighted objective */
   SCIP_RETCODE setWeightedObjective(
      SCIP*                             scip      /**< SCIP solver */
      const std::vector<SCIP_Real>*     weight    /**< vector containing weight for every objective */
      );

   /** creates constraint of the form wCx <= b */
   SCIP_RETCODE createObjectiveConstraint(
      SCIP*                             scip,     /**< SCIP solver */
      SCIP_CONS**                       cons,     /**< pointer for storing the created constraint */
      const std::vector<SCIP_Real>*     weight,   /**< coefficients of cost vectors in constraint */
      SCIP_Real                         rhs       /**< right hand side */
      );

   /** calculate the vector containing the objective value of the current solution
       for every objective */
   std::vector<SCIP_Real>* calculateCost(
      SCIP*              scip,               /**< SCIP solver */
      SCIP_Sol*          sol                 /**< SCIP solution */
      );

   /** calculate the vector containing the objective value of the SCIP primal ray
       for every objective */
   std::vector<SCIP_Real>* calculateCostRay(
      SCIP*              scip                /**< SCIP solver */
      );

   /** returns the number of objectives */
   int getNObjs() const;

   /** returns the list of objective identifiers from mps file */
   const std::vector<std::string>* getObjNames() const;

 private:
   std::map< SCIP_VAR*, std::vector<SCIP_Real>* > cost_columns_; /**< map from SCIP variables to cost vectors */
   std::vector<std::string>                       objnames_;     /**< list of objective identifiers from mps file */
   int                                            nconstraints_; /**< number of created objective constraints */

   /** find the objective index corresponding to the given name */
   int objIndex(
      const char*        objname             /**< identifier of objective in mps file */
      );
};

#endif
