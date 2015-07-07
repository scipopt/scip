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

/**@file   Objectives.cpp
 * @brief  Objective data structure
 * @author Timo Strunk
 *
 * Data structure storing objective data
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <sstream>
#include <algorithm>

#include "scip/cons_linear.h"

#include "LiftedWeightSpaceSolver.h"
#include "Objectives.h"
#include "main.h"

using namespace std;

/** default constructor */
Objectives::Objectives()
   : nconstraints_(0)
{
}

/** default destructor */
Objectives::~Objectives()
{
   for(std::map< SCIP_VAR*, std::vector<SCIP_Real>* >::iterator it = cost_columns_.begin();
       it != cost_columns_.end();
       ++it)
   {
      delete it->second;
   }
}

/** add new objective name */
void Objectives::addObjective(
   const char*           name                /**< identifier of objective in mps file */
   )
{
   objnames_.push_back(std::string(name));
}

/** set objective coefficient corresponding to given variable and objective name */
void Objectives::addCost(
   SCIP_VAR*             var                 /**< pointer to SCIP variable */
   const char*           objname             /**< identifier of objective in mps file */
   SCIP_Real             val                 /**< cost coefficient */
   )
{
   std::vector<SCIP_Real>* column;
   int nobjs;
   int objindex;
   std::map< SCIP_VAR*, std::vector<SCIP_Real>* >::iterator findvar;

   nobjs = getNObjs();
   objindex = objIndex(objname);
   assert(objindex >= 0);
   findvar = cost_columns_.find(var);
   if(findvar == cost_columns_.end())
   {
      column = new std::vector<SCIP_Real>(nobjs, 0.);
      cost_columns_[var] = column;
   }
   else
   {
      column = findvar->second;
   }
   assert(column != NULL);
   assert(column->size() == (unsigned int)nobjs);
   (*column)[objindex] = val;
}

/** change objective function of scip instance to new weighted objective */
SCIP_RETCODE Objectives::setWeightedObjective(
   SCIP*                          scip       /**< SCIP solver */
   const std::vector<SCIP_Real>*  weight     /**< vector containing weight for every objective */
   )
{
   int nobjs;
   std::vector<SCIP_Real>* column;
   SCIP_Real newobj;

   nobjs = getNObjs();

   assert((*weight).size() == (size_t)nobjs);

   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPfreeTransform(scip) );
   }

   for( std::map<SCIP_VAR*, std::vector<SCIP_Real>* >::const_iterator it = cost_columns_.begin();
        it != cost_columns_.end();
        ++it )
   {
      newobj = 0.;
      column = it->second;
      assert( (*column).size() == (size_t)nobjs );

      /* calculate w*C with w chosen from the interior of the weight space */
      for( int j = 0; j < nobjs; ++j )
      {
         newobj +=  (*column)[j] * ((*weight)[j]);
      }
      SCIP_CALL( SCIPchgVarObj(scip, it->first, newobj) );
   }

   return SCIP_OKAY;
}

/** creates constraint of the form wCx <= b */
SCIP_RETCODE Objectives::createObjectiveConstraint(
   SCIP*                             scip,   /**< SCIP solver */
   SCIP_CONS**                       cons,   /**< pointer for storing the created constraint */
   const std::vector<SCIP_Real>*     weight, /**< coefficients of cost vectors in constraint */
   SCIP_Real                         rhs     /**< right hand side */
   )
{
   char* name = new char[32];
   int nvars = cost_columns_.size();
   SCIP_VAR** vars = new SCIP_VAR*[nvars];
   SCIP_Real* vals = new SCIP_Real[nvars];

   int i = 0;

   /* translate objective constraint coefficients to individual variable coefficients */
   for( std::map< SCIP_VAR*, std::vector<SCIP_Real>* >::const_iterator jt = cost_columns_.begin();
        jt != cost_columns_.end();
        ++jt)
   {
      assert(jt->first != NULL);

      vars[i] = jt->first;
      vals[i] = scalar_product(*weight, *(jt->second));

      ++i;
   }

   sprintf(name,"objcons%d",++nconstraints_);
   SCIP_CALL( SCIPcreateConsBasicLinear(
         scip,
         cons,
         name,
         nvars,
         vars,
         vals,
         -SCIP_DEFAULT_INFINITY,
         rhs
         ));

   delete name;
   delete vars;
   delete vals;

   return SCIP_OKAY;
}

/** calculate the vector containing the objective value of the current solution
    for every objective */
std::vector<SCIP_Real>* Objectives::calculateCost(
   SCIP*                 scip,               /**< SCIP solver */
   SCIP_Sol*             sol                 /**< SCIP solution */
   )
{
   int nobjs;
   std::vector<SCIP_Real> * result;

   nobjs = getNObjs();
   result = new std::vector<SCIP_Real>(nobjs, 0.);

   for( std::map<SCIP_VAR*, std::vector<SCIP_Real>* >::const_iterator
           it = cost_columns_.begin();
        it != cost_columns_.end(); ++it )
   {
      for( int j = 0; j < nobjs; ++j )
      {
         (*result)[j] += SCIPgetSolVal(scip, sol, it->first) * it->second->at(j);
      }
   }

   return result;
}

/** calculate the vector containing the objective value of the SCIP primal ray
    for every objective */
std::vector<SCIP_Real>* Objectives::calculateCostRay(
   SCIP*                 scip                /**< SCIP solver */
   )
{
   int                      nobjs  = getNObjs();
   std::vector<SCIP_Real> * result = new std::vector<SCIP_Real>(nobjs, 0.);

   for( std::map<SCIP_VAR*, std::vector<SCIP_Real>* >::const_iterator
           it = cost_columns_.begin();
        it != cost_columns_.end(); ++it )
   {
      if( SCIPhasPrimalRay(scip) )
      {
         /* primal ray is given explicitly */
         for( int j = 0; j < nobjs; ++j )
         {
            (*result)[j] += SCIPgetPrimalRayVal(scip, it->first) * it->second->at(j);
         }
      }
      else
      {
         SCIP_VAR* var = it->first;
         /* find implicit primal ray */
         if(
            ( SCIPvarGetObj(var) < 0 &&
              SCIPvarGetNLocksUp(var) == 0 &&
              SCIPvarGetUbGlobal(var) == SCIPinfinity(scip) ) ||
            ( SCIPvarGetObj(var) > 0 &&
              SCIPvarGetNLocksDown(var) == 0 &&
              SCIPvarGetLbGlobal(var) == - SCIPinfinity(scip) )
            )
         {
            for( int j = 0; j < nobjs; ++j )
            {
               (*result)[j] += it->second->at(j);
            }
            break;
         }
      }
   }

   return result;
}

/* getters and setters */

int Objectives::getNObjs() const
{
   return objnames_.size();
}

/** returns the list of objective identifiers from mps file */
const std::vector<std::string>* Objectives::getObjNames() const
{
   return &objnames_;
}
/** find the objective index corresponding to the given name */
int Objectives::objIndex(
   const char*           objname             /**< identifier of objective in mps file */
   )
{
   int nobjs;

   nobjs = getNObjs();
   for(int i = 0; i < nobjs;  ++i)
   {
      if (std::string(objname) == objnames_[i])
      {
         return i;
      }
   }
   return -1;
}


/* c bindings for objective functions */

extern "C" void c_addCost(
   Objectives*           objectives,         /**< objectives class */
   SCIP_VAR*             var,                /**< pointer to SCIP variable */
   const char*           objname,            /**< identifier of objective in mps file */
   SCIP_Real             val                 /**< cost coefficient */
   )
{
   objectives->addCost(var, objname, val);
}

extern "C" void c_addObjective(
   Objectives*           probdata,           /**< objectives class */
   const char*           name                /**< identifier of objective in mps file */
   )
{
   probdata->addObjective(name);
}

/** factory function for momip data object */
extern "C" Objectives* c_createProbdata()
{
  return new Objectives();
}
