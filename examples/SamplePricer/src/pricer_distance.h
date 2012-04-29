/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pricer_distance.h
 * @brief  p-median pricer plugin
 * @author Joerg Rambau
 * @author Andreas Tuchscherer
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRICER_DISTANCE_H__
#define __SCIP_PRICER_DISTANCE_H__


/***
    Implementation of a very basic pricing plugin for SCIP
***/

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "scip/misc.h"

//#define SCIP_DEBUG

/***
    Generate a derived class of ObjPricer.
***/

#include <string>
using namespace std;
using namespace scip;


class ObjPricerDistance : public ObjPricer
{
public:
  inline ObjPricerDistance(SCIP* scip) : 
     ObjPricer(scip, "Distance_Pricer", 
       "Finds center-customer combinations with negative reduced cost.",
       0,
       TRUE) {}

  inline ~ObjPricerDistance() {}

  inline void init(const int num_points,
		   const int num_centers,
		   const vector<SCIP_CONS*>& setpart_consptr_vector,
		   const vector<SCIP_CONS*>& setpack_consptr_vector,
		   SCIP_CONS*                card_consptr)
  {
    _num_points                = num_points;
    _num_centers               = num_centers;
    _setpart_consptr_vector    = setpart_consptr_vector;
    _setpack_consptr_vector    = setpack_consptr_vector;
    _card_consptr              = card_consptr;
    _nof_iterations            = 0;

    _distances.resize(3);
    for (unsigned int i = 0; i < _distances.size(); ++i) {
      _distances[i].resize(3);
    }
    _distances[0][0] = 0;
    _distances[0][1] = 11;
    _distances[0][2] = 15;
    _distances[1][0] = 11;
    _distances[1][1] = 0;
    _distances[1][2] = 12;
    _distances[2][0] = 15;
    _distances[2][1] = 12;
    _distances[2][2] = 0;
  }


  virtual SCIP_DECL_PRICERINIT(scip_init)
  {
     cout << "Transformed Problem" << endl;
     SCIP_CALL(SCIPprintTransProblem(scip, NULL, NULL, FALSE));
     
     for (unsigned int i = 0; i < _setpart_consptr_vector.size(); ++i) {
       SCIP_CALL(SCIPgetTransformedCons(scip, _setpart_consptr_vector[i], &_setpart_consptr_vector[i]));
     }

     for (unsigned int i = 0; i < _setpack_consptr_vector.size(); ++i) {
       SCIP_CALL(SCIPgetTransformedCons(scip, _setpack_consptr_vector[i], &_setpack_consptr_vector[i]));
     }

     SCIP_CALL(SCIPgetTransformedCons(scip, _card_consptr, &_card_consptr));
          
     return SCIP_OKAY;
   }

  virtual SCIP_DECL_PRICERREDCOST(scip_redcost)
  {  
    // find variable.

#ifdef SCIP_DEBUG
    cout << "Here is scip_redcost(SCIP* scip, SCIP_PRICER* pricer, SCIP_Real lowerbound, SCIP_RESULT* result)" << endl;
    cout << "Current Solution:" << endl;
    SCIP_CALL(SCIPprintTransProblem(scip, NULL, NULL, FALSE));
    // int nofVars = SCIPgetNVars(scip);
//     vector<SCIP_VAR*> varVec(nofVars);
//     for (
//     SCIP_CALL(SCIP
//     SCIP_CALL(SCIPprintBestTransSol(scip, NULL, NULL));
    int nofRows;
    int* basisIndexVec;
    SCIP_CALL(SCIPgetLPRowsData(scip, 0, &nofRows));
    basisIndexVec = new int[nofRows];
    SCIP_CALL(SCIPgetLPBasisInd(scip, basisIndexVec));
    for (int i = 0; i < nofRows; ++i) {
      cout << "basisIndexVec[" << i << "] = " << basisIndexVec[i] << endl;
    }
    delete basisIndexVec;
    for (unsigned int i = 0; i < _var_name_vec.size(); ++i) {
      SCIP_VAR* variable = SCIPfindVar(scip, _var_name_vec[i].c_str());
      double value = SCIPgetSolVal(scip, 0, variable);
      cout << "Variable " << _var_name_vec[i] << " = " << value << endl;
    }	   
    
#endif
    *result = SCIP_SUCCESS;
    double best_reduced_cost(SCIPinfinity(scip));
    double best_primal_cost(SCIPinfinity(scip));
    vector<int> best_sol_vec;
    for (int pot_center = 0; pot_center < _num_points; ++pot_center) {
      double total_primal_cost(0);
      double total_reduced_cost(0);
      vector<int> new_sol_vec;
      new_sol_vec.push_back(pot_center);
      
#ifdef SCIP_DEBUG
      cout << "computing reduced cost in scip_redcost for pot_center " << pot_center << ": " << endl;
#endif
      
      for (int customer = 0; customer < _num_points; ++customer) {
	const double new_reduced_cost(
				      _distances[pot_center][customer]
				      - SCIPgetDualsolSetppc(scip, _setpart_consptr_vector.at(customer))
				      - SCIPgetDualsolSetppc(scip, _setpack_consptr_vector.at(pot_center)));
#ifdef SCIP_DEBUG
	cout << "_distances[pot_center][customer] = " 
		  << _distances[pot_center][customer] << endl;
	cout << "SCIPgetDualsolSetppc(scip, _setpart_consptr_vector.at(customer)) = " 
		  << SCIPgetDualsolSetppc(scip, _setpart_consptr_vector.at(customer)) << endl;
	cout << "SCIPgetDualsolSetppc(scip, _setpack_consptr_vector.at(pot_center)) = " 
		  << SCIPgetDualsolSetppc(scip, _setpack_consptr_vector.at(pot_center)) << endl;
	cout << "reduced cost for pot_center " << pot_center
		  << " and customer " << customer << " = " << new_reduced_cost << endl;
#endif

	if (SCIPisNegative(scip, new_reduced_cost)) {
	  new_sol_vec.push_back(customer);
	  total_reduced_cost += new_reduced_cost;
	  total_primal_cost += _distances[pot_center][customer];
	}
      }

      total_reduced_cost -= SCIPgetDualsolLinear(scip, _card_consptr);

#ifdef SCIP_DEBUG
      cout << "SCIPgetDualsolLinear(scip, _card_consptr) = " 
                << SCIPgetDualsolLinear(scip, _card_consptr) << endl;
      cout << "total reduced cost for pot_center " << pot_center << " = " << total_reduced_cost << endl;
      cout << "total primal  cost for pot_center " << pot_center << " = " << total_primal_cost << endl;
#endif

      if (total_reduced_cost < best_reduced_cost) {
	best_reduced_cost = total_reduced_cost;
	best_primal_cost  = total_primal_cost;
	best_sol_vec.swap(new_sol_vec);
      }
    }
    
    if (!SCIPisNegative(scip, best_reduced_cost)) {
      return SCIP_OKAY;
    }

    SCIP_VAR* var;

    char var_name[255];
    char name_part[255];
    SCIPsnprintf(var_name, 255, "%d", best_sol_vec[0]);
    strcat(var_name, "-");
    for (unsigned int i = 1; i < best_sol_vec.size(); ++i) {
      SCIPsnprintf(name_part, 255, "%d", best_sol_vec[i]);
      strcat(var_name, name_part);
    } 
    
#ifdef SCIP_DEBUG
    cout << "NEW VARIABLE = " << var_name << endl << endl;
#endif
    _var_name_vec.push_back(string(var_name));

    SCIP_CALL(SCIPcreateVar(scip,
			     &var,
			     var_name,
			     0.0,
                             SCIPinfinity(scip), // 1.0
                             best_primal_cost,
                             SCIP_VARTYPE_INTEGER, // SCIP_VARTYPE_BINARY
			     false,
                             false,
			     0,
			     0,
			     0,
			     0,
			     0));

    SCIP_CALL(SCIPaddPricedVar(scip, var, 1));

    SCIP_CALL(SCIPaddCoefSetppc(scip,
				 _setpack_consptr_vector.at(best_sol_vec[0]),
				 var));
    for (unsigned int i = 1; i < best_sol_vec.size(); ++i) {
      SCIP_CALL(SCIPaddCoefSetppc(scip,
				   _setpart_consptr_vector.at(best_sol_vec[i]),
				   var));
    }

    SCIP_CALL(SCIPaddCoefLinear(scip,
				 _card_consptr,
				 var,
				 1));

    SCIP_CALL(SCIPreleaseVar(scip, &var));

    return SCIP_OKAY;
  }
  
  virtual SCIP_DECL_PRICERFARKAS(scip_farkas)
  {

#ifdef SCIP_DEBUG
    cout << "Here is scip_farkas(SCIP* scip, SCIP_PRICER* pricer)" << endl;
#endif

    double best_primal_cost(SCIPinfinity(scip));
    vector<int> best_sol_vec;
    for (int pot_center = 0; pot_center < _num_points; ++pot_center) {
      double total_primal_cost(0);
      vector<int> new_sol_vec;
      new_sol_vec.push_back(pot_center);

#ifdef SCIP_DEBUG
      cout << "computing primal cost for one center in scip_farkas for pot_center " << pot_center 
		<< ":" << endl;
#endif

      for (int customer = 0; customer < _num_points; ++customer) {
	new_sol_vec.push_back(customer);
	total_primal_cost += _distances[pot_center][customer];

#ifdef SCIP_DEBUG
	cout << "_distances[pot_center][customer] = " 
		  << _distances[pot_center][customer] << endl;
	cout << "primal cost for pot_center " << pot_center
		  << " and customer " << customer << " = " << _distances[pot_center][customer] << endl;
#endif
      }

#ifdef SCIP_DEBUG
      cout << "total primal cost for pot_center " << pot_center << " = " << total_primal_cost << endl;
#endif

      if (total_primal_cost < best_primal_cost) {
	best_primal_cost  = total_primal_cost;
	best_sol_vec.swap(new_sol_vec);
      }
    }

    SCIP_VAR* var;

    char var_name[255];
    char name_part[255];
    SCIPsnprintf(var_name, 255, "%d", best_sol_vec[0]);
    strcat(var_name, "-");
    for (unsigned int i = 1; i < best_sol_vec.size(); ++i) {
      SCIPsnprintf(name_part, 255, "%d", best_sol_vec[i]);
      strcat(var_name, name_part);
    } 
    
    _var_name_vec.push_back(string(var_name));
    
    SCIP_CALL(SCIPcreateVar(scip,
			     &var,
			     var_name,
			     0.0,
                             SCIPinfinity(scip), // 1.0
                             best_primal_cost,
                             SCIP_VARTYPE_INTEGER, // SCIP_VARTYPE_BINARY,
			     false,
                             false,
			     0,
			     0,
			     0,
			     0,
			     0));

    SCIP_CALL(SCIPaddVar(scip, var));

    SCIP_CALL(SCIPaddCoefSetppc(scip,
				 _setpack_consptr_vector.at(best_sol_vec[0]),
				 var));
    
    for (unsigned int i = 1; i < best_sol_vec.size(); ++i) {
      SCIP_CALL(SCIPaddCoefSetppc(scip,
				   _setpart_consptr_vector.at(best_sol_vec[i]),
				   var));
    }

    SCIP_CALL(SCIPaddCoefLinear(scip,
				 _card_consptr,
				 var,
				 1));
    
    SCIP_CALL(SCIPreleaseVar(scip, &var));

    return SCIP_OKAY;    
  }
  
private:
  int                           _num_points;
  int                           _num_centers;
  vector< vector<double> > _distances;
  vector< SCIP_CONS* >          _setpart_consptr_vector;
  vector< SCIP_CONS* >          _setpack_consptr_vector;
  SCIP_CONS*                    _card_consptr;
  int                           _nof_iterations;
  vector<string>           _var_name_vec;
};


#endif

