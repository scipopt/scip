/**
 * @file main.cpp
 * @brief Solution Checker main file
 *
 * @author Domenico Salvagnin
 */

#include "model.h"
#include "mpsinput.h"
#include "gmputils.h"

#include <stdlib.h>
#include <string.h>
#include <iostream>

void usage()
{
   std::cout << "Usage: solchecker filename.mps[.gz] solution.sol [linear_tol int_tol]" << std::endl;
}

int main (int argc, char const *argv[])
{
   if( argc < 3 )
   {
      usage();
      return -1;
   }

   // read model
   Model* model = new Model;
   MpsInput* mpsi = new MpsInput;
   
   bool success = mpsi->readMps(argv[1], model);
   std::cout << "Read MPS: " << success << std::endl;
   if( !success ) return -1;
   
   std::cout << "MIP has " << model->numVars() << " vars and " << model->numConss() << " constraints" << std::endl;

   // read solution
   success = model->readSol(argv[2]);
   std::cout << "Read SOL: " << success << std::endl;
   if( !success ) return -1;
   if( model->hasObjectiveValue ) std::cout << "Objective value computed by solver: " << model->objectiveValue.toDouble() << std::endl;
   std::cout << std::endl;
   
   // read tolerances and other arguments
   Rational linearTolerance(1, 10000);
   Rational intTolerance(linearTolerance);
   std::string solverStatus("unknown");
   std::string exactStatus("unknown");
   Rational exactObjVal;
   
   int posArg = 0;
   for( int k = 3; k < argc; k++)
   {
      if( strncmp(argv[k], "exactStatus=", 12) == 0 )
      {
         exactStatus = std::string(argv[k] + 12);
         if( (exactStatus != "unknown") && (exactStatus != "infeasible") && (exactStatus != "feasible") ) exactStatus = "unknown";
      }
      else if( strncmp(argv[k], "exactObjVal=", 12) == 0 )
      {
         exactObjVal.fromString(argv[k] + 12);
      }
      else if( strncmp(argv[k], "solverStatus=", 13) == 0 )
      {
         solverStatus = std::string(argv[k] + 13);
         if( (solverStatus != "unknown") && (solverStatus != "stopped") && (solverStatus != "solved") ) solverStatus = "unknown";
      }
      else
      {
         if( posArg == 0 )
         {
            linearTolerance.fromString(argv[k]);
            intTolerance.fromString(argv[k]);
         }
         else
         {
            intTolerance.fromString(argv[k]);
         }
         posArg++;
      }
   }
   std::cout << "Integrality tolerance: " << intTolerance.toString() << std::endl;
   std::cout << "Linear tolerance: " << linearTolerance.toString() << std::endl;
   std::cout << "Objective tolerance: " << linearTolerance.toString() << std::endl;
   std::cout << std::endl;
   //std::cout << "SolverStatus: " << solverStatus << std::endl;
   //std::cout << "ExactStatus: " << exactStatus << std::endl;
   //std::cout << "ExactObjVal: " << exactObjVal.toString() << std::endl;
   //std::cout << std::endl;
   
   // check!
   bool intFeas;
   bool linFeas;
   bool obj;
   model->check(intTolerance, linearTolerance, intFeas, linFeas, obj);
   
   std::cout << "Check SOL: Integrality " << intFeas 
      << " Constraints " << linFeas
      << " Objective " << obj << std::endl;
   
   // stats on violations
   Rational intViol;
   Rational linearViol;
   Rational objViol;
   model->maxViolations(intViol, linearViol, objViol);
   std::cout << "Maximum violations: Integrality " << intViol.toDouble()
      << " Constraints " << linearViol.toDouble()
      << " Objective " << objViol.toDouble() << std::endl;
   
   // check w.r.t exact solver
   //bool exactFeas;
   //bool exactObj;
   //model->checkWrtExact(solverStatus, exactStatus, exactObjVal, linearTolerance, exactFeas, exactObj);
   //std::cout << "Check w.r.t. exact information: Feasibility " << exactFeas
   //   << " Objective " << exactObj << std::endl;
   
   // clean up
   delete mpsi;
   delete model;

   return 0;
}
