/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   main.cpp
 * @brief  main file
 * @author Sebastian Schenker
 *
 * File containing main function; command line args are processed and
 * solver is run with given command line args.
 */

#include "LiftedWeightSpaceSolver.h"
#include "PolySCIPConfig.h"
#include "tclap/CmdLine.h"

#include <cstddef> //std::size_t
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::size_t;
using std::string;

namespace check {
  bool fileDoesNotExist(const string& filename) {
    std::ifstream file(filename.data());
    return !file;
  }
}

namespace cmd {
  using namespace TCLAP;
  
  void getCmdParams(bool& beVerbose, bool& writeSols, long& timeLimit,
		    string& writePath, string& paramFile, string& probFile, 
		    int argc, char** argv) {
    /* set up command line argument parsing */
    string executableName(EXECUTABLE_NAME), 
      versionNumber = std::to_string(POLYSCIP_VERSION_MAJOR) +
      "." + std::to_string(POLYSCIP_VERSION_MINOR);
    CmdLine cmd(executableName,' ', versionNumber);
    cmd.setExceptionHandling(false); /*set internal exception handling to false*/
    SwitchArg beVerboseArg("v", "verbose", "verbose PolySCIP cmd line output ", false);
    SwitchArg writeSolsArg("w", "writeSols", "write solutions to file", false);
    ValueArg<long> timeLimitArg("t", "timeLimit", 
				"time limit in seconds for total computation time",
				false, std::numeric_limits<long>::max(), "seconds");
    ValueArg<string> writePathArg("W", "writeSolsPath", 
				  "PATH for -w; if only -w is given the default path is ./", 
				  false, "./", "PATH");
    ValueArg<string> paramFileArg("p", "paramSets", "parameter settings file for SCIP", 
				  false, "", "paramFile.set");
    UnlabeledValueArg<string> probFileArg("probFile", "problem file in MOP format", 
					  true, "", "problemFile.mop");
    cmd.add(probFileArg);
    cmd.add(writePathArg);
    cmd.add(writeSolsArg);
    cmd.add(paramFileArg);
    cmd.add(timeLimitArg);
    cmd.add(beVerboseArg);
    cmd.parse(argc, argv);

    beVerbose = beVerboseArg.getValue();
    writeSols = writeSolsArg.getValue();
    timeLimit = timeLimitArg.getValue();    
    writePath = writePathArg.getValue();
    paramFile = paramFileArg.getValue();
    probFile = probFileArg.getValue();
  }
}

int main(int argc, char* argv[]) {
  try {
    bool beVerbose, writeSols;
    long timeLimit;
    string writePath, paramFile, probFile;
    cmd::getCmdParams(beVerbose, writeSols, timeLimit,
		      writePath, paramFile, probFile, argc, argv);
    LiftedWeightSpaceSolver solver; // the actual solver
      
    /* check time limit cmd arg */
    if (timeLimit <= 0)
      throw TCLAP::ArgException("Non-positive time limit given");
    else if (timeLimit != std::numeric_limits<long>::max())
      solver.setTimeLimit(timeLimit);

    /* check parameter file cmd arg and read settings */
    if (!paramFile.empty()) {
      if (check::fileDoesNotExist(paramFile)) 
    	throw TCLAP::ArgException("Given parameter file does not exist");
      else
    	solver.readParamSettings(paramFile);
    }

    /* check problem file cmd arg; start PolySCIP */ 
    if (check::fileDoesNotExist(probFile)) {
      throw TCLAP::ArgException("Given problem file does not exist");
    }
    else {
      SCIP_CALL( solver.readProblem(probFile, beVerbose) );
      SCIP_CALL( solver.initialise(beVerbose) );
      if (solver.foundNewUnbounded()) 
	solver.printRayInfo(cout, false, true, 
			    "#Ray(s) found in initial phase:\n");
      if (solver.foundNewOptimum()) 
	solver.printSolInfo(cout, true, true, 
			    "#Point found in initial phase:\n");

      if (solver.hasUncheckedWeight())
	cout << "\n#Weight space phase:\n";
      while (solver.hasUncheckedWeight()) {
    	SCIP_CALL( solver.solveUncheckedWeight(beVerbose) );
	if (solver.foundNewUnbounded()) 
	  solver.printRayInfo(cout, true, true);
	else if (solver.foundNewOptimum())
	  solver.printSolInfo(cout, true, true);
      }
      solver.deleteWeaklyNondomPointsAndCorrespSols();
      cout << "\n#Finished. Total computation time in seconds: " 
	   << solver.getTotalDuration() << endl;
    }
    

    /* write solutions to file */
    if (writeSols) {
      size_t prefix = probFile.find_last_of("/"), //separate path/ and filename.mop
	suffix = probFile.find_last_of("."),      //separate filename and .mop
	startInd = (prefix == string::npos) ? 0 : prefix+1,
	endInd = (suffix != string::npos) ? suffix : string::npos;
      string solFileName = "solutions_" + 
	probFile.substr(startInd,endInd-startInd)+".txt";
      if (writePath.back() != '/')
	writePath.push_back('/');
      std::ofstream solfs(writePath+solFileName);
      if (solfs.is_open()) {
	solfs << "# considered objectives = [";
	for (unsigned i=0; i<solver.getNObjs(); ++i) {
	  solver.printObjName(solfs,i);
	  solfs << ((i+1 < solver.getNObjs()) ? ", " : "]\n");
	}
	string prelim = string("# Rays in objective space and ") + 
	  string("corresp. non-zero vars in feasible space\n");
	solver.printRayInfo(solfs, false, false, prelim);
	prelim = string("# Non-dominated points in objective space ") +
	  string("and corresp. non-zero vars in feasible space\n");
	solver.printSolInfo(solfs, false, false, prelim);
	solfs.close();
	cout << "#Solution file " << solFileName 
	     << " written to: " << writePath << endl;
      }
      else
	cerr << "ERROR writing solution file\n.";
    }
  } //end try {
  catch (TCLAP::ArgException &e) {
    cerr << "ERROR: " << e.error() << " " << e.argId() << endl;
  }
  catch (TCLAP::ExitException &e) {
    cout << "\n";
  }
  return 0;
}
