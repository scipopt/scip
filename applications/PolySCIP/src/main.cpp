/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file main.cpp
 * @brief PolySCIP main file
 * @author Sebastian Schenker
 *
 */

#include <iostream>

#include "polyscip.h"
#include "scip/def.h"
#include "tclap/ArgException.h"

using polyscip::Polyscip;
using Status = polyscip::Polyscip::PolyscipStatus; ///< abbreviation
using std::cout;

/**
 * Main function
 * @param argc Argument count
 * @param argv Argument vector
 * @return 0
 */
int main(int argc, char** argv) {
    try {
        cout << "Starting PolySCIP...\n";
        Polyscip polyscip(argc, (const char *const *) argv);
        SCIP_CALL( polyscip.readProblem() );
        SCIP_CALL( polyscip.computeNondomPoints() );
        polyscip.printStatus();
        auto poly_status = polyscip.getStatus();
        if (poly_status == Status::Finished ||
            poly_status == Status::TimeLimitReached) {
            cout << "Number of bounded results: " << polyscip.numberOfBoundedResults() << "\n";
            cout << "Number of unbounded results: " << polyscip.numberofUnboundedResults() << "\n";
            if (polyscip.writeResults())
                polyscip.writeResultsToFile();
            else
                polyscip.printResults();
        }
        assert (!polyscip.dominatedPointsFound());
    }
    catch (TCLAP::ArgException& e) {
        std::cerr << "ERROR: " << e.error() << " " << e.argId() << "\n";
    }
    catch (TCLAP::ExitException& e) {
        std::cerr << "\n";
    }
    return 0;
}