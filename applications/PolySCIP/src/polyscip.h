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

/** @brief  PolySCIP solver class
 *
 * The PolySCIP solver class.
 */

#ifndef POLYSCIP_SRC_POLYSCIP_H_INCLUDED
#define POLYSCIP_SRC_POLYSCIP_H_INCLUDED

#include <cstdlib>
#include <iostream>
#include <ostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "cmd_line_args.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_polyhedron.h"

namespace polyscip {

    class Polyscip {
    public:
        Polyscip(int argc, const char *const *argv);

        ~Polyscip();

        SCIP_RETCODE readProblem();

        SCIP_RETCODE computeNondomPoints();

        void printSupportedResults(std::ostream& os = std::cout, bool withSolution = true);

    private:

        enum class PolyscipStatus {
            Unsolved, InitPhase, WeightSpacePhase, CompUnsupportedPhase, Finished, TimeLimitReached,
        };

        bool filenameIsOkay(const std::string &filename);

        /** Computes first non-dominated point and initializes
         * the weight space polyhedron or finds out that there is no non-dominated point
         * @return true if first non-dom point was found and weight space polyhedron initialized;
         * false otherwise
         */
        SCIP_RETCODE initWeightSpace();

        SCIP_RETCODE computeUnitWeightOutcomes();

        SCIP_RETCODE setWeightedObjective(const WeightType& weight);

        SCIP_RETCODE solve();

        SCIP_STATUS separateINFORUNBD(const WeightType& weight, bool with_presolving = true);

        SCIP_RETCODE handleNonOptNonUnbdStatus(SCIP_STATUS status);

        SCIP_RETCODE handleOptimalStatus(bool check_if_new_result=false);

        SCIP_RETCODE handleUnboundedStatus(bool check_if_new_result=false);

        bool outcomeIsNew(const OutcomeType& outcome, bool outcome_is_bounded) const;

        void addResult(bool check_if_new_result, bool outcome_is_bounded = false, SCIP_SOL* primal_sol = nullptr);

        void computeNonRedundantObjectives();

        bool objIsRedundant(size_t index) const;

        /** Computes the supported solutions/rays and corresponding non-dominated points */
        SCIP_RETCODE computeSupported();

        /** Computes the unsupported solutions and corresponding non-dominated points */
        void computeUnsupported() = delete;

        void printSol(const SolType& sol, std::ostream& os);

        /** Prints given ray to given output stream */
        void printRay(const OutcomeType &ray, std::ostream &os = std::cout);

        /** Prints given point to given output stream */
        void printPoint(const OutcomeType &point, std::ostream& os);

        CmdLineArgs cmd_line_args_;
        PolyscipStatus polyscip_status_;
        SCIP *scip_;
        /**< objective sense of given problem */
        SCIP_Objsense obj_sense_;
        /**< number of objectives */
        //std::size_t no_objs_;
        std::vector<std::size_t> non_redundant_objs_;
        /**< clock measuring the time needed for the entire program */
        SCIP_CLOCK* clock_total_;

        std::unique_ptr<WeightSpacePolyhedron> weight_space_poly_;
        ResultContainer supported_;
        ResultContainer unsupported_;
        ResultContainer unbounded_;
    };

}

#endif //POLYSCIP_SRC_POLYSCIP_H_INCLUDED
