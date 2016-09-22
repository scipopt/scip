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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "cmd_line_args.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_polyhedron.h"

namespace polyscip {

    class TwoDProj {
    public:
        explicit TwoDProj(const Result& res, std::size_t first, std::size_t second);
        ValueType getFirst() const {return proj_.first;}
        ValueType getSecond() const {return proj_.second;}
        bool operator<(TwoDProj other) const;
        bool dominates(double epsilon, const TwoDProj& other) const;
    private:
        std::pair<ValueType, ValueType> proj_;
    };

    class NondomTwoDProjections {
    public:
        explicit NondomTwoDProjections(double epsilon,
                                       const ResultContainer& supported,
                                       const ResultContainer& unsupported,
                                       std::size_t first,
                                       std::size_t second);
        bool finished() const;
        bool update(const Result& res, std::size_t first, std::size_t second);
        void add(const Result& res, std::size_t first, std::size_t second);

    private:
        std::map<TwoDProj, std::vector<OutcomeType>> nondom_projections_;
        std::map<TwoDProj, std::vector<OutcomeType>>::iterator current_;
    };

    class Polyscip {
    public:
        explicit Polyscip(int argc, const char *const *argv);

        ~Polyscip();

        SCIP_RETCODE readProblem();

        SCIP_RETCODE computeNondomPoints();

        bool writeResults() const {return cmd_line_args_.writeResults();};

        void writeResultsToFile() const;

        void writeFileForVertexEnumeration() const;

        void printResults(std::ostream &os = std::cout) const;

        void printStatus(std::ostream& os = std::cout) const;

        bool dominatedPointsFound() const;

    private:

        enum class PolyscipStatus {
            Unsolved, InitPhase, WeightSpacePhase, CompUnsupportedPhase, Finished, TimeLimitReached, Error
        };

        using ValPair = std::pair<ValueType, ValueType>;
        using ValPairMap = std::map<ValPair, std::vector<OutcomeType>>;

        bool filenameIsOkay(const std::string &filename);

        //bool ValPairCmp(Polyscip::ValPair p1, Polyscip::ValPair p2);

        /** Computes first non-dominated point and initializes
         * the weight space polyhedron or finds out that there is no non-dominated point
         * @return true if first non-dom point was found and weight space polyhedron initialized;
         * false otherwise
         */
        SCIP_RETCODE initWeightSpace();

        SCIP_RETCODE computeUnitWeightOutcomes();

        void deleteWeaklyNondomSupportedResults();

        /* Return true if other element exists which dominates 'it' or has objective values coinciding with 'it
         */
        bool isDominatedOrEqual(ResultContainer::const_iterator it,
                                ResultContainer::const_iterator beg,
                                ResultContainer::const_iterator end) const;


        SCIP_RETCODE setWeightedObjective(const WeightType& weight);

        SCIP_RETCODE solve();

        SCIP_STATUS separateINFORUNBD(const WeightType& weight, bool with_presolving = true);

        SCIP_RETCODE handleNonOptNonUnbdStatus(SCIP_STATUS status);

        SCIP_RETCODE handleOptimalStatus();
        SCIP_RETCODE handleOptimalStatus(const WeightType& weight,
                                         ValueType current_opt_val);


        SCIP_RETCODE handleUnboundedStatus(bool check_if_new_result=false);

        bool outcomeIsNew(const OutcomeType& outcome, bool outcome_is_bounded) const;

        Result getResult(bool outcome_is_bounded = false, SCIP_SOL *primal_sol = nullptr);

        Result getOptimalResult();

        void printObjective(std::size_t obj_no,
                            const std::vector<int>& nonzero_indices,
                            const std::vector<SCIP_Real>& nonzero_vals,
                            std::ostream& os = std::cout) const;

        bool objIsRedundant(const std::vector<int>& begin_nonzeros,
                            const std::vector< std::vector<int> >& obj_to_nonzero_indices,
                            const std::vector< std::vector<SCIP_Real> >& obj_to_nonzero_values,
                            std::size_t index) const;

        /** Computes the supported solutions/rays and corresponding non-dominated points */
        SCIP_RETCODE computeSupported();

        /** Computes the unsupported solutions and corresponding non-dominated points */
        SCIP_RETCODE computeUnsupported();

        ValPairMap getProjectedNondomPoints(std::size_t obj_1, std::size_t obj_2) const;

        SCIP_RETCODE solveWeightedTchebycheff(SCIP_VAR* new_var,
                                              const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                              const std::vector<std::vector<ValueType>>& orig_vals,
                                              const std::pair<std::size_t, std::size_t>& considered_objs,
                                              ValPairMap nondom_projected_points);


        /** create contraint: new_var  - beta_i* vals \cdot vars >= - beta_i * ref_point[i]
         */
        SCIP_CONS* createNewVarTransformCons(SCIP_VAR *new_var,
                                             const std::vector<SCIP_VAR *> &orig_vars,
                                             const std::vector<ValueType> &orig_vals,
                                             const ValueType &rhs,
                                             const ValueType &beta_i);

        /** create constraint: lhs <= c_i^T x <= rhs*
         * @param new_var
         * @param orig_vars
         * @param orig_vals
         * @param rhs
         * @param beta_i
         * @return
         */
        SCIP_CONS* createObjValCons(const std::vector<SCIP_VAR *>& vars,
                                    const std::vector<ValueType>& vals,
                                    const ValueType& lhs,
                                    const ValueType& rhs);

        SCIP_RETCODE computeNewResult(SCIP_VAR *new_var,
                                      SCIP_CONS *cons1,
                                      SCIP_CONS *cons2,
                                      const ValueType &new_lhs,
                                      const ValueType &new_rhs,
                                      ResultContainer& new_results);

        void deleteVarNameFromResult(SCIP_VAR* var, Result& res) const;

        void printSol(const SolType& sol, std::ostream& os) const;

        /** Prints given point to given output stream */
        void outputOutcome(const OutcomeType &outcome, std::ostream& os, const std::string desc ="") const;

        bool lhsLessEqualrhs(const ValPair &lhs, const ValPair &rhs) const;

        CmdLineArgs cmd_line_args_;
        PolyscipStatus polyscip_status_;
        SCIP* scip_;
        /**< objective sense of given problem */
        SCIP_Objsense obj_sense_;
        /**< number of objectives */
        std::size_t no_objs_;
        //std::vector<std::size_t> considered_objs_;
        /**< clock measuring the time needed for the entire program */
        SCIP_CLOCK* clock_total_;

        std::unique_ptr<WeightSpacePolyhedron> weight_space_poly_;
        ResultContainer supported_;
        ResultContainer unsupported_;
        ResultContainer unbounded_;

    };

}

#endif //POLYSCIP_SRC_POLYSCIP_H_INCLUDED
