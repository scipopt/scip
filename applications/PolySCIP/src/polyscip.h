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

/**
 * @brief  PolySCIP solver class
 * @author Sebastian Schenker
 *
 * Defines PolySCIP solver class.
 */

#ifndef POLYSCIP_SRC_POLYSCIP_H_INCLUDED
#define POLYSCIP_SRC_POLYSCIP_H_INCLUDED

#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
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
        explicit TwoDProj(const OutcomeType& outcome, std::size_t first, std::size_t second);
        ValueType getFirst() const {return proj_.first;}
        ValueType getSecond() const {return proj_.second;}
        friend std::ostream &operator<<(std::ostream& os, const TwoDProj& proj);
    private:
        std::pair<ValueType, ValueType> proj_;

    };

    class NondomProjections {
    public:
        using ProjMap = std::map<TwoDProj, ResultContainer, std::function<bool(const TwoDProj&, const TwoDProj&)>>;
        explicit NondomProjections(double epsilon,
                                   const ResultContainer& supported,
                                   std::size_t first,
                                   std::size_t second);

        friend std::ostream &operator<<(std::ostream& os, const NondomProjections& nd_proj);
        bool epsilonDominates(const TwoDProj& lhs, const TwoDProj& rhs) const;
        bool finished() const;
        void update();
        void update(TwoDProj proj, Result res);
        std::vector<OutcomeType> getNondomProjOutcomes() const;
        TwoDProj getLeftProj() const {return current_->first;};
        TwoDProj getRightProj() const {return std::next(current_)->first;};
        TwoDProj getLastProj() const {return std::prev(end(nondom_projections_))->first;};


    private:
        ProjMap::iterator add(TwoDProj proj, Result res);

        double epsilon_;
        ProjMap nondom_projections_;
        ProjMap::iterator current_;
    };

    class RectangularBox {
    public:
        using Interval = std::pair<ValueType, ValueType>;
        explicit RectangularBox(const std::vector<Interval>& box);
        explicit RectangularBox(std::vector<Interval>&& box);
        friend std::ostream &operator<<(std::ostream& os, const RectangularBox& box);
        bool isSupersetOf(const RectangularBox &other) const;
        bool isSubsetOf(const RectangularBox &other) const;
        bool isDisjointFrom(const RectangularBox &other) const;
        bool isFeasible(double epsilon) const;
        std::vector<RectangularBox> getDisjointPartsFrom(double delta, const RectangularBox &other) const;
        std::size_t size() const;
        Interval getInterval(std::size_t index) const;

    private:
        Interval getIntervalIntersection(std::size_t index, const RectangularBox& other) const;

        RectangularBox(std::vector<Interval>::const_iterator first_beg, std::vector<Interval>::const_iterator first_end,
                       Interval second,
                       std::vector<Interval>::const_iterator third_beg, std::vector<Interval>::const_iterator third_end);

        constexpr static double epsilon_ = 0;
        std::vector<Interval> box_;
    };

    class Polyscip {
    public:
        enum class PolyscipStatus {
            Unsolved, // initial status after calling public constructor
            ProblemRead, // status after problem instance was read successfully
            LexOptPhase, // status after lexicographic optimal results were computed
            WeightSpacePhase, // status while results of weight space polyhedron are computed
            TwoProjPhase, // status while computing 2-projection non-dominated results
            Finished, // status if problem was solved successfully
            TimeLimitReached, // status if given time limit was reached
            Error // status if an error occured
        };

        using ObjPair = std::pair<std::size_t, std::size_t>;

        explicit Polyscip(int argc, const char *const *argv);

        ~Polyscip();

        SCIP_RETCODE readProblem();

        SCIP_RETCODE computeNondomPoints();

        bool writeResults() const {return cmd_line_args_.writeResults();};

        void writeResultsToFile() const;

        void printResults(std::ostream &os = std::cout) const;

        void printStatus(std::ostream& os = std::cout) const;

        PolyscipStatus getStatus() const;

        std::size_t numberOfBoundedResults() const;
        std::size_t numberofUnboundedResults() const;

        bool dominatedPointsFound() const;

        ResultContainer::const_iterator supportedCBegin() {return bounded_.cbegin();};
        ResultContainer::const_iterator supportedCEnd() {return bounded_.cend();};

    private:

        bool filenameIsOkay(const std::string &filename);

        // computes lexicographic optimal results
        SCIP_RETCODE computeLexicographicOptResults(std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                                    std::vector<std::vector<ValueType>>& orig_vals);

        SCIP_RETCODE computeLexicographicOptResult(std::size_t obj,
                                              std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                              std::vector<std::vector<ValueType>>& orig_vals);



        /* Return true if other element exists which epsilonDominates 'it' or has objective values coinciding with 'it
         */
        bool isDominatedOrEqual(ResultContainer::const_iterator it,
                                ResultContainer::const_iterator beg,
                                ResultContainer::const_iterator end) const;

        SCIP_RETCODE setWeightedObjective(const WeightType& weight);

        SCIP_RETCODE solve();

        SCIP_STATUS separateINFORUNBD(const WeightType& weight, bool with_presolving = true);

        SCIP_RETCODE handleNonOptNonUnbdStatus(SCIP_STATUS status);

        //SCIP_RETCODE handleOptimalStatus();
        SCIP_RETCODE handleOptimalStatus(const WeightType& weight,
                                         ValueType current_opt_val);


        SCIP_RETCODE handleUnboundedStatus(bool check_if_new_result=false);

        static bool outcomesCoincide(const OutcomeType& a, const OutcomeType& b, double epsilon);

        bool outcomeIsNew(const OutcomeType& outcome, bool outcome_is_bounded) const;

        bool outcomeIsNew(const OutcomeType& outcome,
                          ResultContainer::const_iterator beg,
                          ResultContainer::const_iterator end) const;

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
        SCIP_RETCODE computeWeightSpaceResults();

        /** Computes the unsupported solutions and corresponding non-dominated points */
        SCIP_RETCODE computeTwoProjResults(const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                           const std::vector<std::vector<ValueType>>& orig_vals);

        // return proj nondominated outcomes
        std::vector<OutcomeType> solveWeightedTchebycheff(const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                              const std::vector<std::vector<ValueType>>& orig_vals,
                                              std::size_t obj_1,
                                              std::size_t obj_2);

        std::vector<RectangularBox> computeDisjointBoxes(std::list<RectangularBox>&& feasible_boxes) const;

        std::list<RectangularBox> computeFeasibleBoxes(
                const std::map<ObjPair, std::vector<OutcomeType>> &proj_nondom_outcomes,
                const std::vector<std::vector<SCIP_VAR *>> &orig_vars,
                const std::vector<std::vector<ValueType>> &orig_vals);

        ResultContainer computeNondomPointsInBox(const RectangularBox& box,
                                                 const std::vector<std::vector<SCIP_VAR *>>& orig_vars,
                                                 const std::vector<std::vector<ValueType>>& orig_vals);

        bool boxResultIsDominated(const OutcomeType& outcome,
                                  const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                  const std::vector<std::vector<ValueType>>& orig_vals);

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

        SCIP_RETCODE computeNondomProjResult(SCIP_CONS* obj_val_cons1,
                                             SCIP_CONS* obj_val_cons2,
                                             ValueType obj_val_cons1_rhs,
                                             ValueType obj_val_cons2_rhs,
                                             std::size_t obj_1,
                                             std::size_t obj_2,
                                             ResultContainer &results);

        bool unboundedResultsExist() const {return !unbounded_.empty();};

        void printSol(const SolType& sol, std::ostream& os) const;

        /** Prints given point to given output stream */
        void outputOutcome(const OutcomeType &outcome, std::ostream& os, const std::string desc ="") const;

        /** Constructor for 4-objective case
        */
        explicit Polyscip(const CmdLineArgs& cmd_line_args,
                          SCIP* scip,
                          SCIP_Objsense obj_sense,
                          std::pair<std::size_t, std::size_t> objs_to_be_ignored,
                          SCIP_CLOCK* clock_total) = delete;

        explicit Polyscip(const CmdLineArgs& cmd_line_args,
                          SCIP *scip,
                          std::size_t no_objs,
                          SCIP_CLOCK *clock_total);

        CmdLineArgs cmd_line_args_;
        PolyscipStatus polyscip_status_;
        SCIP* scip_;
        /**< objective sense of given problem */
        SCIP_Objsense obj_sense_;
        /**< number of objectives */
        std::size_t no_objs_;
        /**< clock measuring the time needed for the entire program */
        SCIP_CLOCK* clock_total_;

        bool only_weight_space_phase_;
        bool is_lower_dim_prob_;
        bool is_sub_prob_;

        std::unique_ptr<WeightSpacePolyhedron> weight_space_poly_;
        ResultContainer bounded_;
        ResultContainer unbounded_;
    };

}

#endif //POLYSCIP_SRC_POLYSCIP_H_INCLUDED
