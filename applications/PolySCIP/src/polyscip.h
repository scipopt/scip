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
 * @file polyscip.h
 * @brief PolySCIP solver class
 * @author Sebastian Schenker
 *
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

    /**
     * @class TwoDProj
     * @brief Class representing a two-dimensional projection of an outcome
     */
    class TwoDProj {
    public:

        /**
         * Default constructor
         * @param outcome Corresponding outcome to take two-dimensional projection of
         * @param first First (objective) index of outcome to consider for projection
         * @param second Second (objective) index of outcome to consider for projection
         */
        explicit TwoDProj(const OutcomeType& outcome,
                          std::size_t first,
                          std::size_t second);

        /**
         * Get first projection value
         * @return First value of projection
         */
        ValueType getFirst() const {return proj_.first;}

        /**
         * Get second projection value
         * @return Second value of projection
         */
        ValueType getSecond() const {return proj_.second;}

        /**
         * Ostream operator
         * @param os Output stream
         * @param proj Projection to write to stream
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream& os, const TwoDProj& proj);

    private:
        std::pair<ValueType, ValueType> proj_; ///< Pair of projection values

    };

    /**
     * @class NondomProjections
     * @brief Class representing non-dominated projections
     */
    class NondomProjections {
    public:
        using ProjMap = std::map<TwoDProj, ResultContainer, std::function<bool(const TwoDProj&, const TwoDProj&)>>; ///< Container for non-dominated projections

        /**
         * Default constructor
         * @param epsilon Error value for comparisons
         * @param supported Results to take non-dominated projections
         * @param first First (objective) index to consider for projection
         * @param second Second (objective) index to consider for projection
         */
        explicit NondomProjections(double epsilon,
                                   const ResultContainer& supported,
                                   std::size_t first,
                                   std::size_t second);

        /**
         * Ostream operator
         * @param os Output stream
         * @param nd_proj Non-dominated projections to write to stream
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream& os, const NondomProjections& nd_proj);

        /**
         * lhs-Projection epsilonDominates rhs-Projection if lhs.first - epsilon < rhs.first && lhs.second - epsilon < rhs.second
         * @param lhs lhs-Projection
         * @param rhs rhs-Projection
         * @return true if lhs-Projection epsilon-dominated rhs-Projection; false otherwise
         */
        bool epsilonDominates(const TwoDProj& lhs,
                              const TwoDProj& rhs) const;

        /**
         * Indicates that all stored projections are investigated
         * @return true if all stored projections have been investigated; false otherwise
         */
        bool finished() const;

        /**
         * Advances current_ iterator
         */
        void update();

        /**
         * Incorporates a new projection and corresponding result into non-dominated projections
         * @param proj Projection to incorporated
         * @param res Corresponding result of projection
         */
        void update(TwoDProj proj, Result res);

        /**
         * Get outcomes corresponding to non-dominated projections
         * @return Vector of outcomes corresponding to non-dominated projections
         */
        std::vector<OutcomeType> getNondomProjOutcomes() const;

        /**
         * Get projection to be investigated
         * @return Two-dimensional projection
         */
        TwoDProj getLeftProj() const {return current_->first;};

        /**
         * Get right neighbour of projection to be investigated
         * @return Two-dimensional projection
         */
        TwoDProj getRightProj() const {return std::next(current_)->first;};

        /**
         * Get projection with maximal value for second index
         * @return Two-dimensional projection
         */
        TwoDProj getLastProj() const {return std::prev(end(nondom_projections_))->first;};


    private:
        /**
         * Add projection and corresponding result to non-dominated projections
         * @param proj Projection to add
         * @param res Corresponding result of projections
         * @return Iterator pointing to proj
         */
        ProjMap::iterator add(TwoDProj proj,
                              Result res);

        double epsilon_; ///< Epsilon value used in fct 'epsilonDominates'
        ProjMap nondom_projections_; ///< Container for non-dominated projections
        ProjMap::iterator current_; ///< Currently investigated projection
    };

    /**
     * @class RectangularBox
     * @brief A rectangular box R = [a_1,e_1) x ... x [a_k,e_k) is a k-ary Cartesian product of half-open intervals
     */
    class RectangularBox {
    public:
        using Interval = std::pair<ValueType, ValueType>; ///< Interval I = [a,b)

        /**
         * Copy constructor
         * @param box RectangularBox to copy
         */
        explicit RectangularBox(const std::vector<Interval>& box);

        /**
         * Move constructor
         * @param box RectangularBox to copy
         */
        explicit RectangularBox(std::vector<Interval>&& box);

        /**
         * Ostream operator
         * @param os Output stream to write to
         * @param box Box to write to stream
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream& os, const RectangularBox& box);

        /**
         * Indicates whether given box is subset
         * @param other Box to compare
         * @return true if 'other' is subset; false otherwise
         */
        bool isSupersetOf(const RectangularBox &other) const;

        /**
         * Indicates whether given box is superset
         * @param other Box to compare
         * @return true if 'other' is superset; false otherwise
         */
        bool isSubsetOf(const RectangularBox &other) const;

        /**
         * Indicates whether given box is disjoint
         * @param other Box to compare
         * @return true if 'other' is disjoint; false otherwise
         */
        bool isDisjointFrom(const RectangularBox &other) const;

        /**
         * Indicates whether a_i + epsilon > e_i for all i
         * @param epsilon Value to add to left interval limit
         * @return true if a_i + epsilon <= e_i for all i; false otherwise
         */
        bool isFeasible(double epsilon) const;

        /**
         * Makes disjoint rectangular boxes with respect to given box
         * @param delta Feasibility threshold
         * @param other Box to compare to
         * @return Container of disjoint boxes
         */
        std::vector<RectangularBox> getDisjointPartsFrom(double delta, const RectangularBox &other) const;

        /**
         * Get get number of intervals
         * @return Dimension of rectangular box
         */
        std::size_t size() const;

        /**
         * Get interval of box
         * @param index Corresponding interval index
         * @return Interval corresponding to index
         */
        Interval getInterval(std::size_t index) const;

        /**
         * Indicates whether outcome dominates entire box
         * @param outcome Outcome to compare to
         * @return true if given outcome dominates entire box; false otherwise
         */
        bool isDominated(const OutcomeType& outcome) const;

    private:
        /**
         * Get interval intersection with respect to given dimension and given box
         * @param index Interval index to take intersection
         * @param other Box to consider intersection
         * @return Interval intersection
         */
        Interval getIntervalIntersection(std::size_t index,
                                         const RectangularBox& other) const;

        /**
         * Constructor: constructs box first_beg x ... x (first_end-1) x second x third_bex x ... x (third_end-1)
         * @param first_beg Iterator referring to interval
         * @param first_end Iterator referring to past-the-end interval
         * @param second Middle interval
         * @param third_beg Iterator referring to interval
         * @param third_end Iterator referring to past-the-end interval
         */
        RectangularBox(std::vector<Interval>::const_iterator first_beg,
                       std::vector<Interval>::const_iterator first_end,
                       Interval second,
                       std::vector<Interval>::const_iterator third_beg,
                       std::vector<Interval>::const_iterator third_end);

        std::vector<Interval> box_; ///< Container storing interval of rectangular box
    };

    /**
     * @class Polyscip
     * @brief Class for PolySCIP solver functions
     */
    class Polyscip {
    public:

        /**
         * Different statuses of PolySCIP solver
         */
        enum class PolyscipStatus {
            Unsolved, ///< Initial status after calling public constructor
            ProblemRead, ///< Status after problem instance was read successfully
            LexOptPhase, ///< Status after lexicographic optimal results were computed
            WeightSpacePhase, ///< Status while results of weight space polyhedron are computed
            TwoProjPhase, ///< Status while computing 2-dimensional non-dominated projection results
            Finished, ///< Status if problem was solved successfully
            TimeLimitReached, ///< Status if given time limit was reached
            Error ///< Status if an error occured
        };

        using ObjPair = std::pair<std::size_t, std::size_t>; ///< Pair of objectives indices

        /**
         * Default constructor
         * @param argc Argument count
         * @param argv Argument vector
         */
        explicit Polyscip(int argc,
                          const char *const *argv);

        /**
         * Destructor
         */
        ~Polyscip();

        /**
         * Read multi-objective problem file
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE readProblem();

        /**
         * Compute non-dominated points of given problem
         * @attention readProblem() needs to be called before
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE computeNondomPoints();

        /**
         * Indicates whether results shold be written to a file
         * @return true if results should be written to a file; otherwise false
         */
        bool writeResults() const {return cmd_line_args_.writeResults();};

        /**
         * Write results to file named 'solutions_name-of-problem-file.txt'
         */
        void writeResultsToFile() const;

        /**
         * Print results
         * @param os Output stream to print to
         */
        void printResults(std::ostream &os = std::cout) const;

        /**
         * Print PolySCIP status
         * @param os Output stream to print to
         */
        void printStatus(std::ostream& os = std::cout) const;

        /**
         * Get PolySCIP status
         * @return Current PolySCIP status
         */
        PolyscipStatus getStatus() const;

        /**
         * Get number of bounded results
         * @return Number of computed bounded results
         */
        std::size_t numberOfBoundedResults() const;

        /**
         * Get number of unbounded results
         * @return Number of computed unbounded results
         */
        std::size_t numberofUnboundedResults() const;

        /**
         * Indicates whether dominated results were computed
         * @return true if dominated bounded results were computed
         */
        bool dominatedPointsFound() const;

        /**
         * Get iterator to beginning of bounded results
         * @return Const_iterator to beginning of bounded results
         */
        ResultContainer::const_iterator boundedCBegin() {return bounded_.cbegin();};

        /**
         * Get iterator to past-the-end of bounded results
         * @return Const_iterator to past-the-end of bounded results
         */
        ResultContainer::const_iterator boundedCEnd() {return bounded_.cend();};

    private:

        /**
         * Check whether file can be opened
         * @param filename Name of file to open
         * @return true if corresponding file can be opened; false otherwise
         */
        bool filenameIsOkay(const std::string &filename);

        /**
         * Compute lexicographic optimal results
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE computeLexicographicOptResults(std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                                    std::vector<std::vector<ValueType>>& orig_vals);

        /**
         * Compute lexicographic optimal result with given objective having highest preference
         * @param obj Objective with highest preference
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE computeLexicographicOptResult(std::size_t obj,
                                              std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                              std::vector<std::vector<ValueType>>& orig_vals);


        /**
         * Checks whether results corresponding to given iterator is dominated or equal to other given elements
         * @param it Const_iterator corresponding to result
         * @param beg_it Const_iterator to beginning of result container
         * @param end_it Const_iterator to past-the-end of result container
         * @return true if result given by it is dominated or equal to other given results; false otherwise
         */
        bool isDominatedOrEqual(ResultContainer::const_iterator it,
                                ResultContainer::const_iterator beg_it,
                                ResultContainer::const_iterator end_it) const;

        /**
         * Set weighted objective: weight * (c_1,...,c_k) \\cdot x
         * @param weight Weight
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE setWeightedObjective(const WeightType& weight);

        /**
         * Solves currently considered SCIP instance
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE solve();

        /**
         * Resolve INFORUNBD SCIP status to either infeasible or unbounded
         * @param weight Weight yielding INFORUNBD status
         * @param with_presolving Indicates whether presolving should be used or not
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_STATUS separateINFORUNBD(const WeightType& weight,
                                      bool with_presolving = true);

        /**
         * Handle SCIP status that is neither optimal nor unbounded
         * @param status Current SCIP status
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE handleNonOptNonUnbdStatus(SCIP_STATUS status);

        /**
         * Handle unbounded SCIP status
         * @param check_if_new_result Indicates whether to check if computed results is already known
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE handleUnboundedStatus(bool check_if_new_result=false);

        /**
         * Indicates whether given outcomes coincide within some epsilon error
         * @param a First outcome to compare
         * @param b Second outcome to compare
         * @param epsilon Allowed error
         * @return true if outcomes coincides; false otherwise
         */
        static bool outcomesCoincide(const OutcomeType& a,
                                     const OutcomeType& b,
                                     double epsilon);

        /**
         * Indicates whether given outcome was not computed before
         * @param outcome Outcome to check
         * @param outcome_is_bounded Indicates whether given outcome is bounded or unbounded
         * @return true if given outcome was not computed before; false otherwise
         */
        bool outcomeIsNew(const OutcomeType& outcome,
                          bool outcome_is_bounded) const;

        /**
         * Indicates whether given outcome is new with respect to other given results
         * @param outcome Outcome to check
         * @param beg Const_iterator to beginning of result container
         * @param last Const_iterator to past-the-end of result container
         * @return true if given outcome does not coincide with outcomes; false otherwise
         */
        bool outcomeIsNew(const OutcomeType& outcome,
                          ResultContainer::const_iterator beg,
                          ResultContainer::const_iterator last) const;

        /**
         * Get computed result
         * @param outcome_is_bounded Indicates whether previous computation yielded unbounded status
         * @param primal_sol Corresponding SCIP primal solution pointer if previous computation yielded optimal status
         * @return Result type
         */
        Result getResult(bool outcome_is_bounded = false,
                         SCIP_SOL* primal_sol = nullptr);

        /**
         * Get bounded optimal result
         * @return Result type
         */
        Result getOptimalResult();

        /**
         * Print objective
         * @param obj_no Corresponding index of objective
         * @param nonzero_indices Indices of variables with non-zero coefficients
         * @param nonzero_vals Corresponding non-zero coefficient variable values
         * @param os Output stream to write to
         */
        void printObjective(std::size_t obj_no,
                            const std::vector<int>& nonzero_indices,
                            const std::vector<SCIP_Real>& nonzero_vals,
                            std::ostream& os = std::cout) const;

        /**
         * Indicates whether objective given by index is redundant
         * @param begin_nonzeros begin_nonzeros[i+1] = begin_nonzeros[i] + obj_probdata->getNumberNonzeroCoeffs(i)
         * @param obj_to_nonzero_indices indices of non-zero variables for each objective
         * @param obj_to_nonzero_values non-zero variables for each objective
         * @param index index of objective to check
         * @return true if checked objective is redundant; false otherwise
         */
        bool objIsRedundant(const std::vector<int>& begin_nonzeros,
                            const std::vector<std::vector<int>>& obj_to_nonzero_indices,
                            const std::vector<std::vector<SCIP_Real>>& obj_to_nonzero_values,
                            std::size_t index) const;

        /**
         * Compute non-dominated extreme point results
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE computeWeightSpaceResults();

        /**
         * Compute bounded non-dominated extreme points for objective for which unbounded ray exits
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE computeBoundedNondomResultsForUnbdObjs();

        /**
         * Compute non-dominated points which are not lexicographically optimal
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE computeNonLexicographicNondomResults(const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                                          const std::vector<std::vector<ValueType>>& orig_vals);

        /**
         * Compute non-dominated points via subproblems with weighted Tchebycheff norm
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @param obj_1 Index of first considered objective
         * @param obj_2 Index of second considered objective
         * @return Container of non-dominated outcomes which are also non-dominated for projection onto obj_1 and obj_2
         */
        std::vector<OutcomeType> solveWeightedTchebycheff(const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                                          const std::vector<std::vector<ValueType>>& orig_vals,
                                                          std::size_t obj_1,
                                                          std::size_t obj_2);

        /**
         * Compute disjoint rectangular boxes from given feasible rectangular boxes
         * @param feasible_boxes List of feasible boxes
         * @return Vector of disjoint feasible rectangular boxes
         */
        std::vector<RectangularBox> computeDisjointBoxes(std::list<RectangularBox>&& feasible_boxes) const;

        /**
         * Compute feasible rectangular boxes
         * @param proj_nondom_outcomes Non-dominated outcomes which are non-dominated for objective pair
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @return List of feasible rectangular boxes
         */
        std::list<RectangularBox> computeFeasibleBoxes(
                const std::map<ObjPair, std::vector<OutcomeType>> &proj_nondom_outcomes,
                const std::vector<std::vector<SCIP_VAR *>> &orig_vars,
                const std::vector<std::vector<ValueType>> &orig_vals);

        /**
         * Compute locally non-dominated results in given rectangular box
         * @param box Rectangular box
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @return Container with locally non-dominated results
         */
        ResultContainer computeNondomPointsInBox(const RectangularBox& box,
                                                 const std::vector<std::vector<SCIP_VAR *>>& orig_vars,
                                                 const std::vector<std::vector<ValueType>>& orig_vals);

        /**
         * Indicates whether given outcome is globally dominated
         * @param outcome Outcome to check for dominance
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @return true if given outcome is dominated; false otherwise
         */
        bool boxResultIsDominated(const OutcomeType& outcome,
                                  const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                  const std::vector<std::vector<ValueType>>& orig_vals);


        /**
         * Create constraint: new_var - beta_i * orig_vals \\cdot orig_vars >= - beta_i * rhs
         * @param new_var Non-original variable
         * @param orig_vars Container storing original problem variables with non-zero coefficients for each objective
         * @param orig_vals Container storing original non-zero objective coefficients for each objective
         * @param rhs rhs value
         * @param beta_i coefficient
         * @return Pointer to corresponding SCIP constraint
         */
        SCIP_CONS* createNewVarTransformCons(SCIP_VAR *new_var,
                                             const std::vector<SCIP_VAR *> &orig_vars,
                                             const std::vector<ValueType> &orig_vals,
                                             const ValueType &rhs,
                                             const ValueType &beta_i);

        /**
         * Create constraint: lhs <= vals \\cdot vars <= rhs
         * @param vars Considered variables
         * @param vals Considered coefficient values
         * @param lhs lhs value
         * @param rhs rhs value
         * @return Pointer to corresponding SCIP constraint
         */
        SCIP_CONS* createObjValCons(const std::vector<SCIP_VAR *>& vars,
                                    const std::vector<ValueType>& vals,
                                    const ValueType& lhs,
                                    const ValueType& rhs);

        /**
         * Computes non-dominated point which fulfills: obj_val_cons1 = obj_val_cons1_rhs and obj_val_cons2 = obj_val_cons2_rhs
         * @param obj_val_cons1 First constraint to consider
         * @param obj_val_cons2 Second constraint to consider
         * @param obj_val_cons1_rhs Corresponding rhs of first constraint
         * @param obj_val_cons2_rhs Corresponding rhs of second constraint
         * @param obj_1 Considered objective index corresponding to first constraint
         * @param obj_2 Considered objective index corresponding to second constraint
         * @param results Container to store computed non-dominated result
         * @return SCIP_OKAY if everything worked; otherwise a suitable error code is passed
         */
        SCIP_RETCODE computeNondomProjResult(SCIP_CONS* obj_val_cons1,
                                             SCIP_CONS* obj_val_cons2,
                                             ValueType obj_val_cons1_rhs,
                                             ValueType obj_val_cons2_rhs,
                                             std::size_t obj_1,
                                             std::size_t obj_2,
                                             ResultContainer &results);

        /**
         * Indicates whether unbounded results were computed
         * @return true if unbounded results were computed; false otherwise
         */
        bool unboundedResultsExist() const {return !unbounded_.empty();};

        /**
         * Print solution
         * @param sol Solution to print
         * @param os Output stream to write to
         */
        void printSol(const SolType& sol,
                      std::ostream& os) const;

        /**
         * Print outcome
         * @param outcome Outcome to print
         * @param os Output stream to write to
         * @param desc Description to print before given outcome
         */
        void outputOutcome(const OutcomeType &outcome,
                           std::ostream& os,
                           const std::string desc ="") const;

        /**
         * Constructor
         * @param cmd_line_args Command line parameter object
         * @param scip SCIP pointer
         * @param no_objs Number of considered objective
         * @param clock_total Clock measuring total computation time
         */
        explicit Polyscip(const CmdLineArgs& cmd_line_args,
                          SCIP* scip,
                          std::size_t no_objs,
                          SCIP_CLOCK *clock_total);

        CmdLineArgs cmd_line_args_; ///< Object containing command line parameter information
        PolyscipStatus polyscip_status_; ///< Current PolySCIP status
        SCIP* scip_; ///< SCIP pointer
        SCIP_Objsense obj_sense_; ///< Objective sense of given problem
        std::size_t no_objs_; ///< Considered number of objectives
        SCIP_CLOCK* clock_total_; ///< Clock measuring the time needed for the entire computation
        bool only_weight_space_phase_; ///< Indicates whether only non-dominated extreme points should be computed
        bool is_sub_prob_; ///< Indicates whether PolySCIP instance belongs to subproblem
        std::unique_ptr<WeightSpacePolyhedron> weight_space_poly_; ///< Pointer holding weight space polyhedron object
        ResultContainer bounded_; ///< Container storing bounded non-dominated results
        ResultContainer unbounded_; ///< Container storing unbounded non-dominated results
        std::vector<std::size_t> unbd_orig_objs_; ///< Container storing objectives indices for which unbounded rays were found
    };

}

#endif //POLYSCIP_SRC_POLYSCIP_H_INCLUDED
