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

#include "polyscip.h"

#include <algorithm> //std::transform, std::max
#include <array>
#include <cmath> //std::abs
#include <cstddef> //std::size_t
#include <fstream>
#include <functional> //std::plus
#include <iomanip> //std::set_precision
#include <iostream>
#include <iterator> //std::advance
#include <limits>
#include <list>
#include <ostream>
#include <memory> //std::addressof
#include <numeric> //std::inner_product
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "polytope_representation.h"
#include "scip/scip.h"
#include "objscip/objscipdefplugins.h"
#include "cmd_line_args.h"
#include "global_functions.h"
#include "polyscip_types.h"
#include "prob_data_objectives.h"
#include "ReaderMOP.h"
#include "weight_space_polyhedron.h"

using std::addressof;
using std::array;
using std::begin;
using std::cout;
using std::end;
using std::ostream;
using std::size_t;
using std::string;
using std::vector;

namespace polyscip {

    using DDMethod = polytoperepresentation::DoubleDescriptionMethod;

    Polyscip::Polyscip(int argc, const char *const *argv)
            : cmd_line_args_(argc, argv),
              polyscip_status_(PolyscipStatus::Unsolved),
              scip_(nullptr),
              obj_sense_(SCIP_OBJSENSE_MINIMIZE), // default objective sense is minimization
              clock_total_(nullptr)
    {
        if (cmd_line_args_.hasTimeLimit() && cmd_line_args_.getTimeLimit() <= 0)
            throw std::domain_error("Invalid time limit.");
        if (cmd_line_args_.hasParameterFile() && !filenameIsOkay(cmd_line_args_.getParameterFile()))
            throw std::invalid_argument("Invalid parameter settings file.");
        if (!filenameIsOkay(cmd_line_args_.getProblemFile()))
            throw std::invalid_argument("Invalid problem file.");

        SCIPcreate(&scip_);
        assert (scip_ != nullptr);
        SCIPincludeDefaultPlugins(scip_);
        SCIPincludeObjReader(scip_, new ReaderMOP(scip_), TRUE);
        SCIPcreateClock(scip_, addressof(clock_total_));
        if (cmd_line_args_.hasParameterFile())
            SCIPreadParams(scip_, cmd_line_args_.getParameterFile().c_str());
    }

    Polyscip::~Polyscip() {
        SCIPfreeClock(scip_, addressof(clock_total_));
        SCIPfree(addressof(scip_));
    }


    SCIP_RETCODE Polyscip::computeNondomPoints() {
        SCIP_CALL( SCIPstartClock(scip_, clock_total_) );
        SCIP_CALL( computeSupported() );

        if (cmd_line_args_.withUnsupported() && polyscip_status_ == PolyscipStatus::CompUnsupportedPhase) {
            SCIP_CALL( computeUnsupported() );
        }
        deleteWeaklyNondomResults();
        SCIP_CALL( SCIPstopClock(scip_, clock_total_) );
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::computeUnitWeightOutcomes() {
        polyscip_status_ = PolyscipStatus::InitPhase;
        auto cur_opt_vals = OutcomeType(considered_objs_.size(), std::numeric_limits<ValueType>::max());
        auto weight = WeightType(considered_objs_.size(),0.);
        for (size_t unit_weight_index=0; unit_weight_index!=considered_objs_.size(); ++unit_weight_index) {
            if (polyscip_status_ != PolyscipStatus::InitPhase)
                break;
            auto supported_size_before = supported_.size();
            weight[unit_weight_index] = 1.;
            SCIP_CALL(setWeightedObjective(weight));
            SCIP_CALL(solve());
            auto scip_status = SCIPgetStatus(scip_);
            if (scip_status == SCIP_STATUS_INFORUNBD)
                scip_status = separateINFORUNBD(weight);

            if (scip_status == SCIP_STATUS_OPTIMAL) {
                SCIP_CALL( handleOptimalStatus(weight, cur_opt_vals[unit_weight_index]) );
            }
            else if (scip_status == SCIP_STATUS_UNBOUNDED) {
                SCIP_CALL( handleUnboundedStatus(true) );
            }
            else {
                SCIP_CALL( handleNonOptNonUnbdStatus(scip_status) );
            }

            if (supported_size_before < supported_.size()) {
                std::transform(begin(cur_opt_vals),
                               end(cur_opt_vals),
                               begin(supported_.back().second),
                               begin(cur_opt_vals),
                               [](ValueType val1, ValueType val2){return std::min<ValueType>(val1, val2);});
            }
            weight[unit_weight_index] = 0.;
        }
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::computeUnsupported() {
        if (SCIPisTransformed(scip_))
            SCIP_CALL( SCIPfreeTransform(scip_) );
        // change objective values of existing variabless to zero
        auto vars = SCIPgetOrigVars(scip_);
        auto no_vars = SCIPgetNOrigVars(scip_);
        for (auto i=0; i<no_vars; ++i) {
            SCIP_CALL( SCIPchgVarObj(scip_, vars[i], 0.) );
        }
        // add new variable with objective value = 1 (for transformed Tchebycheff norm objective)
        SCIP_VAR* z = nullptr;
        SCIP_CALL( SCIPcreateVarBasic(scip_,
                                      addressof(z),
                                      "z",
                                      -SCIPinfinity(scip_),
                                      SCIPinfinity(scip_),
                                      1,
                                      SCIP_VARTYPE_CONTINUOUS) );
        assert (z != nullptr);
        SCIP_CALL( SCIPaddVar(scip_, z) );

        // get variables (excluding new variable z) with nonzero objective coefficients
        auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        auto nonzero_obj_orig_vars = vector<vector<SCIP_VAR*>>{}; // excluding new variable z
        auto nonzero_obj_orig_vals = vector<vector<ValueType>>{}; // excluding objective value of new variable z

        for (auto obj_ind : considered_objs_) {
            nonzero_obj_orig_vars.push_back(obj_probdata->getNonZeroCoeffVars(obj_ind));
            assert (!nonzero_obj_orig_vars.empty());
            auto nonzero_obj_vals = vector<ValueType>{};
            std::transform(nonzero_obj_orig_vars.back().cbegin(),
                           nonzero_obj_orig_vars.back().cend(),
                           std::back_inserter(nonzero_obj_vals),
                           [obj_ind, obj_probdata](SCIP_VAR *var) { return obj_probdata->getObjCoeff(var, obj_ind); });
            nonzero_obj_orig_vals.push_back(std::move(nonzero_obj_vals));
        }

        // sort supported bounded results lexicographically and compute unsupported points for successive pairs of supported results
        std::sort(begin(supported_), end(supported_), [](const Result& res1, const Result& res2){return res1.second < res2.second;});
        auto res = supported_.cbegin();
        while (res != std::prev(supported_.cend()) && polyscip_status_==PolyscipStatus::CompUnsupportedPhase) {
            auto& pred = res->second;
            auto& succ = std::next(res)->second;
            SCIP_CALL( computeUnsupported(z, nonzero_obj_orig_vars, nonzero_obj_orig_vals, pred, succ) );
            ++res;
        }
        if (polyscip_status_ == PolyscipStatus::CompUnsupportedPhase)
            polyscip_status_ = PolyscipStatus::Finished;

        // clean up
        if (SCIPisTransformed(scip_))
            SCIP_CALL( SCIPfreeTransform(scip_) );
        SCIP_Bool var_deleted {FALSE};
        SCIP_CALL( SCIPdelVar(scip_, z, addressof(var_deleted)) );
        assert (var_deleted);
        SCIP_CALL( SCIPreleaseVar(scip_, addressof(z)) );

        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::computeUnsupported(SCIP_VAR* new_var,
                                              vector<vector<SCIP_VAR*>>& orig_vars,
                                              vector<vector<ValueType>>& orig_vals,
                                              const OutcomeType& pred,
                                              const OutcomeType& succ) {
        assert (pred.size() == succ.size());
        auto ref_point = OutcomeType(pred.size(),0.);
        auto upper_point = OutcomeType(pred.size(), 0.);
        std::transform(pred.cbegin(), pred.cend(),
                       succ.cbegin(), begin(ref_point),
                       [](ValueType a, ValueType  b){return std::min(a,b);});
        std::transform(pred.cbegin(), pred.cend(),
                       succ.cbegin(), begin(upper_point),
                       [](ValueType a, ValueType b){return std::max(a,b);});

        // create and add constraints for each objective
        auto constraints = vector<SCIP_CONS*>{};
        assert (orig_vars.size() == orig_vals.size());
        for (size_t i=0; i<orig_vars.size(); ++i) {
            SCIP_CONS* cons = nullptr;
            auto cons_name = "obj_cons_" + std::to_string(i);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip_,
                                                 addressof(cons),
                                                 cons_name.data(),
                                                 global::narrow_cast<int>(orig_vars[i].size()),
                                                 orig_vars[i].data(),
                                                 orig_vals[i].data(),
                                                 ref_point[i],
                                                 upper_point[i]) );
            assert (cons != nullptr);
            SCIP_CALL( SCIPaddCons(scip_, cons) );
            constraints.push_back(cons);
        }
        std::cout << "entering Tchebycheff...";
        solveWeightedTchebycheff(new_var, orig_vars, orig_vals, pred, succ, ref_point);
        std::cout << "...finished.\n";
        // release and delete constraints
        if (SCIPisTransformed(scip_))
            SCIP_CALL( SCIPfreeTransform(scip_) );
        for (auto cons : constraints) {
            SCIP_CALL( SCIPdelCons(scip_, cons) );
            SCIP_CALL( SCIPreleaseCons(scip_, addressof(cons)) );
        }

        return SCIP_OKAY;
    }

    SCIP_CONS* Polyscip::createObjFctCons(SCIP_VAR* new_var,
                                          const vector<SCIP_VAR*>& orig_vars,
                                          const vector<ValueType>& orig_vals,
                                          const ValueType& rhs,
                                          const ValueType& beta_i) {
        auto vars = vector<SCIP_VAR*>(begin(orig_vars), end(orig_vars));
        auto vals = vector<ValueType>(orig_vals.size(), 0.);
        std::transform(begin(orig_vals),
                       end(orig_vals),
                       begin(vals),
                       [beta_i](ValueType val){return -beta_i*val;});
        vars.push_back(new_var);
        vals.push_back(1.);

        SCIP_CONS* cons = nullptr;
        // add contraint new_var  - beta_i* vals \cdot vars >= - beta_i * ref_point[i]
        SCIPcreateConsBasicLinear(scip_,
                                  addressof(cons),
                                  "max_abs_value_transf_cons",
                                  global::narrow_cast<int>(vars.size()),
                                  vars.data(),
                                  vals.data(),
                                  -beta_i*rhs,
                                  SCIPinfinity(scip_));
        assert (cons != nullptr);
        return cons;
    }

    SCIP_RETCODE Polyscip::solveWeightedTchebycheff(SCIP_VAR* new_var,
                                                    const vector<vector<SCIP_VAR*>>& orig_vars,
                                                    const vector<vector<ValueType>>& orig_vals,
                                                    const OutcomeType& pred,
                                                    const OutcomeType& succ,
                                                    const OutcomeType& ref_point) {
        assert (pred.size() == succ.size());
        assert (ref_point.size() == pred.size());
        assert (ref_point.size() == orig_vars.size());
        assert (orig_vars.size() == orig_vals.size());
        assert (pred.size() > 1);
        auto size = ref_point.size();

        auto found_unsupported = ResultContainer{};

        auto beta_space = std::list<betaT>{};
        // add initial element to beta_space
        auto beta = betaT{};
        beta.push_back({1,1});
        for (size_t i=1; i<size; ++i) {
            beta.push_back({cmd_line_args_.getEpsilon(), std::numeric_limits<ValueType>::max()});
        }
        beta_space.push_back(std::move(beta));

        // while beta_space is not empty
        while (!beta_space.empty() && polyscip_status_==PolyscipStatus::CompUnsupportedPhase) {
            auto beta = beta_space.back();
            beta_space.pop_back();
            auto constraints = vector<SCIP_CONS*>{};
            for (size_t i=0; i<size; ++i) {
                assert (beta[i].second - beta[i].first > cmd_line_args_.getEpsilon());
                auto obj_cons = createObjFctCons(new_var,
                                                 orig_vars[i],
                                                 orig_vals[i],
                                                 ref_point[i],
                                                 beta[i].first);
                SCIP_CALL( SCIPaddCons(scip_, obj_cons) );
                constraints.push_back(obj_cons);
            }

            // solve problem
            std::cout << "solving new tchebycheff...";
            SCIP_CALL( solve() );
            std::cout << "finished.\n";

            // check status and result
            auto scip_status = SCIPgetStatus(scip_);
            if (scip_status == SCIP_STATUS_OPTIMAL) {
                auto opt_val = SCIPgetPrimalbound(scip_);
                assert ( opt_val >= 0.);
                auto new_result = getOptimalResult();
                auto upper_beta_values = computeUpperBetaValues(opt_val, new_result.second, ref_point);

                //todo incorporate upper beta vals etc

                // check whether new result coincides with known results

            }
            else if (scip_status == SCIP_STATUS_TIMELIMIT) {
                polyscip_status_ = PolyscipStatus::TimeLimitReached;
            }
            else {
                string error_msg = "unexpected solution status in solveWeightedTchebycheff: " +
                                   std::to_string(SCIPgetStatus(scip_)) + "\n";
                throw std::runtime_error(error_msg);
            }

            // release and delete constraints
            if (SCIPisTransformed(scip_))
                SCIP_CALL( SCIPfreeTransform(scip_) );
            for (auto cons : constraints) {
                SCIP_CALL( SCIPdelCons(scip_, cons) );
                SCIP_CALL( SCIPreleaseCons(scip_, addressof(cons)) );
            }
        }
        return SCIP_OKAY;
    }

    std::vector<ValueType> Polyscip::computeUpperBetaValues(SCIP_Real opt_value,
                                                            const OutcomeType& result,
                                                            const OutcomeType& ref_point) const {
        assert (result.size() == ref_point.size());
        auto upper_values = vector<ValueType>{};
        upper_values.push_back(1.);
        for (size_t i=1; i<ref_point.size(); ++i) {
            assert (SCIPisPositive(scip_, result[i]-ref_point[i]));
            upper_values.push_back(global::narrow_cast<ValueType>(opt_value / (result[i]-ref_point[i]) + cmd_line_args_.getEpsilon()));
        }
        return upper_values;
    }

    SCIP_RETCODE Polyscip::initWeightSpace() {
        SCIP_CALL( computeUnitWeightOutcomes() ); // computes optimal outcomes for all unit weights
        if (polyscip_status_ == PolyscipStatus::InitPhase) {
            if (supported_.empty()) {
                polyscip_status_ = PolyscipStatus::Finished; // all outcomes for unit weights are unbounded
            }
            else {
                auto v_rep = DDMethod(scip_, no_all_objs_, supported_, unbounded_);
                v_rep.computeVRep_Var1();
                std::cout << "No of initial vertices = " << v_rep.size() << "\n";
                std::cout << "Starting initializing WSP...";
                weight_space_poly_ = global::make_unique<WeightSpacePolyhedron>(scip_,
                                                                                considered_objs_.size(),
                                                                                v_rep.moveVRep(),
                                                                                v_rep.moveHRep());
                std::cout << "...finished.\n";
                assert (weight_space_poly_->hasValidSkeleton(considered_objs_.size()));
                polyscip_status_ = PolyscipStatus::WeightSpacePhase;
            }
        }
        return SCIP_OKAY;
    }



    SCIP_STATUS Polyscip::separateINFORUNBD(const WeightType& weight, bool with_presolving) {
        if (!with_presolving)
            SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, TRUE);
        auto zero_weight = WeightType(considered_objs_.size(), 0.);
        setWeightedObjective(zero_weight);
        solve(); // re-compute with zero objective
        if (!with_presolving)
            SCIPsetPresolving(scip_, SCIP_PARAMSETTING_DEFAULT, TRUE);
        auto status = SCIPgetStatus(scip_);
        setWeightedObjective(weight); // re-set to previous objective
        if (status == SCIP_STATUS_INFORUNBD) {
            if (with_presolving)
                separateINFORUNBD(weight, false);
            else
                throw std::runtime_error("INFORUNBD Status for problem with zero objective and no presolving.\n");
        }
        else if (status == SCIP_STATUS_UNBOUNDED) {
            throw std::runtime_error("UNBOUNDED Status for problem with zero objective.\n");
        }
        else if (status == SCIP_STATUS_OPTIMAL) { // previous problem was unbounded
            return SCIP_STATUS_UNBOUNDED;
        }
        return status;
    }


    SCIP_RETCODE Polyscip::handleNonOptNonUnbdStatus(SCIP_STATUS status) {
        assert (status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_UNBOUNDED);
        if (status == SCIP_STATUS_INFORUNBD) {
            throw std::runtime_error("INFORUNBD Status unexpected at this stage.\n");
        }
        else if (status == SCIP_STATUS_TIMELIMIT) {
            polyscip_status_ = PolyscipStatus::TimeLimitReached;
        }
        else {
            polyscip_status_ = PolyscipStatus::Finished;
        }
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::handleUnboundedStatus(bool check_if_new_result) {
        if (!SCIPhasPrimalRay(scip_)) {
            SCIP_CALL( SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, TRUE) );
            if (SCIPisTransformed(scip_))
                SCIP_CALL( SCIPfreeTransform(scip_) );
            SCIP_CALL( solve() );
            SCIP_CALL( SCIPsetPresolving(scip_, SCIP_PARAMSETTING_DEFAULT, TRUE) );
            if (SCIPgetStatus(scip_) != SCIP_STATUS_UNBOUNDED)
                throw std::runtime_error("Status UNBOUNDED expected.\n");
            if (!SCIPhasPrimalRay(scip_))
                throw std::runtime_error("Existence of primal ray expected.\n");
        }
        auto result = getResult(false);
        if (!check_if_new_result || outcomeIsNew(result.second, false)) {
            unbounded_.push_back(std::move(result));
        }
        else {
            global::print(result.second, "Outcome: [", "]");
            cout << "not added to results.\n";
        }
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::handleOptimalStatus(const WeightType& weight,
                                               ValueType current_opt_val) {
        auto best_sol = SCIPgetBestSol(scip_);
        SCIP_SOL *finite_sol{nullptr};
        SCIP_Bool same_obj_val{FALSE};
        SCIP_CALL(SCIPcreateFiniteSolCopy(scip_, addressof(finite_sol), best_sol, addressof(same_obj_val)));

        if (!same_obj_val) {
            auto diff = std::abs(SCIPgetSolOrigObj(scip_, best_sol) -
                                 SCIPgetSolOrigObj(scip_, finite_sol));
            if (diff > 1.0e-5) {
                std::cerr << "absolute value difference after calling SCIPcreateFiniteSolCopy: " << diff << "\n";
                SCIP_CALL(SCIPfreeSol(scip_, addressof(finite_sol)));
                throw std::runtime_error("SCIPcreateFiniteSolCopy: unacceptable difference in objective values.");
            }
        }
        assert (finite_sol != nullptr);
        auto result = getResult(true, finite_sol);

        assert (weight.size() == result.second.size());
        auto weighted_outcome = std::inner_product(weight.cbegin(),
                                                   weight.cend(),
                                                   result.second.cbegin(),
                                                   0.);

        if (SCIPisLT(scip_, weighted_outcome, current_opt_val)) {
            supported_.push_back(std::move(result));
        }
        else {
            global::print(result.second, "Outcome: [", "]");
            cout << "not added to results.\n";
        }

        SCIP_CALL(SCIPfreeSol(scip_, addressof(finite_sol)));
        return SCIP_OKAY;
    }

    Result Polyscip::getResult(bool outcome_is_bounded, SCIP_SOL *primal_sol) {
        SolType sol;
        auto outcome = OutcomeType(considered_objs_.size(),0.);
        auto no_vars = SCIPgetNOrigVars(scip_);
        auto vars = SCIPgetOrigVars(scip_);
        auto objs_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        for (auto i=0; i<no_vars; ++i) {
            auto var_sol_val = outcome_is_bounded ? SCIPgetSolVal(scip_, primal_sol, vars[i]) :
                               SCIPgetPrimalRayVal(scip_, vars[i]);

            if (!SCIPisZero(scip_, var_sol_val)) {
                sol.emplace_back(SCIPvarGetName(vars[i]), var_sol_val);
                auto var_obj_vals = OutcomeType(considered_objs_.size(), 0.);
                for (auto index : considered_objs_) {
                    var_obj_vals[index] = objs_probdata->getObjVal(vars[i], index, var_sol_val);
                }
                std::transform(begin(outcome), end(outcome),
                               begin(var_obj_vals),
                               begin(outcome),
                               std::plus<ValueType>());

            }
        }
        return {sol, outcome};
    }



    Result Polyscip::getOptimalResult() {
        auto best_sol = SCIPgetBestSol(scip_);
        SCIP_SOL *finite_sol{nullptr};
        SCIP_Bool same_obj_val{FALSE};
        auto retcode = SCIPcreateFiniteSolCopy(scip_, addressof(finite_sol), best_sol, addressof(same_obj_val));
        if (retcode != SCIP_OKAY)
            throw std::runtime_error("SCIPcreateFiniteSolCopy: return code != SCIP_OKAY.\n");
        if (!same_obj_val) {
            auto diff = std::abs(SCIPgetSolOrigObj(scip_, best_sol) -
                                 SCIPgetSolOrigObj(scip_, finite_sol));
            if (diff > 1.0e-5) {
                std::cerr << "absolute value difference after calling SCIPcreateFiniteSolCopy: " << diff << "\n";
                SCIPfreeSol(scip_, addressof(finite_sol));
                throw std::runtime_error("SCIPcreateFiniteSolCopy: unacceptable difference in objective values.");
            }
        }
        assert (finite_sol != nullptr);
        auto new_result = getResult(true, finite_sol);
        SCIPfreeSol(scip_, addressof(finite_sol));
        return new_result;
    }


    bool Polyscip::outcomeIsNew(const OutcomeType& outcome, bool outcome_is_bounded) const {
        auto beg_it = outcome_is_bounded ? begin(supported_) : begin(unbounded_);
        auto end_it = outcome_is_bounded ? end(supported_) : end(unbounded_);
        return std::find_if(beg_it, end_it, [&outcome](const Result& res){return outcome == res.second;}) == end_it;
    }

    SCIP_RETCODE Polyscip::solve() {
        if (cmd_line_args_.hasTimeLimit()) { // set SCIP timelimit
            auto remaining_time = std::max(cmd_line_args_.getTimeLimit() -
                                           SCIPgetClockTime(scip_, clock_total_), 0.);
            SCIP_CALL(SCIPsetRealParam(scip_, "limits/time", remaining_time));
        }
        SCIP_CALL( SCIPsolve(scip_) );    // actual SCIP solver call
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::setWeightedObjective(const WeightType& weight){
        if (SCIPisTransformed(scip_))
            SCIP_CALL( SCIPfreeTransform(scip_) );
        auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        assert (obj_probdata != nullptr);
        auto vars = SCIPgetOrigVars(scip_);
        auto no_vars = SCIPgetNOrigVars(scip_);
        for (auto i=0; i<no_vars; ++i) {
            auto val = obj_probdata->getWeightedObjVal(vars[i], weight, considered_objs_);
            SCIP_CALL( SCIPchgVarObj(scip_, vars[i], val) );
        }
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::computeSupported() {
        SCIP_CALL( initWeightSpace() );
        if (polyscip_status_ == PolyscipStatus::WeightSpacePhase) {
            std::cout << "Starting weight space phase...\n";
            while (weight_space_poly_->hasUntestedWeight()) {
                auto untested_weight = weight_space_poly_->getUntestedWeight();
                std::cout << "new weight: ";
                global::print(untested_weight);
                SCIP_CALL( setWeightedObjective(untested_weight) );
                SCIP_CALL( solve() );
                auto scip_status = SCIPgetStatus(scip_);
                if (scip_status == SCIP_STATUS_INFORUNBD)
                    scip_status = separateINFORUNBD(untested_weight);
                if (scip_status == SCIP_STATUS_OPTIMAL) {
                    //if (SCIPisLT(scip_, SCIPgetPrimalbound(scip_), weight_space_poly_->getUntestedVertexWOV(untested_weight))) {
                    auto supported_size_before = supported_.size();
                    SCIP_CALL(handleOptimalStatus(untested_weight,
                                                  weight_space_poly_->getUntestedVertexWOV(
                                                          untested_weight))); //might add bounded result to supported_
                    if (supported_size_before < supported_.size()) {

                        std::cout << "incorporating new outcome...";
                        weight_space_poly_->incorporateNewOutcome(scip_,
                                                                  untested_weight,
                                                                  supported_.back().second); // was added by handleOptimalStatus()
                        std::cout << "...finished.\n";
                    }
                    else {
                        std::cout << "incorporating old outcome...";
                        weight_space_poly_->incorporateKnownOutcome(untested_weight);
                        std::cout << "...finished.\n";
                    }
                    //}
                }
                else if (scip_status == SCIP_STATUS_UNBOUNDED) {
                    SCIP_CALL( handleUnboundedStatus() ); //adds unbounded result to unbounded_
                    std::cout << "incorporating unbounded outcome...";
                    weight_space_poly_->incorporateNewOutcome(scip_,
                                                              untested_weight,
                                                              unbounded_.back().second, // was added by handleUnboundedStatus()
                                                              true);
                    std::cout << "...finished.\n";
                }
                else {
                    SCIP_CALL( handleNonOptNonUnbdStatus(scip_status) ); //polyscip_status_ is set to finished or time limit reached
                    return SCIP_OKAY;
                }
            }
            std::cout << "...finished.\n";
            if (SCIPgetNOrigContVars(scip_) == SCIPgetNOrigVars(scip_)) {   //check whether there exists integer variables
                polyscip_status_ = PolyscipStatus::Finished;
            }
            else {
                polyscip_status_ = PolyscipStatus::CompUnsupportedPhase;
            }
        }
        return SCIP_OKAY;
    }

    void Polyscip::printResults(ostream &os, bool withSolutions) const {
        os << "Number of supported bounded results: " << supported_.size() << "\n";
        for (const auto& result : supported_) {
            printPoint(result.second, os);
            if (withSolutions)
                printSol(result.first, os);
            os << "\n";
        }
        os << "Number of supported unbounded results: " << unbounded_.size() << "\n";
        for (const auto& result : unbounded_) {
            printRay(result.second, os);
            if (withSolutions)
                printSol(result.first, os);
            os << "\n";
        }
        os << "Number of unsupported bounded results: " << unsupported_.size() << "\n";
        for (const auto& result : unsupported_) {
            printPoint(result.second, os);
            if (withSolutions)
                printSol(result.first, os);
            os << "\n";
        }
    }

    void Polyscip::printSol(const SolType& sol, ostream& os) const {
        os << " Non-zero solution variables: ";
        for (const auto& elem : sol)
            os << elem.first << "=" << elem.second << " ";
    }

    void Polyscip::printPoint(const OutcomeType& point, ostream& os) const {
        global::print(point, "Point = [", "]", os);
    }

    void Polyscip::printRay(const OutcomeType& ray, ostream& os) const {
        global::print(ray, "Ray = [", "]", os);
    }

    bool Polyscip::filenameIsOkay(const string& filename) {
        std::ifstream file(filename.c_str());
        return file.good();
    }

    /*void Polyscip::computeNonRedundantObjectives(bool printObjectives) {
        auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        assert (obj_probdata != nullptr);
        auto vars = SCIPgetOrigVars(scip_);
        auto no_objs = obj_probdata->getNoAllObjs();
        auto begin_nonzeros = vector<int>(no_objs,0);
        for (size_t i = 0; i<no_objs-1; ++i)
            begin_nonzeros[i+1] = global::narrow_cast<int>(begin_nonzeros[i] + obj_probdata->getNumberNonzeroCoeffs(i));

        auto obj_to_nonzero_inds = vector< vector<int> >{};
        auto obj_to_nonzero_vals = vector< vector<SCIP_Real> >{};
        for (size_t obj_index=0; obj_index<no_objs; ++obj_index) {
            auto nonzero_vars = obj_probdata->getNonZeroCoeffVars(obj_index);
            auto size = nonzero_vars.size();
            assert (size > 0);
            auto nonzero_inds = vector<int>(size, 0);
            std::transform(begin(nonzero_vars),
                           end(nonzero_vars),
                           begin(nonzero_inds),
                           [](SCIP_VAR* var){return SCIPvarGetProbindex(var);});
            std::sort(begin(nonzero_inds), end(nonzero_inds));

            auto nonzero_vals = vector<SCIP_Real>(size, 0.);
            std::transform(begin(nonzero_inds),
                           end(nonzero_inds),
                           begin(nonzero_vals),
                           [&](int var_ind){return obj_probdata->getObjCoeff(vars[var_ind], obj_index);});

            if (printObjectives)
                printObjective(obj_index, nonzero_inds, nonzero_vals);

            obj_to_nonzero_inds.push_back(std::move(nonzero_inds));
            obj_to_nonzero_vals.push_back(std::move(nonzero_vals));
        }

        considered_objs_.push_back(0);
        for (size_t obj_no=1; obj_no<no_objs; ++obj_no) {
            if (!objIsRedundant(begin_nonzeros,
                               obj_to_nonzero_inds,
                               obj_to_nonzero_vals,
                               obj_no))
                considered_objs_.push_back(obj_no);
            else
                std::cout << "objective no: " << obj_no << " is redundant.\n";
        }
    }*/

    void Polyscip::printObjective(size_t obj_no,
                                  const std::vector<int>& nonzero_indices,
                                  const std::vector<SCIP_Real>& nonzero_vals,
                                  ostream& os) const {
        assert (!nonzero_indices.empty());
        auto size = nonzero_indices.size();
        assert (size == nonzero_vals.size());
        auto obj = vector<SCIP_Real>(global::narrow_cast<size_t>(SCIPgetNOrigVars(scip_)), 0);
        for (size_t i=0; i<size; ++i)
            obj[nonzero_indices[i]] = nonzero_vals[i];
        global::print(obj, std::to_string(obj_no) + ". obj: [", "]", os);
        os << "\n";
    }

    bool Polyscip::objIsRedundant(const vector<int>& begin_nonzeros,
                                  const vector< vector<int> >& obj_to_nonzero_indices,
                                  const vector< vector<SCIP_Real> >& obj_to_nonzero_values,
                                  size_t checked_obj) const {
        bool is_redundant = false;
        auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));

        assert (obj_probdata != nullptr);
        assert (checked_obj >= 1 && checked_obj < obj_probdata->getNoAllObjs());

        SCIP_LPI* lpi;
        auto retcode = SCIPlpiCreate(addressof(lpi), nullptr, "check objective redundancy", SCIP_OBJSEN_MINIMIZE);
        if (retcode != SCIP_OKAY)
            throw std::runtime_error("no SCIP_OKAY for SCIPlpiCreate\n.");

        auto no_cols = global::narrow_cast<int>(checked_obj);
        auto obj = vector<SCIP_Real>(checked_obj, 1.);
        auto lb = vector<SCIP_Real>(checked_obj, 0.);
        auto ub = vector<SCIP_Real>(checked_obj, SCIPlpiInfinity(lpi));
        auto no_nonzero = begin_nonzeros.at(checked_obj);

        auto beg = vector<int>(begin(begin_nonzeros), begin(begin_nonzeros)+checked_obj);
        auto ind = vector<int>{};
        ind.reserve(global::narrow_cast<size_t>(no_nonzero));
        auto val = vector<SCIP_Real>{};
        val.reserve(global::narrow_cast<size_t>(no_nonzero));
        for (size_t i=0; i<checked_obj; ++i) {
            ind.insert(end(ind), begin(obj_to_nonzero_indices[i]), end(obj_to_nonzero_indices[i]));
            val.insert(end(val), begin(obj_to_nonzero_values[i]), end(obj_to_nonzero_values[i]));
        }

        auto no_rows = SCIPgetNOrigVars(scip_);
        auto vars = SCIPgetOrigVars(scip_);
        auto lhs = vector<SCIP_Real>(global::narrow_cast<size_t>(no_rows), 0.);
        for (auto i=0; i<no_rows; ++i)
            lhs[i] = obj_probdata->getObjCoeff(vars[i], checked_obj);
        auto rhs = vector<SCIP_Real>(lhs);

        retcode =  SCIPlpiLoadColLP(lpi,
                                    SCIP_OBJSEN_MINIMIZE,
                                    no_cols,
                                    obj.data(),
                                    lb.data(),
                                    ub.data(),
                                    nullptr,
                                    no_rows,
                                    lhs.data(),
                                    rhs.data(),
                                    nullptr,
                                    no_nonzero,
                                    beg.data(),
                                    ind.data(),
                                    val.data());

        if (retcode != SCIP_OKAY)
            throw std::runtime_error("no SCIP_OKAY for SCIPlpiLoadColLP\n");

        //SCIPlpiWriteLP(lpi, "redundancy_check.lp");

        retcode = SCIPlpiSolvePrimal(lpi);
        if (retcode != SCIP_OKAY)
            throw std::runtime_error("no SCIP_OKAY for SCIPlpiSolvePrimal\n");

        if (SCIPlpiIsPrimalFeasible(lpi)) {
            is_redundant = true;
        }
        else {
            assert (SCIPlpiIsPrimalInfeasible(lpi));
        }

        retcode = SCIPlpiFree(addressof(lpi));
        if (retcode != SCIP_OKAY)
            throw std::runtime_error("no SCIP_OKAY for SCIPlpiFree\n");

        return is_redundant;
    }

    SCIP_RETCODE Polyscip::readProblem() {
        auto filename = cmd_line_args_.getProblemFile();
        SCIP_CALL( SCIPreadProb(scip_, filename.c_str(), "mop") );
        auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        assert (obj_probdata != nullptr);
        no_all_objs_ = obj_probdata->getNoAllObjs();

        if (cmd_line_args_.beVerbose() || cmd_line_args_.checkForRedundantObjs()) {
            auto vars = SCIPgetOrigVars(scip_);
            auto begin_nonzeros = vector<int>(no_all_objs_, 0);
            for (size_t i = 0; i<no_all_objs_-1; ++i)
                begin_nonzeros[i+1] = global::narrow_cast<int>(begin_nonzeros[i] + obj_probdata->getNumberNonzeroCoeffs(i));

            auto obj_to_nonzero_inds = vector< vector<int> >{};
            auto obj_to_nonzero_vals = vector< vector<SCIP_Real> >{};
            for (size_t obj_ind=0; obj_ind<no_all_objs_; ++obj_ind) {
                auto nonzero_vars = obj_probdata->getNonZeroCoeffVars(obj_ind);
                auto size = nonzero_vars.size();
                if (size == 0)
                    throw std::runtime_error(std::to_string(obj_ind) + ". objective is zero objective!");
                auto nonzero_inds = vector<int>(size, 0);
                std::transform(begin(nonzero_vars),
                               end(nonzero_vars),
                               begin(nonzero_inds),
                               [](SCIP_VAR* var){return SCIPvarGetProbindex(var);});
                std::sort(begin(nonzero_inds), end(nonzero_inds));

                auto nonzero_vals = vector<SCIP_Real>(size, 0.);
                std::transform(begin(nonzero_inds),
                               end(nonzero_inds),
                               begin(nonzero_vals),
                               [&](int var_ind){return obj_probdata->getObjCoeff(vars[var_ind], obj_ind);});

                if (cmd_line_args_.beVerbose())
                    printObjective(obj_ind, nonzero_inds, nonzero_vals);

                if (cmd_line_args_.checkForRedundantObjs()) {
                    obj_to_nonzero_inds.push_back(std::move(nonzero_inds));
                    obj_to_nonzero_vals.push_back(std::move(nonzero_vals));
                }
                else {
                    considered_objs_.push_back(obj_ind);
                }
            }

            if (cmd_line_args_.checkForRedundantObjs()) {
                considered_objs_.push_back(0);
                for (size_t obj_no=1; obj_no<no_all_objs_; ++obj_no) {
                    if (!objIsRedundant(begin_nonzeros,
                                        obj_to_nonzero_inds,
                                        obj_to_nonzero_vals,
                                        obj_no))
                        considered_objs_.push_back(obj_no);
                    else
                        cout << "Objective no: " << obj_no << " is redundant.\n";
                }
            }
        }
        else {
            for (size_t obj_ind=0; obj_ind<no_all_objs_; ++obj_ind)
                considered_objs_.push_back(obj_ind);
        }

        if (SCIPgetObjsense(scip_) == SCIP_OBJSENSE_MAXIMIZE) {
            obj_sense_ = SCIP_OBJSENSE_MAXIMIZE;
            // internally we treat problem as min problem and negate objective coefficients
            SCIPsetObjsense(scip_, SCIP_OBJSENSE_MINIMIZE);
            obj_probdata->negateAllCoeffs();
        }
        if (cmd_line_args_.beVerbose()) {
            cout << "Objective sense: ";
            if (obj_sense_ == SCIP_OBJSENSE_MAXIMIZE)
                cout << "MAXIMIZE\n";
            else
                cout << "MINIMIZE\n";
            cout << "Number of considered objectives: " << considered_objs_.size() << "\n";
        }
        return SCIP_OKAY;
    }

    void Polyscip::writeFileForVertexEnumeration() const {
        auto prob_file = cmd_line_args_.getProblemFile();
        size_t prefix = prob_file.find_last_of("/"), //separate path/ and filename.mop
                suffix = prob_file.find_last_of("."),      //separate filename and .mop
                start_ind = (prefix == string::npos) ? 0 : prefix + 1,
                end_ind = (suffix != string::npos) ? suffix : string::npos;
        string file_name = prob_file.substr(start_ind, end_ind - start_ind) + ".ine";
        std::ofstream solfs(file_name);
        if (solfs.is_open()) {
            solfs << "WeightSpacePolyhedron\n";
            solfs << "H-representation\n";
            solfs << "begin\n";
            solfs << supported_.size() + unbounded_.size() + no_all_objs_ << " " << considered_objs_.size() + 1 << " rational\n";
            for (const auto& elem : supported_) {
                global::print(elem.second, "0 ", " -1\n", solfs);
            }
            for (const auto& elem : unbounded_) {
                global::print(elem.second, "0 ", " 0", solfs);
            }
            for (size_t i=0; i<no_all_objs_; ++i) {
                auto ineq = vector<unsigned>(no_all_objs_, 0);
                ineq[i] = 1;
                global::print(ineq, "0 ", " 0\n", solfs);
            }
            solfs << "end\n";
            solfs.close();
        }
        else
            cout << "ERROR writing vertex enumeration file\n.";
    }

    void Polyscip::writeSupportedResults() const {
        auto prob_file = cmd_line_args_.getProblemFile();
        size_t prefix = prob_file.find_last_of("/"), //separate path/ and filename.mop
                suffix = prob_file.find_last_of("."),      //separate filename and .mop
                start_ind = (prefix == string::npos) ? 0 : prefix + 1,
                end_ind = (suffix != string::npos) ? suffix : string::npos;
        string file_name = "solutions_" +
                           prob_file.substr(start_ind, end_ind - start_ind) + ".txt";
        auto write_path = cmd_line_args_.getWritePath();
        if (write_path.back() != '/')
            write_path.push_back('/');
        std::ofstream solfs(write_path + file_name);
        if (solfs.is_open()) {
            printResults(solfs);
            solfs.close();
            cout << "#Solution file " << file_name
            << " written to: " << write_path << "\n";
        }
        else
            cout << "ERROR writing solution file\n.";
    }

    bool Polyscip::isDominatedOrEqual(ResultContainer::const_iterator it) const {
        for (auto curr = supported_.cbegin(); curr != supported_.cend(); ++curr) {
            if (it == curr)
                continue;
            else if (std::equal(begin(curr->second),
                                end(curr->second),
                                begin(it->second),
                                std::less_equal<ValueType>()))
                return true;
        }
        return false;
    }


    void Polyscip::deleteWeaklyNondomResults() {
        auto it = begin(supported_);
        while (it != end(supported_)) {
            if (isDominatedOrEqual(it)) {
                cout << "Deleting weakly non-dominated point.\n";
                it = supported_.erase(it);
            }
            else {
                ++it;
            }
        }
    }

}
