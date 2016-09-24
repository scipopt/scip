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
#include <functional> //std::plus, std::negate, std::function
#include <iomanip> //std::set_precision
#include <iostream>
#include <iterator> //std::advance
#include <limits>
#include <list>
#include <ostream>
#include <map>
#include <memory> //std::addressof, std::unique_ptr
#include <numeric> //std::inner_product
#include <stdexcept>
#include <string>
#include <type_traits>
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
using std::list;
using std::ostream;
using std::pair;
using std::size_t;
using std::string;
using std::vector;

namespace polyscip {

    using DDMethod = polytoperepresentation::DoubleDescriptionMethod;

    TwoDProj::TwoDProj(const OutcomeType& outcome, std::size_t first, std::size_t second)
            : proj_(outcome[first], outcome[second])
    {}

    bool TwoDProj::operator<(TwoDProj other) const {
        return ((this->getFirst() < other.getFirst()) ||
                (this->getFirst() == other.getFirst() && this->getSecond() < other.getSecond()));
    }

    bool TwoDProj::dominates(double eps, const TwoDProj& other) const {
        assert (eps >= 0);
        return (this->getFirst()-eps < other.getFirst() && this->getSecond()-eps < other.getSecond());
    }

    std::ostream& operator<<(std::ostream& os, const TwoDProj& p) {
        os << "Proj = [" << p.proj_.first << ", " << p.proj_.second << "]";
        return os;
    }

    NondomProjections::ProjMap::iterator NondomProjections::add(TwoDProj proj, const OutcomeType& outcome) {
        auto ret = nondom_projections_.emplace(std::move(proj), vector<OutcomeType>{outcome});
        if (!ret.second) {
            ret.first->second.push_back(outcome);
        }
        return ret.first;
    }

    NondomProjections::NondomProjections(double eps,
                                         const ResultContainer &supported,
                                         const ResultContainer &unsupported,
                                         std::size_t first,
                                         std::size_t second)
            : epsilon_(eps)
    {
        assert (first != second);
        assert (!supported.empty());
        assert (first < supported.front().second.size() && second < supported.front().second.size());
        for (const auto& res : supported) {
            add(TwoDProj(res.second, first, second), res.second);
        }
        for (const auto& res : unsupported) {
            add(TwoDProj(res.second, first, second), res.second);
        }
        auto it = begin(nondom_projections_);
        while (it!=std::prev(end(nondom_projections_))) {
            auto next = std::next(it);
            if (it->first.dominates(epsilon_, next->first)) {
                nondom_projections_.erase(next);
            }
            else {
                ++it;
            }
        }
        assert (!nondom_projections_.empty());
        current_ = begin(nondom_projections_);
    }

    void NondomProjections::update() {
        assert (current_ != std::prev(end(nondom_projections_)) && current_ != end(nondom_projections_));
        ++current_;
    }

    void NondomProjections::update(TwoDProj proj, const OutcomeType& outcome) {
        assert (current_ != std::prev(end(nondom_projections_)) && current_ != end(nondom_projections_));
        if (proj.dominates(epsilon_, current_->first)) {
            auto end_del = std::next(current_);
            while (proj.dominates(epsilon_, end_del->first)) {
                ++end_del;
            }
            nondom_projections_.erase(current_, end_del);
            current_ = this->add(std::move(proj), outcome);
        }
        else if (proj.dominates(epsilon_, std::next(current_)->first)) {
            auto end_del = current_;
            std::advance(end_del, 2);
            while (proj.dominates(epsilon_, end_del->first)) {
                ++end_del;
            }
            nondom_projections_.erase(std::next(current_), end_del);
            this->add(std::move(proj), outcome);
        }
        else {
            auto ret = this->add(std::move(proj), outcome);
            assert (current_ == std::prev(ret)); // assert that new element was inserted right after current_
        }
    }

    bool NondomProjections::finished() const {
        assert (current_ != end(nondom_projections_));
        return current_ == std::prev(end(nondom_projections_));
    }

    Polyscip::Polyscip(int argc, const char *const *argv)
            : cmd_line_args_(argc, argv),
              polyscip_status_(PolyscipStatus::Unsolved),
              scip_(nullptr),
              obj_sense_(SCIP_OBJSENSE_MINIMIZE), // default objective sense is minimization
              no_objs_(0),
              clock_total_(nullptr),
              is_subproblem_(false)
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

    Polyscip::Polyscip(const CmdLineArgs& cmd_line_args,
                       SCIP *scip,
                       SCIP_Objsense obj_sense,
                       pair<size_t, size_t> objs_to_be_ignored,
                       SCIP_CLOCK *clock_total)
            : cmd_line_args_{cmd_line_args},
              polyscip_status_{PolyscipStatus::Unsolved},
              scip_{scip},
              obj_sense_{obj_sense},
              clock_total_{clock_total},
              is_subproblem_(true)
    {
        auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        obj_probdata->ignoreObjectives(objs_to_be_ignored.first, objs_to_be_ignored.second);
        no_objs_ = obj_probdata->getNoObjs();
    }

    Polyscip::~Polyscip() {
        if (is_subproblem_) {
            auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
            obj_probdata->unignoreObjectives();
        }
        else {
            SCIPfreeClock(scip_, addressof(clock_total_));
            SCIPfree(addressof(scip_));
        }
    }


    void Polyscip::printStatus(std::ostream& os) const {
        os << "Number of extremal supported bounded results: " << supported_.size() << "\n";
        os << "Number of supported unbounded results: " << unbounded_.size() << "\n";
        os << "Number of non-extremal bounded results: " << unsupported_.size() << "\n";
        switch(polyscip_status_) {
            case PolyscipStatus::CompUnsupportedPhase:
                os << "PolySCIP Status: ComputeUnsupportedPhase\n";
                break;
            case PolyscipStatus::Error:
                os << "PolySCIP Status: Error\n";
                break;
            case PolyscipStatus::Finished:
                os << "PolySCIP Status: Successfully finished\n";
                break;
            case PolyscipStatus::InitPhase:
                os << "PolySCIP Status: InitPhase\n";
                break;
            case PolyscipStatus::TimeLimitReached:
                os << "PolySCIP Status: TimeLimitReached\n";
                break;
            case PolyscipStatus::Unsolved:
                os << "PolySCIP Status: Unsolved\n";
                break;
            case PolyscipStatus::WeightSpacePhase:
                os << "PolySCIP Status: WeightSpacePhase\n";
                break;
        }
    }

    SCIP_RETCODE Polyscip::computeNondomPoints() {
        SCIP_CALL( SCIPstartClock(scip_, clock_total_) );
        SCIP_CALL( computeSupported() );
        deleteWeaklyNondomSupportedResults();
        if (polyscip_status_ == PolyscipStatus::CompUnsupportedPhase) {
            SCIP_CALL( computeUnsupported() );
        }
        SCIP_CALL( SCIPstopClock(scip_, clock_total_) );
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::computeUnitWeightOutcomes() {
        polyscip_status_ = PolyscipStatus::InitPhase;
        auto cur_opt_vals = OutcomeType(no_objs_, std::numeric_limits<ValueType>::max());
        auto weight = WeightType(no_objs_, 0.);
        for (size_t unit_weight_index=0; unit_weight_index!=no_objs_; ++unit_weight_index) {
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

    /*bool Polyscip::ValPairCmp(Polyscip::ValPair p1, Polyscip::ValPair p2) {
        return (SCIPisLT(scip_, p1.first, p2.first) ||
                (SCIPisEQ(scip_, p1.first, p2.first) && SCIPisLT(scip_, p1.second, p2.second)));
    }*/

    /*Polyscip::ValPairMap Polyscip::getProjectedNondomPoints(std::size_t obj_1, std::size_t obj_2) const {
        auto projected_points = ValPairMap{};
        for (const auto& sup : supported_) {
            auto proj = ValPair(sup.second[obj_1], sup.second[obj_2]);
            if (projected_points.count(proj) == 0)
                projected_points.emplace(proj, vector<OutcomeType>{});
            projected_points[proj].push_back(sup.second);
        }
        for (const auto& unsup : unsupported_) {
            auto proj = ValPair(unsup.second[obj_1], unsup.second[obj_2]);
            if (projected_points.count(proj) == 0)
                projected_points.emplace(proj, vector<OutcomeType>{});
            projected_points[proj].push_back(unsup.second);
        }

        for (auto it=std::next(projected_points.begin()); it!=projected_points.end(); ++it) {
            auto prev_it = std::prev(it);
            if (SCIPisGE(scip_, it->first.second, prev_it->first.second)) {
                assert (SCIPisLE(scip_, prev_it->first.first, it->first.first));
                projected_points.erase(it);
                it = prev_it;
            }
        }
        return projected_points;
    }*/


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
        auto nonzero_obj_orig_vars = vector<vector<SCIP_VAR*>>{};
        auto nonzero_obj_orig_vals = vector<vector<ValueType>>{};

        for (size_t obj_ind=0; obj_ind < no_objs_; ++obj_ind) {
            nonzero_obj_orig_vars.push_back(obj_probdata->getNonZeroCoeffVars(obj_ind)); // excluding new variable z
            assert (!nonzero_obj_orig_vars.empty());
            auto nonzero_obj_vals = vector<ValueType>{};
            std::transform(nonzero_obj_orig_vars.back().cbegin(),
                           nonzero_obj_orig_vars.back().cend(),
                           std::back_inserter(nonzero_obj_vals),
                           [obj_ind, obj_probdata](SCIP_VAR *var) { return obj_probdata->getObjCoeff(var, obj_ind); });
            nonzero_obj_orig_vals.push_back(std::move(nonzero_obj_vals)); // excluding objective value of new variable z
        }

        // consider all (k over 2 ) combinations of considered objective functions
        for (size_t obj_1=0; obj_1!=no_objs_-1; ++obj_1) {
            for (auto obj_2=obj_1+1; obj_2!=no_objs_; ++obj_2) {
                if (polyscip_status_ == PolyscipStatus::CompUnsupportedPhase) {
                    if (cmd_line_args_.beVerbose()) {
                        std::cout << "Considered objective projection: obj_1=" << obj_1 << ", obj_2=" << obj_2 << "\n";
                    }
                    /*solveWeightedTchebycheff(z,
                                             nonzero_obj_orig_vars,
                                             nonzero_obj_orig_vals,
                                             {obj_1, obj_2},
                                             std::move(getProjectedNondomPoints(obj_1, obj_2)));*/
                    solveWeightedTchebycheff(z,
                                             nonzero_obj_orig_vars,
                                             nonzero_obj_orig_vals,
                                             obj_1, obj_2);
                }
            }
        }

        // clean up
        if (SCIPisTransformed(scip_))
            SCIP_CALL( SCIPfreeTransform(scip_) );
        SCIP_Bool var_deleted = FALSE;
        SCIP_CALL( SCIPdelVar(scip_, z, addressof(var_deleted)) );
        assert (var_deleted);
        SCIP_CALL( SCIPreleaseVar(scip_, addressof(z)) );

        if (polyscip_status_ == PolyscipStatus::CompUnsupportedPhase)
            polyscip_status_ = PolyscipStatus::Finished;

        return SCIP_OKAY;
    }

    SCIP_CONS* Polyscip::createNewVarTransformCons(SCIP_VAR *new_var,
                                                   const vector<SCIP_VAR *> &orig_vars,
                                                   const vector<ValueType> &orig_vals,
                                                   const ValueType &rhs,
                                                   const ValueType &beta_i) {
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
                                  "new_variable_transformation_constraint",
                                  global::narrow_cast<int>(vars.size()),
                                  vars.data(),
                                  vals.data(),
                                  -beta_i*rhs,
                                  SCIPinfinity(scip_));
        assert (cons != nullptr);
        return cons;
    }

    /** create constraint:
     *
     * @param orig_vars
     * @param orig_vals
     * @param lhs
     * @param rhs
     * @return
     */
    SCIP_CONS* Polyscip::createObjValCons(const vector<SCIP_VAR *>& vars,
                                          const vector<ValueType>& vals,
                                          const ValueType& lhs,
                                          const ValueType& rhs) {
        SCIP_CONS* cons = nullptr;
        std::remove_const<const vector<SCIP_VAR*>>::type non_const_vars(vars);
        std::remove_const<const vector<ValueType>>::type non_const_vals(vals);
        SCIPcreateConsBasicLinear(scip_,
                                  addressof(cons),
                                  "lhs <= c_i^T x <= rhs",
                                  global::narrow_cast<int>(vars.size()),
                                  non_const_vars.data(),
                                  non_const_vals.data(),
                                  lhs,
                                  rhs);
        assert (cons != nullptr);
        return cons;
    }

    SCIP_RETCODE Polyscip::computeNondomResult(SCIP_VAR *var_z,
                                               SCIP_CONS *cons1,
                                               SCIP_CONS *cons2,
                                               const ValueType &rhs_cons1,
                                               const ValueType &rhs_cons2,
                                               ResultContainer &results) {
        // set new objective value constraints
        SCIP_CALL( SCIPchgRhsLinear(scip_, cons1, rhs_cons1) );
        SCIP_CALL( SCIPchgRhsLinear(scip_, cons2, rhs_cons2) );
        // set new objective function
        SCIP_CALL( setWeightedObjective(WeightType(no_objs_, 1.)) );
        assert( SCIPvarGetObj(var_z) == 0. );

        // solve auxiliary problem
        SCIP_CALL( solve() );
        auto scip_status = SCIPgetStatus(scip_);
        if (scip_status == SCIP_STATUS_OPTIMAL) {
            auto nondom_result = getOptimalResult();
            deleteVarNameFromResult(var_z, nondom_result);
            results.push_back(std::move(nondom_result));
        }
        else if (scip_status == SCIP_STATUS_TIMELIMIT) {
            polyscip_status_ = PolyscipStatus::TimeLimitReached;
        }
        else {
            std::cout << "unexpected SCIP status in computeNewResults: " +
                         std::to_string(SCIPgetStatus(scip_)) + "\n";
            polyscip_status_ = PolyscipStatus::Error;
        }

        // unset objective function
        SCIP_CALL( setWeightedObjective( WeightType(no_objs_, 0.)) );
        SCIP_CALL( SCIPchgVarObj(scip_, var_z, 1.) );
        return SCIP_OKAY;
    }

    /*bool Polyscip::lhsLessEqualrhs(const ValPair &lhs, const ValPair &rhs) const {
        if (lhs.first - cmd_line_args_.getEpsilon() < rhs.first && lhs.second - cmd_line_args_.getEpsilon() < rhs.second)
            return true;
        else
            return false;
    }*/

    SCIP_RETCODE Polyscip::solveWeightedTchebycheff(SCIP_VAR* new_var,
                                                    const std::vector<std::vector<SCIP_VAR*>>& orig_vars,
                                                    const std::vector<std::vector<ValueType>>& orig_vals,
                                                    std::size_t obj_1,
                                                    std::size_t obj_2) {

        assert (orig_vars.size() == orig_vals.size());
        assert (orig_vals.size() == no_objs_);
        assert (obj_1 < no_objs_ && obj_2 < no_objs_);

        auto nondom_projs = NondomProjections(cmd_line_args_.getEpsilon(),
                                              supported_,
                                              unsupported_,
                                              obj_1,
                                              obj_2);

        auto new_results = ResultContainer{};
        auto last_proj = nondom_projs.getLastProj();
        while (!nondom_projs.finished() && polyscip_status_ == PolyscipStatus::CompUnsupportedPhase) {
            auto left_proj = nondom_projs.getLeftProj();
            auto right_proj = nondom_projs.getRightProj();
            assert (left_proj.getFirst() < right_proj.getFirst());
            assert (left_proj.getSecond() > last_proj.getSecond());

            auto cons = vector<SCIP_CONS*>{};
            // create constraint pred.first <= c_{objs.first} \cdot x <= succ.first
            cons.push_back(createObjValCons(orig_vars[obj_1],
                                            orig_vals[obj_1],
                                            left_proj.getFirst(),
                                            right_proj.getFirst()));
            // create constraint optimal_val_objs.second <= c_{objs.second} \cdot x <= pred.second
            cons.push_back(createObjValCons(orig_vars[obj_2],
                                            orig_vals[obj_2],
                                            last_proj.getSecond(),
                                            left_proj.getSecond()));

            auto ref_point = std::make_pair(left_proj.getFirst() - 1., last_proj.getSecond() - 1.);
            // set beta = (beta_1,beta_2) s.t. pred and succ are both on the norm rectangle defined by beta
            auto beta_1 = 1.0;
            auto beta_2 = (right_proj.getFirst() - ref_point.first) / (left_proj.getSecond() - ref_point.second);
            // create constraint with respect to beta_1
            cons.push_back(createNewVarTransformCons(new_var,
                                                     orig_vars[obj_1],
                                                     orig_vals[obj_1],
                                                     ref_point.first,
                                                     beta_1));
            // create constraint with respect to beta_2
            cons.push_back(createNewVarTransformCons(new_var,
                                                     orig_vars[obj_2],
                                                     orig_vals[obj_2],
                                                     ref_point.second,
                                                     beta_2));
            for (auto c : cons) {
                SCIP_CALL(SCIPaddCons(scip_, c));
            }
            SCIP_CALL(solve());
            auto scip_status = SCIPgetStatus(scip_);
            if (scip_status == SCIP_STATUS_OPTIMAL) {
                assert (SCIPisGE(scip_, SCIPgetPrimalbound(scip_), 0.));
                auto res = getOptimalResult();
                auto proj = TwoDProj(res.second, obj_1, obj_2);

                if (left_proj.dominates(cmd_line_args_.getEpsilon(), proj) ||
                    right_proj.dominates(cmd_line_args_.getEpsilon(), proj)) {
                    nondom_projs.update();
                }
                else {
                    computeNondomResult(new_var,
                                        cons.front(),
                                        *std::next(begin(cons)),
                                        proj.getFirst(),
                                        proj.getSecond(),
                                        new_results);
                    auto nd_proj = TwoDProj(new_results.back().second, obj_1, obj_2);
                    nondom_projs.update(std::move(nd_proj), new_results.back().second);
                }
            }
            else if (scip_status == SCIP_STATUS_TIMELIMIT) {
                polyscip_status_ = PolyscipStatus::TimeLimitReached;
            }
            else {
                std::cout << "unexpected SCIP status in solveWeightedTchebycheff: " +
                             std::to_string(SCIPgetStatus(scip_)) + "\n";
                polyscip_status_ = PolyscipStatus::Error;
            }

            // release and delete constraints
            if (SCIPisTransformed(scip_))
                SCIP_CALL(SCIPfreeTransform(scip_));
            for (auto c : cons) {
                SCIP_CALL(SCIPdelCons(scip_, c));
                SCIP_CALL(SCIPreleaseCons(scip_, addressof(c)));
            }
        }
        if (no_objs_ > 3) {

            std::unique_ptr<Polyscip> subpoly( new Polyscip(cmd_line_args_,
                                                            scip_,
                                                            obj_sense_,
                                                            std::make_pair(obj_1, obj_2),
                                                            clock_total_) );
            subpoly->computeNondomPoints();
            subpoly->printResults();
            subpoly.reset();
        }
        for (const auto& res : new_results)
            unsupported_.push_back(res);
        return SCIP_OKAY;

    }



    /*SCIP_RETCODE Polyscip::solveWeightedTchebycheff(SCIP_VAR* new_var,
                                                    const vector<vector<SCIP_VAR*>>& orig_vars,
                                                    const vector<vector<ValueType>>& orig_vals,
                                                    const pair<size_t, size_t>& objs,
                                                    ValPairMap nondom_proj) {
        assert (!nondom_proj.empty());
        assert (orig_vars.size() == orig_vals.size());
        assert (orig_vals.size() == no_objs_);

        if (nondom_proj.size() > 1) {

            auto new_results = ResultContainer{};
            auto pred_it = nondom_proj.begin();

            std::cout << "SIZE = " << nondom_proj.size() << "\n";

            while (pred_it != std::prev(end(nondom_proj)) && polyscip_status_==PolyscipStatus::CompUnsupportedPhase) {

                auto pred = pred_it->first;
                auto succ_it = std::next(pred_it);
                auto succ = succ_it->first;
                auto last = std::prev(end(nondom_proj))->first;

                assert (pred.first < succ.first);
                assert (pred.second > last.second);

                auto cons = vector<SCIP_CONS *>{};
                // create constraint pred.first <= c_{objs.first} \cdot x <= succ.first
                cons.push_back(createObjValCons(orig_vars[objs.first],
                                                        orig_vals[objs.first],
                                                        pred.first,
                                                        succ.first));
                // create constraint optimal_val_objs.second <= c_{objs.second} \cdot x <= pred.second
                cons.push_back(createObjValCons(orig_vars[objs.second],
                                                        orig_vals[objs.second],
                                                        last.second,
                                                        pred.second));

                auto ref_point = std::make_pair(pred.first - 1., last.second - 1.);
                // set beta = (beta_1,beta_2) s.t. pred and succ are both on the norm rectangle defined by beta
                auto beta_1 = 1.0;
                auto beta_2 = (succ.first - ref_point.first) / (pred.second - ref_point.second);
                // create constraint with respect to beta_1
                cons.push_back(createNewVarTransformCons(new_var,
                                                         orig_vars[objs.first],
                                                         orig_vals[objs.first],
                                                         ref_point.first,
                                                         beta_1));
                // create constraint with respect to beta_2
                cons.push_back(createNewVarTransformCons(new_var,
                                                         orig_vars[objs.second],
                                                         orig_vals[objs.second],
                                                         ref_point.second,
                                                         beta_2));
                for (auto c : cons) {
                    SCIP_CALL(SCIPaddCons(scip_, c));
                }

                SCIP_CALL(solve());
                auto scip_status = SCIPgetStatus(scip_);
                if (scip_status == SCIP_STATUS_OPTIMAL) {
                    assert (SCIPisGE(scip_, SCIPgetPrimalbound(scip_), 0.));
                    auto res = getOptimalResult();
                    auto proj = std::make_pair(res.second[objs.first], res.second[objs.second]);


                    if (lhsLessEqualrhs(pred, proj) || lhsLessEqualrhs(succ, proj)) {
                        ++pred_it;
                    }
                    else {
                        computeNondomResult(new_var,
                                         cons.front(),
                                         *std::next(begin(cons)),
                                         proj.first,
                                         proj.second,
                                         new_results);
                        auto new_projection = std::make_pair(new_results.back().second[objs.first],
                                                             new_results.back().second[objs.second]);

                        std::cout << "NEW PROJ = " << new_projection.first << ", " << new_projection.second << "\n";

                        auto new_elem_it = nondom_proj.insert({new_projection, vector<OutcomeType>{new_results.back().second}});
                        assert (new_elem_it.second);
                        assert (new_elem_it.first != begin(nondom_proj)); // begin(nondom_proj) should never be dominated
                        // delete projections which are dominated by new_projection
                        if (lhsLessEqualrhs(new_projection, pred)) {
                            auto del_end_it = pred_it;
                            while (lhsLessEqualrhs(new_projection, del_end_it->first))
                                        ++del_end_it;
                            assert (del_end_it != std::prev(end(nondom_proj))); // --end(nondom_proj) should never be dominated
                            nondom_proj.erase(pred_it, del_end_it);
                            pred_it = new_elem_it.first;
                        }
                        else if (lhsLessEqualrhs(new_projection, succ)) {
                            auto del_end_it = succ_it;
                            while (lhsLessEqualrhs(new_projection, del_end_it->first))
                                ++del_end_it;
                            assert (del_end_it != std::prev(end(nondom_proj))); // --end(nondom_proj) should never be dominated
                            nondom_proj.erase(succ_it, del_end_it);
                        }
                    }
                }
                else if (scip_status == SCIP_STATUS_TIMELIMIT) {
                    polyscip_status_ = PolyscipStatus::TimeLimitReached;
                }
                else {
                    std::cout << "unexpected SCIP status in solveWeightedTchebycheff: " +
                                 std::to_string(SCIPgetStatus(scip_)) + "\n";
                    polyscip_status_ = PolyscipStatus::Error;
                }

                // release and delete constraints
                if (SCIPisTransformed(scip_))
                    SCIP_CALL(SCIPfreeTransform(scip_));
                for (auto c : cons) {
                    SCIP_CALL(SCIPdelCons(scip_, c));
                    SCIP_CALL(SCIPreleaseCons(scip_, addressof(c)));
                }
            }
            for (const auto& res : new_results)
                unsupported_.push_back(res);
        }
        return SCIP_OKAY;
    }*/


    void Polyscip::deleteVarNameFromResult(SCIP_VAR* var, Result& res) const {
        string name = SCIPvarGetName(var);
        auto pos = std::find_if(begin(res.first),
                                 end(res.first),
                     [&name](const std::pair<std::string, ValueType>& p){return name == p.first;});
        if (pos != end(res.first)) {
            res.first.erase(pos);
        }
    }

    SCIP_RETCODE Polyscip::initWeightSpace() {
        SCIP_CALL( computeUnitWeightOutcomes() ); // computes optimal outcomes for all unit weights
        if (polyscip_status_ == PolyscipStatus::InitPhase) {
            if (supported_.empty()) {
                polyscip_status_ = PolyscipStatus::Finished; // all outcomes for unit weights are unbounded
            }
            else {
                auto v_rep = DDMethod(scip_, no_objs_, supported_, unbounded_);
                v_rep.computeVRep_Var1();
                std::cout << "No of initial vertices = " << v_rep.size() << "\n";
                std::cout << "Starting initializing WSP...";
                weight_space_poly_ = global::make_unique<WeightSpacePolyhedron>(scip_,
                                                                                no_objs_,
                                                                                v_rep.moveVRep(),
                                                                                v_rep.moveHRep());
                std::cout << "...finished.\n";
                assert (weight_space_poly_->hasValidSkeleton(no_objs_));
                polyscip_status_ = PolyscipStatus::WeightSpacePhase;
            }
        }
        return SCIP_OKAY;
    }



    SCIP_STATUS Polyscip::separateINFORUNBD(const WeightType& weight, bool with_presolving) {
        if (!with_presolving)
            SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, TRUE);
        auto zero_weight = WeightType(no_objs_, 0.);
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
        auto outcome = OutcomeType(no_objs_,0.);
        auto no_vars = SCIPgetNOrigVars(scip_);
        auto vars = SCIPgetOrigVars(scip_);
        auto objs_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        for (auto i=0; i<no_vars; ++i) {
            auto var_sol_val = outcome_is_bounded ? SCIPgetSolVal(scip_, primal_sol, vars[i]) :
                               SCIPgetPrimalRayVal(scip_, vars[i]);

            if (!SCIPisZero(scip_, var_sol_val)) {
                sol.emplace_back(SCIPvarGetName(vars[i]), var_sol_val);
                auto var_obj_vals = OutcomeType(no_objs_, 0.);
                for (size_t index=0; index!=no_objs_; ++index) {
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
            auto val = obj_probdata->getWeightedObjVal(vars[i], weight);
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
            if (cmd_line_args_.onlyExtremal() || SCIPgetNOrigContVars(scip_) == SCIPgetNOrigVars(scip_)) { //check whether there exists integer variables
                polyscip_status_ = PolyscipStatus::Finished;
            }
            else {
                polyscip_status_ = PolyscipStatus::CompUnsupportedPhase;
            }
        }
        return SCIP_OKAY;
    }

    void Polyscip::printResults(ostream &os) const {
        for (const auto& result : supported_) {
            if (cmd_line_args_.outputOutcomes())
                outputOutcome(result.second, os);
            if (cmd_line_args_.outputSols())
                printSol(result.first, os);
            os << "\n";
        }
        for (const auto& result : unbounded_) {
            if (cmd_line_args_.outputOutcomes())
                outputOutcome(result.second, os, "Ray = ");
            if (cmd_line_args_.outputSols())
                printSol(result.first, os);
            os << "\n";
        }
        for (const auto& result : unsupported_) {
            if (cmd_line_args_.outputOutcomes())
                outputOutcome(result.second, os);
            if (cmd_line_args_.outputSols())
                printSol(result.first, os);
            os << "\n";
        }
    }

    void Polyscip::printSol(const SolType& sol, ostream& os) const {
        for (const auto& elem : sol)
            os << elem.first << "=" << elem.second << " ";
    }

    void Polyscip::outputOutcome(const OutcomeType &outcome, std::ostream &os, const std::string desc) const {
        if (obj_sense_ == SCIP_OBJSENSE_MAXIMIZE) {
            global::print(outcome, desc + "[ ", "] ", os, std::negate<ValueType>());
        }
        else {
            global::print(outcome, desc + "[ ", "] ", os);
        }
    }

    /*void Polyscip::printRay(const OutcomeType& ray, ostream& os) const {
        if (obj_sense_ == SCIP_OBJSENSE_MAXIMIZE) {
            global::print(ray, "Ray = [ ", "]", os, std::negate<ValueType>());
        }
        else {
            global::print(ray, "Ray = [ ", "]", os);
        }
    }*/

    bool Polyscip::filenameIsOkay(const string& filename) {
        std::ifstream file(filename.c_str());
        return file.good();
    }

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
        assert (checked_obj >= 1 && checked_obj < obj_probdata->getNoObjs());

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
        no_objs_ = obj_probdata->getNoObjs();

        auto vars = SCIPgetOrigVars(scip_);
        auto begin_nonzeros = vector<int>(no_objs_, 0);
        for (size_t i = 0; i < no_objs_ - 1; ++i)
            begin_nonzeros[i + 1] = global::narrow_cast<int>(
                    begin_nonzeros[i] + obj_probdata->getNumberNonzeroCoeffs(i));

        auto obj_to_nonzero_inds = vector<vector<int> >{};
        auto obj_to_nonzero_vals = vector<vector<SCIP_Real> >{};
        for (size_t obj_ind = 0; obj_ind < no_objs_; ++obj_ind) {
            auto nonzero_vars = obj_probdata->getNonZeroCoeffVars(obj_ind);
            auto size = nonzero_vars.size();
            if (size == 0)
                throw std::runtime_error(std::to_string(obj_ind) + ". objective is zero objective!");
            auto nonzero_inds = vector<int>(size, 0);
            std::transform(begin(nonzero_vars),
                           end(nonzero_vars),
                           begin(nonzero_inds),
                           [](SCIP_VAR *var) { return SCIPvarGetProbindex(var); });
            std::sort(begin(nonzero_inds), end(nonzero_inds));

            auto nonzero_vals = vector<SCIP_Real>(size, 0.);
            std::transform(begin(nonzero_inds),
                           end(nonzero_inds),
                           begin(nonzero_vals),
                           [&](int var_ind) { return obj_probdata->getObjCoeff(vars[var_ind], obj_ind); });


            if (cmd_line_args_.beVerbose())
                printObjective(obj_ind, nonzero_inds, nonzero_vals);

            obj_to_nonzero_inds.push_back(std::move(nonzero_inds));  // nonzero_inds invalid from now on
            obj_to_nonzero_vals.push_back(std::move(nonzero_vals));  // nonzero_vals invalid from now on

            if (obj_ind > 0 && objIsRedundant(begin_nonzeros, // first objective is always non-redundant
                               obj_to_nonzero_inds,
                               obj_to_nonzero_vals,
                               obj_ind))
                throw std::runtime_error(std::to_string(obj_ind) + ". objective is non-negative linear combination of previous objectives! Only problems with non-redundant objectives will be solved.");
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
            cout << "Number of objectives: " << no_objs_ << "\n";
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
            solfs << supported_.size() + unbounded_.size() + no_objs_ << " " << no_objs_ + 1 << " rational\n";
            for (const auto& elem : supported_) {
                global::print(elem.second, "0 ", " -1\n", solfs);
            }
            for (const auto& elem : unbounded_) {
                global::print(elem.second, "0 ", " 0", solfs);
            }
            for (size_t i=0; i<no_objs_; ++i) {
                auto ineq = vector<unsigned>(no_objs_, 0);
                ineq[i] = 1;
                global::print(ineq, "0 ", " 0\n", solfs);
            }
            solfs << "end\n";
            solfs.close();
        }
        else
            cout << "ERROR writing vertex enumeration file\n.";
    }

    void Polyscip::writeResultsToFile() const {
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

    bool Polyscip::isDominatedOrEqual(ResultContainer::const_iterator it, ResultContainer::const_iterator beg_it, ResultContainer::const_iterator end_it) const {
        for (auto curr = beg_it; curr != end_it; ++curr) {
            if (it == curr)
                continue;
            else if (std::equal(begin(curr->second),
                                end(curr->second),
                                begin(it->second),
                                std::less_equal<ValueType>())) {
                return true;
            }
        }
        return false;
    }


    bool Polyscip::dominatedPointsFound() const {
        auto results = ResultContainer{};
        for (auto& res : supported_)
            results.push_back(res);
        for (auto& res : unsupported_)
            results.push_back(res);

        for (auto cur=begin(results); cur!=end(results); ++cur) {
            if (isDominatedOrEqual(cur, begin(results), end(results)))
                return true;
        }
        return false;
    }

    void Polyscip::deleteWeaklyNondomSupportedResults() {
        auto it = begin(supported_);
        while (it != end(supported_)) {
            if (isDominatedOrEqual(it, begin(supported_), end(supported_))) {
                it = supported_.erase(it);
            }
            else {
                ++it;
            }
        }
    }




}
