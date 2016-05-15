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
#include <cmath> //std::abs
#include <fstream>
#include <functional> //std::plus
#include <iostream>
#include <limits>
#include <ostream>
#include <memory> //std::addressof
#include <numeric> //std::inner_product
#include <stdexcept>
#include <string>
#include <tuple>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "cmd_line_args.h"
#include "global_functions.h"
#include "polyscip_types.h"
#include "prob_data_objectives.h"
#include "ReaderMOP.h"

using std::addressof;
using std::cout;
using std::get;
using std::ostream;
using std::string;
using polyscip::OutcomeType;
using polyscip::ValueType;

namespace polyscip {

    Polyscip::Polyscip(int argc, const char *const *argv)
            : cmd_line_args_(argc, argv),
              polyscip_status_(PolyscipStatus::Unsolved),
              scip_(nullptr),
              obj_sense_(SCIP_OBJSENSE_MINIMIZE), // default objective sense is minimization
              no_objs_(0),
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
        for (auto& result : supported_) {
            auto sol = get<global::toField(ResultField::Solution)>(result);
            SCIPfreeSol(scip_, addressof(sol));
        }
        for (auto& result : unsupported_) {
            auto sol = get<global::toField(ResultField::Solution)>(result);
            SCIPfreeSol(scip_, addressof(sol));
        }
        for (auto& result : unbounded_) {
            auto sol = get<global::toField(ResultField::Solution)>(result);
            SCIPfreeSol(scip_, addressof(sol));
        }
        SCIPfreeClock(scip_, addressof(clock_total_));
        SCIPfree(addressof(scip_));
    }


    SCIP_RETCODE Polyscip::computeNondomPoints() {
        SCIP_CALL( computeSupported() );
        if (!cmd_line_args_.withUnsupported())
            std::cerr << "NOT IMPLEMENTED.\n";
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::initWeightSpace() {
        SCIP_CALL( SCIPstartClock(scip_, clock_total_) );
        decltype(no_objs_) obj_counter{0}; // objCounter has same type as no_objs
        auto initial_weight = WeightType(no_objs_,0.);

        while (polyscip_status_ == PolyscipStatus::Unsolved && obj_counter < no_objs_) {
            initial_weight[obj_counter] = 1.;
            SCIP_CALL( setWeightedObjective(initial_weight) );
            SCIP_CALL( solve() );
            auto scip_status = SCIPgetStatus(scip_);
            if (scip_status == SCIP_STATUS_INFORUNBD) {
                scip_status = separateINFORUNBD(initial_weight);
            }
            SCIP_CALL( handleStatus(scip_status, initial_weight, std::numeric_limits<ValueType>::max()) );
            initial_weight[obj_counter] = 0.;
            ++obj_counter;
        }
        return SCIP_OKAY;
    }

    SCIP_STATUS Polyscip::separateINFORUNBD(const WeightType& weight) {
        auto zero_weight = WeightType(no_objs_,0.);
        setWeightedObjective(zero_weight);
        solve(); // re-compute with zero objective
        auto status = SCIPgetStatus(scip_);
        if (status == SCIP_STATUS_INFORUNBD) {
            throw std::runtime_error("INFORUNBD Status with zero objective.");
        }
        else if (status == SCIP_STATUS_OPTIMAL) { // previous problem was unbounded
            setWeightedObjective(weight); // re-set to previous weighted objective
            return SCIP_STATUS_UNBOUNDED;
        }
        else {
            return status;
        }
    }

    SCIP_RETCODE Polyscip::handleStatus(SCIP_STATUS status, const WeightType& weight,
                                        ValueType current_wov) {
        if (status == SCIP_STATUS_OPTIMAL) {
            std::cout << "OPTIMAL STATUS\n";
            SCIP_CALL( handleOptimalStatus(weight, current_wov) );
        }
        else if (status == SCIP_STATUS_UNBOUNDED) {
            std::cout << "UNBOUNDED STATUS\n";
            //todo
        }
        else if (status == SCIP_STATUS_INFORUNBD) {
            //todo
        }
        else if (status == SCIP_STATUS_INFEASIBLE) {
            polyscip_status_ = PolyscipStatus::Solved;
        }
        else if (status == SCIP_STATUS_TIMELIMIT) {
            polyscip_status_ = PolyscipStatus::TimeLimitReached;
        }
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::handleOptimalStatus(const WeightType& weight,
                                               ValueType current_wov) {
        SCIP_SOL* best_sol = SCIPgetBestSol(scip_);
        SCIP_SOL* finite_sol = nullptr;
        SCIP_Bool same_obj_val = FALSE;
        SCIP_CALL( SCIPcreateFiniteSolCopy(scip_, addressof(finite_sol), best_sol, addressof(same_obj_val)) );
        if (!same_obj_val) {
            auto diff = std::abs(SCIPgetSolOrigObj(scip_, best_sol) -
                                 SCIPgetSolOrigObj(scip_, finite_sol));
            if (diff > 1.0e-5) {
                std::cerr << "absolute value difference after calling SCIPcreateFiniteSolCopy: " << diff << "\n";
                SCIP_CALL( SCIPfreeSol(scip_, addressof(finite_sol)) );
                throw std::runtime_error("SCIPcreateFiniteSolCopy: unacceptable difference in objective values.");
            }
        }
        auto outcome = getOutcome(finite_sol);
        auto new_wov = std::inner_product(begin(weight), end(weight), begin(outcome), 0.);
        if (new_wov < current_wov) {
            supported_.push_back(std::make_tuple(finite_sol, outcome, weight));
        }
        else {
            SCIP_CALL( SCIPfreeSol(scip_, addressof(finite_sol)) );
        }
        return SCIP_OKAY;
    }

    OutcomeType Polyscip::getOutcome(SCIP_SOL* sol) {
        auto objs = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        auto no_vars = SCIPgetNOrigVars(scip_);
        auto vars = SCIPgetOrigVars(scip_);
        auto point = OutcomeType(no_objs_,0.);
        for (auto i=0; i<no_vars; ++i) {
            auto current_var_sol_val = SCIPgetSolVal(scip_, sol, vars[i]);
            auto current_var_obj_vals = OutcomeType(no_objs_,0.);
            for (decltype(no_objs_) j=0; j<no_objs_; ++j)
                current_var_obj_vals[j] = objs->getObjVal(vars[i], j, current_var_sol_val);
            std::transform(begin(point), end(point),
                           begin(current_var_obj_vals),
                           begin(point),
                           std::plus<ValueType>());
        }
        return point;
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
        return SCIP_OKAY;
    }

    void Polyscip::printPoint(const OutcomeType& point, ostream& os) {
        global::print(point, {"Point = "}, os);
    }

    void Polyscip::printRay(const OutcomeType& ray, ostream& os) {
        global::print(ray, {"Ray = "}, os);
    }

    void Polyscip::printWeight(const WeightType& weight, ostream& os) {
        global::print(weight, {"Weight = "}, os);
    }

    bool Polyscip::filenameIsOkay(const string& filename) {
        std::ifstream file(filename.c_str());
        return file.good();
    }

    SCIP_RETCODE Polyscip::readProblem() {
        auto filename = cmd_line_args_.getProblemFile();
        SCIP_CALL( SCIPreadProb(scip_, filename.c_str(), "mop") );
        auto obj_probdata = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        assert (obj_probdata != nullptr);
        no_objs_ = obj_probdata->getNObjs();
        if (SCIPgetObjsense(scip_) == SCIP_OBJSENSE_MAXIMIZE) {
            obj_sense_ = SCIP_OBJSENSE_MAXIMIZE;
            // internally we treat problem as min problem and negate objective coefficients
            SCIPsetObjsense(scip_, SCIP_OBJSENSE_MINIMIZE);
            obj_probdata->negateAllCoeffs();
        }
        if (cmd_line_args_.beVerbose()) {
            cout << "No of objectives: " << no_objs_;
            cout << "\nObjective sense: ";
            if (obj_sense_ == SCIP_OBJSENSE_MAXIMIZE)
                cout << "MAXIMIZE\n";
            else
                cout << "MINIMIZE\n";
        }
        obj_probdata = nullptr;
        return SCIP_OKAY;
    }

}
