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

#include <fstream>
#include <iostream>
#include <ostream>
#include <memory> //std::addressof
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

        SCIPcreate(addressof(scip_));
        SCIPincludeDefaultPlugins(scip_);
        SCIPincludeObjReader(scip_, new ReaderMOP(scip_), TRUE);
        if (cmd_line_args_.hasParameterFile())
            SCIPreadParams(scip_, cmd_line_args_.getParameterFile().c_str());
        SCIPcreateClock(scip_, addressof(clock_total_));
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
        SCIPfree(addressof(scip_));
    }


    void Polyscip::computeNondomPoints() {
        computeSupported();
        if (!cmd_line_args_.withUnsupported())
            computeUnsupported();
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
            if (scip_status == SCIP_STATUS_INFORUNBD)
                scip_status = separateINFORUNBD(initial_weight);
            SCIP_CALL( handleStatus(scip_status) );
            //todo
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

    SCIP_RETCODE Polyscip::handleStatus(SCIP_STATUS status) {
        if (status == SCIP_STATUS_OPTIMAL) {
            //todo
        }
        else if (status == SCIP_STATUS_UNBOUNDED) {
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

    SCIP_RETCODE Polyscip::solve() {
        if (cmd_line_args_.hasTimeLimit()) { // set SCIP timelimit
            auto remaining_time_limit = cmd_line_args_.getTimeLimit() -
                                        SCIPgetClockTime(scip_, clock_total_);
            SCIP_CALL(SCIPsetRealParam(scip_, "limits/time", remaining_time_limit));
        }
        SCIP_CALL( SCIPsolve(scip_) );    // actual SCIP solver call
        return SCIP_OKAY;
    }

    SCIP_RETCODE Polyscip::setWeightedObjective(const WeightType& weight){
        if (SCIPisTransformed(scip_))
            SCIP_CALL( SCIPfreeTransform(scip_) );
        ProbDataObjectives* objs = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        auto vars = SCIPgetOrigVars(scip_);
        auto no_vars = SCIPgetNOrigVars(scip_);
        for (auto i=0; i<no_vars; ++i) {
            auto val = objs->getWeightedObjVal(vars[i], weight);
            SCIP_CALL( SCIPchgVarObj(scip_, vars[i], val) );
        }
        return SCIP_OKAY;
    }

    void Polyscip::computeSupported() {
        initWeightSpace();
        //TODO
    }

    void Polyscip::computeUnsupported() {
        //TODO
        ;
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
        ProbDataObjectives* objs = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
        no_objs_ = objs->getNObjs();
        if (SCIPgetObjsense(scip_) == SCIP_OBJSENSE_MAXIMIZE) {
            obj_sense_ = SCIP_OBJSENSE_MAXIMIZE;
            // internally we treat problem as min problem and negate objective coefficients
            SCIPsetObjsense(scip_, SCIP_OBJSENSE_MINIMIZE);
            objs->negateAllCoeffs();
        }
        if (cmd_line_args_.beVerbose()) {
            cout << "No of objectives: " << no_objs_;
            cout << "\nObjective sense: ";
            if (obj_sense_ == SCIP_OBJSENSE_MAXIMIZE)
                cout << "MAXIMIZE\n";
            else
                cout << "MINIMIZE\n";
        }
        objs = nullptr;
        return SCIP_OKAY;
    }

}
