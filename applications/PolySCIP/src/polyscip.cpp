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

#include "own_make_unique.h"
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "cmd_line_args.h"
#include "global_functions.h"
#include "polyscip_types.h"
//ToDo change ReaderMOP und ProbDataObjs
#include "ProbDataObjectives.h"
#include "ReaderMOP.h"

using std::addressof;
using std::cout;
using std::get;
using std::ostream;
using std::string;

namespace polyscip {

    Polyscip::Polyscip(int argc, const char *const *argv)
            : cmd_line_args_(argc, argv),
              scip_(nullptr),
              obj_sense_(SCIP_OBJSENSE_MINIMIZE), // default objective sense is minimization
              no_objs_(0)
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
        SCIPcreateClock(scip_, addressof(clock_iteration_));
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

    bool Polyscip::initWeightSpace() {
        decltype(no_objs_) obj_counter{0}; // objCounter has same type as no_objs
        auto initial_weight = WeightType(no_objs_,0.);
        SCIP_CALL( SCIPstartClock(scip_, clock_total_) );
        auto found_point = false;
        while (obj_counter < no_objs_) {
            SCIP_CALL( restartClockIteration() );
            //todo
            ++obj_counter;
        }

        return false;
    }

    SCIP_RETCODE Polyscip::restartClockIteration() {
        SCIP_CALL( SCIPresetClock(scip_, clock_iteration_) );
        SCIP_CALL( SCIPstartClock(scip_, clock_iteration_) );
        return SCIP_OKAY;
    }

    void Polyscip::computeSupported() {
        auto point_found = initWeightSpace();
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
        no_objs_ = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_))->getNObjs();
        if (SCIPgetObjsense(scip_) == SCIP_OBJSENSE_MAXIMIZE) { // objective sense of read problem is maximization
            obj_sense_ = SCIP_OBJSENSE_MAXIMIZE;
            //internally we treat problem as min problem and negate objective values
            SCIPsetObjsense(scip_, SCIP_OBJSENSE_MINIMIZE);
            dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_))->negateAllObjCoeffs();
        }
        if (cmd_line_args_.beVerbose()) {
            cout << "No of objectives: " << no_objs_;
            cout << "\nObjective sense: ";
            if (obj_sense_ == SCIP_OBJSENSE_MAXIMIZE)
                cout << "MAXIMIZE\n";
            else
                cout << "MINIMIZE\n";

        }
        return SCIP_OKAY;
    }

}
