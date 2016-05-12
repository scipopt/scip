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
#include <ostream>
#include <stdexcept>
#include <string>

#include "own_make_unique.h"
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "cmd_line_args.h"
#include "polyscip_types.h"
//ToDo change ReaderMOP
#include "ReaderMOP.h"

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
            throw std::domain_error("Invalid time limit.\n");
        if (cmd_line_args_.hasParameterFile() && !filenameIsOkay(cmd_line_args_.getParameterFile()))
            throw std::invalid_argument("Invalid parameter settings file.\n");
        if (!filenameIsOkay(cmd_line_args_.getProblemFile()))
            throw std::invalid_argument("Invalid problem file.\n");

        SCIPcreate(&scip_);
        assert (scip_ != nullptr);
        SCIPincludeDefaultPlugins(scip_);
        SCIPincludeObjReader(scip_, new ReaderMOP(scip_), TRUE);
    }

    void Polyscip::initWeightSpace(const OutcomeType& point,
                                   const std::vector<OutcomeType> &initial_rays,
                                   std::pair<bool, unsigned> unit_weight_info) {
        weight_space_poly_ = own_stl::make_unique<WeightSpacePolyhedron>(no_objs_, point, initial_rays, unit_weight_info);
    }

    Polyscip::~Polyscip() {
        SCIPfree(&scip_);
    }

    void Polyscip::printPoint(const OutcomeType& point, ostream& os) {
        print(point, {"Point = "}, os);
    }

    void Polyscip::printRay(const OutcomeType& ray, ostream& os) {
        print(ray, {"Ray = "}, os);
    }

    void Polyscip::printWeight(const WeightType& weight, ostream& os) {
        print(weight, {"Weight = "}, os);
    }

    bool Polyscip::filenameIsOkay(const string& filename) {
        std::ifstream file(filename.c_str());
        return file.good();
    }


//    SCIP_RETCODE Polyscip::readParamSettings() {
//        auto param_file = cmd_line_args_.getParameterSettingsFile();
//        SCIP_CALL( SCIPreadParams(scip_, param_file.c_str()) );
//        return SCIP_OKAY;
//    }

}
