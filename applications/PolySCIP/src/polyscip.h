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

#include <iostream>
#include <ostream>
#include <memory>
#include <string>
#include <utility> // std::pair
#include <vector>

#include "cmd_line_args.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_polyhedron.h"

namespace polyscip {

    class Polyscip {
    public:
        /**< container type for nondominated points */
        using PointContainer = std::vector <std::pair<OutcomeType, WeightType>>;
        /**< container type for nondominated rays; needs to support: empty()*/
        using RayContainer = std::vector<OutcomeType>;

        Polyscip(int argc, const char *const *argv);
        ~Polyscip();

        void computeSupportedNondomPoints() = delete;
        void computeUnSupportedNondomPoints() = delete;

        void initWeightSpace(const OutcomeType& point,
                             const std::vector<OutcomeType>& initial_rays,
                             std::pair<bool, unsigned> unit_weight_info);

        /** Prints given weight to given output stream
         */
        void printWeight(const WeightType& weight, std::ostream& os = std::cout);

        /** Prints given point to given output stream
         */
        void printPoint(const OutcomeType& point, std::ostream& os = std::cout);

        /** Prints given ray to given output stream
         */
        void printRay(const OutcomeType& ray, std::ostream& os = std::cout);

    private:
        bool filenameIsOkay(const std::string& filename);
        /** Reads SCIP parameter settings */
        SCIP_RETCODE readParamSettings() = delete;

        CmdLineArgs cmd_line_args_;
        SCIP* scip_;
        SCIP_Objsense obj_sense_;                      /**< objective sense of given problem */
        unsigned no_objs_;                             /**< number of objectives */
        std::unique_ptr<WeightSpacePolyhedron> weight_space_poly_;
        PointContainer supported_nondom_points_;
        PointContainer unsupported_nondom_points_;
        RayContainer unbounded_nondom_rays_;
    };

}

#endif //POLYSCIP_SRC_POLYSCIP_H_INCLUDED
