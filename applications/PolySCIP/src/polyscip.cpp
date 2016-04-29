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

#include <ostream>

using std::ostream;

namespace polyscip {

    void Polyscip::printPoint(const OutcomeType& point, ostream& os) {
        print(point, {"Point = "}, os);
    }

    void Polyscip::printRay(const OutcomeType& ray, ostream& os) {
        print(ray, {"Ray = "}, os);
    }

    void Polyscip::printWeight(const WeightType& weight, ostream& os) {
        print(weight, {"Weight = "}, os);
    }

}
