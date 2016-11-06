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

/** @brief  Command line arguments
 *
 * Class holding command line arguments.
 */

#ifndef POLYSCIP_SRC_CMD_LINE_ARGS_H_INCLUDED
#define POLYSCIP_SRC_CMD_LINE_ARGS_H_INCLUDED

#include <string>
#include <limits>

namespace polyscip {

    class CmdLineArgs {
    public:
        using TimeLimitType = long;
        constexpr static auto kTimeLimitInf = std::numeric_limits<TimeLimitType>::max();

        CmdLineArgs(int argc, const char *const *argv);

        bool beVerbose() const {return be_verbose_;};
        bool onlyExtremal() const {return only_extremal_;};
        bool writeResults() const {return write_results_;};
        bool outputOutcomes() const { return output_outcomes_;};
        bool outputSols() const {return output_solutions_;};
        bool hasTimeLimit() const {return time_limit_ != kTimeLimitInf;}
        bool hasParameterFile() const {return !param_file_.empty();};
        TimeLimitType getTimeLimit() const {return time_limit_;};
        double getBeta() const {return beta_;};
        double getDelta() const {return delta_;};
        double getEpsilon() const {return epsilon_;};
        std::string getParameterFile() const {return param_file_;};
        std::string getProblemFile() const {return prob_file_;};
        std::string getWritePath() const {return write_results_path_;};

    private:
        std::string executable_name_;
        std::string version_no_;

        // arguments read from command line
        bool be_verbose_;
        bool only_extremal_;
        bool write_results_;
        bool output_solutions_;
        bool output_outcomes_;
        TimeLimitType time_limit_;
        unsigned beta_;
        double delta_;
        double epsilon_;
        std::string param_file_;
        std::string prob_file_;
        std::string write_results_path_;
    };

}
#endif //POLYSCIP_SRC_CMD_LINE_ARGS_H_INCLUDED