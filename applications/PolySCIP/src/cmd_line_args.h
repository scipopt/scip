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

        bool checkForRedundantObjs() const { return check_for_redundant_objs_;};
        bool beVerbose() const {return be_verbose_;};
        bool withUnsupported() const {return with_unsupported_;};
        bool writeSolutions() const {return write_sols_;};
        bool hasTimeLimit() const {return time_limit_ != kTimeLimitInf;}
        bool hasParameterFile() const {return !param_file_.empty();};
        TimeLimitType getTimeLimit() const {return time_limit_;};
        double getEpsilon() const {return epsilon_;};
        std::string getParameterFile() const {return param_file_;};
        std::string getProblemFile() const {return prob_file_;};
        std::string getWritePath() const {return write_sols_path_;};

    private:
        std::string executable_name_;
        std::string version_no_;

        // arguments read from command line
        bool check_for_redundant_objs_;
        bool be_verbose_;
        bool with_unsupported_;
        bool write_sols_;
        TimeLimitType time_limit_;
        double epsilon_;
        std::string param_file_;
        std::string prob_file_;
        std::string write_sols_path_;
    };

}
#endif //POLYSCIP_SRC_CMD_LINE_ARGS_H_INCLUDED