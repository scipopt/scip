/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file cmd_line_args.h
 * @brief  PolySCIP command line arguments
 * @author Sebastian Schenker
 *
 */

#ifndef POLYSCIP_SRC_CMD_LINE_ARGS_H_INCLUDED
#define POLYSCIP_SRC_CMD_LINE_ARGS_H_INCLUDED

#include <string>
#include <limits>

namespace polyscip {

    /**
     * @class CmdLineArgs
     * @brief Command line arguments for PolySCIP
     */
    class CmdLineArgs {
    public:
        /// Value type used for cmd line parameter -t | - - timeLimit
        using TimeLimitType = long;

        /// Constant used for default value of cmd line parameter -t | - - timeLimit
        constexpr static auto kTimeLimitInf = std::numeric_limits<TimeLimitType>::max();

        /**
         * Constructor
         * @param argc Argument count
         * @param argv Argument vector
         */
        CmdLineArgs(int argc, const char *const *argv);

        /**
         * For cmd line parameter -v | - -verbose
         * @return true if -v was set on cmd line
         */
        bool beVerbose() const {return be_verbose_;};

        /**
         * For cmd line parameter -x | - -extremal
         * @return true if -x was set on cmd line
         */
        bool onlyExtremal() const {return only_extremal_;};

        /**
         * For cmd line parameter -w | - -writeResults
         * @return true if -w was set on cmd line
         */
        bool writeResults() const {return write_results_;};

        /**
         * For cmd line parameter -o | - -noOutcomes
         * @return false if -o was set on cmd line
         */
        bool outputOutcomes() const { return output_outcomes_;};

        /**
         * For cmd line parameter -s | - -noSolutions
         * @return false if -s was set on cmd line
         */
        bool outputSols() const {return output_solutions_;};

        /**
         * For cmd line parameter -t | - -timeLimit
         * @return true if -t "<seconds>" was set on cmd line
         */
        bool hasTimeLimit() const {return time_limit_ != kTimeLimitInf;}

        /**
         * For cmd line parameter -p | - -params
         * @return true if -p <paramFile.set> was set on cmd line
         */
        bool hasParameterFile() const {return !param_file_.empty();};

        /**
         * For cmd line parameter -r | - -round
         * @return corresponding integer if -r <5|10|15> was set on cmd line
         */
        int roundWeightedObjCoeff() const {return round_weighted_obj_coeff_;};

        /**
         * For cmd line parameter -t | - -timeLimit
         * @return corresponding seconds given for parameter -t "<seconds>"
         */
        TimeLimitType getTimeLimit() const {return time_limit_;};

        /**
         * For cmd line parameter -d | - -delta
         * @return corresponding double value given for -d "<double>"
         */
        double getDelta() const {return delta_;};

        /**
         * For cmd line parameter -e | - -epsilon
         * @return corresponding double value given for -e "<double>"
         */
        double getEpsilon() const {return epsilon_;};

        /**
         * For cmd line parameter -p | - -params
         * @return corresponding parameter file given for -p "<paramFile.set>"
         */
        std::string getParameterFile() const {return param_file_;};

        /**
         * For cmd line argument problemFile.mop
         * @return corresponding problem file given for "<problemFile.mop>" on cmd line
         */
        std::string getProblemFile() const {return prob_file_;};

        /**
         * For cmd line parameter -W | - -writeSolsPath
         * @return corresponding path given -W "<PATH>"
         */
        std::string getWritePath() const {return write_results_path_;};

    private:
        std::string executable_name_; ///< Name of the executable binary
        std::string version_no_; ///< Current version number of PolySCIP
        bool be_verbose_; ///< -v parameter value
        bool only_extremal_; ///< -e parameter value
        bool write_results_; ///< -w paramter value
        bool output_solutions_; ///< -s paramter value
        bool output_outcomes_; ///< -o parameter value
        int round_weighted_obj_coeff_; ///< -r parameter value
        TimeLimitType time_limit_; ///< -t parameter value
        double delta_; ///< -d parameter value
        double epsilon_; ///< -e parameter value
        std::string param_file_; ///< -p parameter value
        std::string prob_file_; ///< Name of the given problem file
        std::string write_results_path_; ///< Path where result file is written to
    };

}
#endif //POLYSCIP_SRC_CMD_LINE_ARGS_H_INCLUDED