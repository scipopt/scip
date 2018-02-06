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
 * @file cmd_line_args.cpp
 * @brief Implements PolySCIP command line arguments via TCLAP
 * @author Sebastian Schenker
 *
 */


#include "cmd_line_args.h"

#include <string>
#include <vector>

#include "PolySCIPConfig.h" // defines EXECUTABLE_NAME, POLYSCIP_VERSION_{MAJOR,MINOR}
#include "tclap/CmdLine.h"

using std::string;
using TCLAP::CmdLine;
using TCLAP::SwitchArg;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;

namespace polyscip {

    /**
     * Constructor
     * @param argc Argument count
     * @param argv Argument vector
     */
    CmdLineArgs::CmdLineArgs(int argc, const char *const *argv)
            : executable_name_(EXECUTABLE_NAME)
    {
        version_no_ = std::to_string(POLYSCIP_VERSION_MAJOR) + string(".") + std::to_string(POLYSCIP_VERSION_MINOR);
        CmdLine cmd(executable_name_,' ', version_no_);
        cmd.setExceptionHandling(false); // set internal exception handling
        SwitchArg only_extremal_arg("x", "extremal", "compute only extremal supported non-dominated results", false);
        cmd.add(only_extremal_arg);
        SwitchArg be_verbose_arg("v", "verbose", "verbose PolySCIP cmd line output ", false);
        cmd.add(be_verbose_arg);
        SwitchArg write_results_arg("w", "writeResults", "write results to file; default path is ./", false);
        cmd.add(write_results_arg);
        SwitchArg output_sols_arg("s", "noSolutions", "switching output of solutions off", true);
        cmd.add(output_sols_arg);
        SwitchArg output_outcomes_arg("o", "noOutcomes", "switching output of outcomes off", true);
        cmd.add(output_outcomes_arg);
        std::vector<int> round_vals = {5,10,15};
        TCLAP::ValuesConstraint<int> allowed_vals(round_vals);
        ValueArg<int> round_weighted_obj_coeff_arg("r", "round",
                                          "round weighted objective coefficient in fct 'setWeightedObjective' at r-th decimal position",
                                          false, 0, &allowed_vals);
        cmd.add(round_weighted_obj_coeff_arg);
        ValueArg<TimeLimitType> time_limit_arg("t", "timeLimit",
                                               "time limit in seconds for total SCIP computation time",
                                               false, kTimeLimitInf, "seconds");
        cmd.add(time_limit_arg);
        ValueArg<double> delta_arg("d", "delta", "Delta used in computation of feasible boxes; default value: 0.01",
                                     false, 0.01, "double");
        cmd.add(delta_arg);
        ValueArg<double> epsilon_arg("e", "epsilon", "epsilon used in computation of unsupported points; default value: 1e-3",
                                     false, 1e-3, "double");
        cmd.add(epsilon_arg);
        ValueArg<string> write_sols_path_arg("W", "writeSolsPath",
                                             "PATH for -w",
                                             false, "./", "PATH");
        cmd.add(write_sols_path_arg);
        ValueArg<string> param_file_arg("p", "params", "parameter settings file for SCIP",
                                        false, "", "paramFile.set");
        cmd.add(param_file_arg);
        UnlabeledValueArg<string> prob_file_arg("probFile", "problem file in MOP format",
                                                true, "", "problemFile.mop");
        cmd.add(prob_file_arg);
        cmd.parse(argc, argv);

        be_verbose_ = be_verbose_arg.getValue();
        only_extremal_ = only_extremal_arg.getValue();
        write_results_ = write_results_arg.getValue();
        output_solutions_ = output_sols_arg.getValue();
        output_outcomes_ = output_outcomes_arg.getValue();
        round_weighted_obj_coeff_ = round_weighted_obj_coeff_arg.getValue();
        time_limit_ = time_limit_arg.getValue();
        delta_ = delta_arg.getValue();
        epsilon_ = epsilon_arg.getValue();
        write_results_path_ = write_sols_path_arg.getValue();
        param_file_ = param_file_arg.getValue();
        prob_file_ = prob_file_arg.getValue();
    }

}