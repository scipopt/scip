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

#include "cmd_line_args.h"

#include <string>

#include "PolySCIPConfig.h" // defines EXECUTABLE_NAME, POLYSCIP_VERSION_{MAJOR,MINOR}
#include "tclap/CmdLine.h"

using std::string;
using TCLAP::CmdLine;
using TCLAP::SwitchArg;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;

namespace polyscip {

    CmdLineArgs::CmdLineArgs(int argc, char **argv)
            : executable_name_(EXECUTABLE_NAME)
    {
        version_no_ = std::to_string(POLYSCIP_VERSION_MAJOR) + string(".") + std::to_string(POLYSCIP_VERSION_MINOR);
        CmdLine cmd(executable_name_,' ', version_no_);
        //cmd.setExceptionHandling(false); // set internal exception handling
        SwitchArg compute_unsupported_arg("c", "computeUnsupported", "computes unsupported non-dominated points", true);
        cmd.add(compute_unsupported_arg);
        SwitchArg be_verbose_arg("v", "verbose", "verbose PolySCIP cmd line output ", false);
        cmd.add(be_verbose_arg);
        SwitchArg write_sols_arg("w", "writeSols", "write solutions to file", false);
        cmd.add(write_sols_arg);
        ValueArg<TimeLimitType> time_limit_arg("T", "timeLimit",
                                               "time limit in seconds for total computation time",
                                               false, 0, "seconds");
        cmd.add(time_limit_arg);
        ValueArg<string> write_sols_path_arg("W", "writeSolsPath",
                                             "PATH for -w; if only -w is given the default path is ./",
                                             false, "./", "PATH");
        cmd.add(write_sols_path_arg);
        ValueArg<string> param_file_arg("p", "paramSets", "parameter settings file for SCIP",
                                        false, "", "paramFile.set");
        cmd.add(param_file_arg);
        UnlabeledValueArg<string> prob_file_arg("probFile", "problem file in MOP format",
                                                true, "", "problemFile.mop");
        cmd.add(prob_file_arg);
        cmd.parse(argc, argv);

        compute_unsupported_ = compute_unsupported_arg.getValue();
        be_verbose_ = be_verbose_arg.getValue();
        write_sols_ = write_sols_arg.getValue();
        time_limit_ = time_limit_arg.getValue();
        write_sols_path_ = write_sols_path_arg.getValue();
        param_file_ = param_file_arg.getValue();
        prob_file_ = prob_file_arg.getValue();
    }

}