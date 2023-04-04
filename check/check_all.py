import argparse
import pathlib
from typing import List
from multiprocessing import Pool
import subprocess
import socket
import datetime
import os
import io

__is_local = socket.gethostname() == "TUE020355"


class Configuration:
    def __init__(self, binary : pathlib.Path, result : pathlib.Path,
                 testset : pathlib.Path, settings : str, instance : pathlib.Path, seed : int = 0,) -> None:
        self.binary : pathlib.Path = binary
        self.result : pathlib.Path = result
        self.testset : pathlib.Path = testset
        self.settings : pathlib.Path = settings
        self.instance : pathlib.Path = instance
        self.seed : int = seed

    def __repr__(self) -> str:
        return f"Configuration[{self.testset.stem},{self.settings.stem},{self.instance.stem},s{self.seed}]"

    def get_info(self):
        return [self.binary, self.result, self.testset, self.settings, self.instance, self.seed]

def run(conf : Configuration):
    print(f"{datetime.datetime.now()} | {os.getpid():5d} | Running {conf}")
    basename = f"result_{conf.testset.stem}_{conf.settings.stem}.{conf.instance.stem}.s{conf.seed}"
    outfile = conf.result / f"{basename}.out"
    errfile = conf.result / f"{basename}.err"
    setfile = conf.result / f"{basename}.set"

    with open(outfile, "w", buffering=io.DEFAULT_BUFFER_SIZE) as outf, \
            open(errfile, "w", buffering=io.DEFAULT_BUFFER_SIZE) as errf:
        # write configuration information
        outf.write(f"@-1 {datetime.datetime.now()}" + os.linesep)
        for idx, val in enumerate(conf.get_info()):
            outf.write(f"@{idx:02d} {val}" + os.linesep)
            errf.write(f"@{idx:02d} {val}" + os.linesep)
        outf.flush()
        errf.flush()

        args = []
        # slurm settings
        if __is_local:
            outf.write("@   Running without slurm (local tests)")
        else:
            outf.write("@   Running slurm")
            args.extend([
                "srun",
                "--cpu-bind=cores",
                "--ntasks=1",
                "--gres=cpu:4",
                "--cpu-freq=highm1",
                "--exclusive",
                "--cpus-per-task=4",
                "--mem-per-cpu=2000M",  # Not 2G (2048M) but 2000M to have some slack.
                "--time=01:30:00"  # Reserve half an hour of slack, in case shutting down takes more time.
            ])

        # run binary
        args.extend([
            conf.binary,
            "-c", f"set load {conf.settings}",
            "-c", "set limits time 3600",
            "-c", "set limits nodes 2100000000",
            "-c", "set limits memory 7360",  # 8GB = hardmemlimit = 1.1 * memlimit + 100, use 7360MB as soft limit.
            "-c", "set timing clocktype 1",
            "-c", "set display freq 10000",
            "-c", f"set save {setfile}",
            "-c", f"read {conf.instance}",
            "-c", "optimize",
            "-c", "display statistics",
            "-c", "q"
        ])

        # run the process, write stdout and stderr to the out and err-files
        status = subprocess.run(
            args=tuple(args),
            stdout=outf, stderr=errf, bufsize=io.DEFAULT_BUFFER_SIZE
        )
        outf.flush()
        errf.flush()

        # write finished information
        outf.write(os.linesep)
        errf.write(os.linesep)
        outf.write(f"@{idx+1:02d} {datetime.datetime.now()}" + os.linesep)
        errf.write(f"@{idx+1:02d} {datetime.datetime.now()}" + os.linesep)
        outf.flush()
        errf.flush()

    print(f"{datetime.datetime.now()} | {os.getpid():5d} | Finished {conf} with status {status.returncode}")
    try:
        status.check_returncode()
    except subprocess.CalledProcessError as e:
        print(f"{datetime.datetime.now()} | {os.getpid():5d} | BEGIN error output")
        print(e.output)
        print(f"{datetime.datetime.now()} | {os.getpid():5d} | END error output")


def extend_configurations(configurations, bin_path, settings, testsets, seed):
    testset : pathlib.Path
    instance : pathlib.Path
    setting : pathlib.Path
    for testset in testsets:
        with open(testset, "r") as f:
            for instance_str in f:
                instance_nice = instance_str.strip()
                if not instance_nice:
                    continue
                instance = check_folder / pathlib.Path(instance_str.strip())
                if not instance.exists():
                    print(f"Instance {instance.name} of {testset.name} does not exist.")
                    print(f"  {testset.name} -> {instance}")
                    continue
                for setting in settings:
                    configurations.append(
                        Configuration(
                            binary=bin_path, result=results_folder,
                            testset=testset, settings=setting, instance=instance, seed=seed
                        )
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Run testset")
    parser.add_argument("--ncores", type=int, default=4)
    parser.add_argument("--noise_dosage", type=bool, action=argparse.BooleanOptionalAction)
    parser.add_argument("--covering_designs", type=bool, action=argparse.BooleanOptionalAction)
    parser.add_argument("--miplib", type=bool, action=argparse.BooleanOptionalAction)
    parser.add_argument("--minlplib", type=bool, action=argparse.BooleanOptionalAction)
    parser.add_argument("--split_num", type=int, default=1, help="split the configurations over multiple runs")
    parser.add_argument("--split_id", type=int, default=0, help="split ids are from 0 to split_num")
    args = parser.parse_args()

    check_folder = pathlib.Path("check/")
    assert check_folder.exists(), "check folder does not exist"

    results_folder = check_folder / pathlib.Path("results")
    results_folder.mkdir(exist_ok=True)

    # Determine configurations
    configurations : List[Configuration] = []

    base_path_unified : pathlib.Path
    base_path_master : pathlib.Path
    if socket.gethostname() == "TUE020355":
        base_path_master = base_path_unified = pathlib.Path(".")
    else:
        base_path_unified = pathlib.Path("/home/mjohannesva/scip-unified")
        base_path_master = pathlib.Path("/home/mjohannesva/scip-8443db21")

    if args.noise_dosage:
        for seed in range(3):
            testsets = [
                base_path_unified / check_folder / pathlib.Path("testset/noise_dosage.test")
            ]

            # Master runs
            bin_path = base_path_master / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_nosym.set")
                # Noise dosage instances are nonbinary, so the standard settings cannot hold.
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)

            # Unified runs
            bin_path = base_path_unified / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_dyno.set"),
                base_path_master / pathlib.Path("settings/perf_dyno_first.set")
                # Noise dosage only has orbitopes -> only 'dyno' is considered.
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)

            # Sherali SHCs
            testsets = [
                base_path_unified / check_folder / pathlib.Path("testset/noise_dosage_sherali.test")
            ]
            bin_path = base_path_master / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_nosym.set")
                # Completely remove all symmetries from the problem. So, just run nosym.
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)


    if args.covering_designs:
        for seed in range(3):
            testsets = [
                base_path_unified / check_folder / pathlib.Path("testset/covering_designs_benchmark.test")
            ]

            # Master runs
            bin_path = base_path_master / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_nosym.set")
                # Covering design instances are nonbinary, so the standard settings cannot yield reductions
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)

            # Unified runs
            bin_path = base_path_unified / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_dynof.set"),
                base_path_master / pathlib.Path("settings/perf_dynofs.set"),
                base_path_master / pathlib.Path("settings/perf_dyns.set")
                # Covering designs are not orbitopal.
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)


    if args.miplib:
        for seed in range(1):  # only 1 seed!!
            testsets = [
                base_path_unified / check_folder / pathlib.Path("testset/miplib_benchmark_symmetries.test")
            ]

            # Master runs
            bin_path = base_path_master / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_nosym.set"),
                base_path_master / pathlib.Path("settings/perf_us1.set"),
                base_path_master / pathlib.Path("settings/perf_us2.set"),
                base_path_master / pathlib.Path("settings/perf_us3.set")
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)

            # Unified runs
            bin_path = base_path_unified / pathlib.Path("bin/scip")
            settings = [
                base_path_unified / pathlib.Path("settings/perf_dyno.set"),
                base_path_unified / pathlib.Path("settings/perf_dynof.set"),
                base_path_unified / pathlib.Path("settings/perf_dyns.set"),
                base_path_unified / pathlib.Path("settings/perf_dynofo.set"),
                base_path_unified / pathlib.Path("settings/perf_dynos.set"),
                base_path_unified / pathlib.Path("settings/perf_dynofs.set"),
                base_path_unified / pathlib.Path("settings/perf_dynofos.set")
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)


    if args.minlplib:
        raise NotImplementedError()

    # Split-run configuration. We run only for configurations[k::nruns] for k in range(nruns).
    if args.split_num != 1:
        assert args.split_num > 1, "Can only split in positive numbers"
        assert 0 <= args.split_id < args.split_num, "Split ID not in range"
        print(f"# Splitting configurations. ID {args.split_id}/{args.split_num}")
        configurations = list(configurations[args.split_id::args.split_num])

    # print configuration list
    print("# Configurations")
    for conf in configurations:
        print(conf)
    print(f"Number: {len(configurations)}")
    print()

    if __is_local:
        print("# Detected local environment, run without slurm")
    else:
        print("# Using slurm commands")
    print()

    # start processes
    with Pool(processes=(args.ncores // 4) - 1) as pool:  # The hypervising process blocks other runs, so maximally 31.
        pool.map(run, configurations)
