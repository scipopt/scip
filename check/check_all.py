import argparse
import pathlib
from typing import List
from multiprocessing import Pool
import subprocess
import socket


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

def run(conf : Configuration):
    print(f"Running {conf}")
    basename = f"result_{conf.testset.stem}_{conf.settings.stem}.{conf.instance.stem}.s{conf.seed}"
    outfile = conf.result / f"{basename}.out"
    errfile = conf.result / f"{basename}.err"
    setfile = conf.result / f"{basename}.set"

    with open(outfile, "w") as outf, open(errfile, "w") as errf:
        status = subprocess.run(
            args=(
                # # slurm settings
                # "srun", "--cpu-bind=cores",
                # "--ntasks=1",
                # "--gres=cpu:1",
                # "--cpu-freq=highm1",
                # "--exclusive",
                # run binary
                conf.binary,
                "-c", f"set load {conf.settings}",
                "-c", "set limits time 7200",
                "-c", "set limits nodes 2100000000",
                "-c", "set limits memory 4096",
                "-c", "set timing clocktype 1",
                "-c", "set display freq 10000",
                "-c", f"set save {setfile}",
                "-c", f"read {conf.instance}",
                "-c", "optimize",
                "-c", "display statistics",
                "-c", "q"
            ),
            stdout=outf, stderr=errf
        )

    print(f"Finished {conf} with status {status.returncode}")


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
    parser.add_argument("-ncores", type=int, default=4)
    parser.add_argument("-variant", type=str, default="noise_dosage")
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

    if args.variant == "noise_dosage":
        for seed in range(3):
            testsets = [
                base_path_unified / check_folder / pathlib.Path("testset/short.test")
            ]

            # Master runs
            bin_path = base_path_master / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_nosym.set"),
                base_path_master / pathlib.Path("settings/perf_us1.set")
                # Noise dosage only has orbitopes -> polyhedral; us2,3 are not needed.
            ]
            extend_configurations(configurations, bin_path, settings, testsets, seed)

            # Unified runs
            bin_path = base_path_unified / pathlib.Path("bin/scip")
            settings = [
                base_path_master / pathlib.Path("settings/perf_dyno.set")
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


    # print configuration list
    print("# Configurations")
    for conf in configurations:
        print(conf)
    print()

    # start processes
    with Pool(processes=args.ncores) as pool:
        pool.map(run, configurations)
