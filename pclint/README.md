How to run pclint (for the first time)
---

In order for PC Lint Plus to be able to configure your compiler, create a virtual python environment and install the packages regex and pyyaml

    virtualenv -p python3 venv
    . venv/bin/activate
    pip install regex pyyaml

then run `make pclint` with the variable `PCLINTCONFIG` pointing to the file `pclp_config.py` of your PC Lint Plus installation

    make pclint PCLINTCONFIG="/path/to/pc-lint-plus/pclp_config.py"

This process will create two config files in the `pclint` folder.
After their creation you will be able to run just `make pclint`.

The file `pclint.out` will contain the test report.
