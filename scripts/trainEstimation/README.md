How to train custom tree size estimation for SCIP                         {#TRAINESTIMATION}
=================================================

While watching SCIP's periodic output, users may sometimes wonder to what extent the
search process is already finished.
As of version 7.0, SCIP features a new display column "compl." that, in many cases, serves as
a much better progress bar of the search than the traditional gap.
This column tries to approximate the fraction of completed nodes, the so-called _search completion_,
as closely as possible.
By default, the display is based on simple search tree and gap statistics that are collected during
the search.
While most of these statistics by themselves show a significant improvement over the classical gap,
they can also be combined into even more accurate approximations of the search completion.

SCIP provides two ways to combine the individual statistics:

1. a linear regression that uses the two values "tree weight" and "SSG" and is guaranteed to be monotone.
2. a regression forest on all tree statistics

Especially the second method requires **careful training to the instances of interest**. Therefore,
SCIP does not come with a built-in regression forest, but can read custom regression forests
that have been trained on user data.

This tutorial shows how to use the scripts in the directory "scripts/trainEstimation/" of SCIP
to train and use custom regressions from user data.


Prerequisites: Installation of necessary R packages
---------------------------------------------------

The training is performed by the R script "train.R". It depends on the availability of some R packages.
Please make sure that you have the following packages installed, ideally in their newest versions.

- readr
- magrittr
- ggplot2
- dplyr
- knitr
- rpart
- randomForest
- reshape2

You can install packages from within an R session interactively by executing

```{.R}
install.packages(
    c(
        "readr",
        "magrittr",
        "ggplot2",
        "dplyr",
        "knitr",
        "rpart",
        "randomForest",
        "reshape2"
    )
)
```

You will be prompted to specify a CRAN-mirror before downloading the package resources.
The above packages are well-maintained standard packages in the R universe. If you
experience errors during the installation, please refer to the troubleshooting
in the respective package documentation.

If the installation was succesful, you should be able to test the scripts on our provided test data:

```
./run_training.sh testdata/
```

Th directory "testdata/" contains only a handful of example log files to verify that the R packages have been set up successfully.
For good training results in a practical scenario, it is recommended to provide at least 50-100 such log files.
Smaller test beds can be enriched, for example, by running SCIP multiple times per instance with different random seed initializations.

At successful termination, the training summarizes the training in several new files in the output directory "output/".
The output of the script is explained further below.


Creation of SCIP Log File Output
--------------------------------

The second step consists of producing meaningful training data in the form of SCIP output on instances of interest.
The required additional output can be enabled using the settings file "periodic_report.set" in this directory for SCIP.
The Log files must be stored in a common directory used as argument for "run_training.sh", one log file per instance.
Log files must have the file extension ".out".

The following example shows how to create a new subdirectory "mydata/" and produce a log file for the instance "bell5.mps"
under Linux.

```
mkdir mydata/
scip -s periodic_report.set -f ../../check/instances/MIP/bell5.mps -l mydata/scip-logfile.out
```

The training only considers instances that could be solved, and discards all instances with trees that are too small.
Therefore, it should be ensured that at the end of the data collection, there are instances with interesting trees not too small
that could be completely solved.

If the data set consists of many instances, make sure to also read [how to run automated tests with SCIP](@ref TEST)
of the SCIP documentation for one possibility to automate the data collection.
Unless the `OUTPUTDIR=results` flag is modified, the necessary log files are then collected in the subdirectory "check/results/"
under the SCIP root directory.
Please note that the above scripts create an additional log file after the runs are finished, in which
log files for each instance are concatenated.
Please (re-)move this file, which can be easily identified by the prefix "check." before proceeding.

Invocation of run_training.sh
-----------------------------

Now that the data set has been collected in a subdirectory "mydata/", we can invoke

```
./run_training.sh mydata/
```

to obtain information of the approximation/estimation accuracy of the different estimation methods of SCIP on the newly created data set.
All information is stored under the default subdirectory "output/".
The name of the out directory can be changed by providing it as second argument to run_training.sh

```
./run_training.sh mydata/ /path/to/my/output/
```

Interpretation of the Output
----------------------------

After the training has finished, a comparison of different tree size estimation methods is printed to the console. In addition,
the (user-specified) output directory contains the following files.

- monotone.set, a SCIP settings file with linear regression coefficients.
- rf_model.rfcsv, a trained regression forest
- searchcompletion_mse.pdf, a bar chart that compares the search completion accuracy.

In addition, some intermediate data pipeline results are also stored in the output directory in the form of CSV files.

### Output of run_training.sh

At termination, "run_training.sh" outputs the estimation accuracy for a total of 13 different techniques to estimate the tree size.
An example table could look as follows:

|   |Method          |    n|   MSE| MeanRatio| 2Accurate| 3Accurate| 4Accurate|Group            |
|:--|:---------------|----:|-----:|---------:|---------:|---------:|---------:|:----------------|
|13 |Random.Forest   | 3012| 0.009|     1.540|     0.824|     0.883|     0.910|Learned          |
|12 |linear.monotone | 3012| 0.049|     2.074|     0.707|     0.822|     0.863|Learned          |
|1  |treeweight      | 3012| 0.074|     2.191|     0.686|     0.791|     0.841|SearchCompletion |
|6  |wbe             | 3012|   NaN|     2.349|     0.655|     0.777|     0.824|Custom           |
|11 |tree-weight     | 3012|   NaN|     2.395|     0.651|     0.759|     0.810|Forecast         |
|8  |leaf-frequency  | 3012|   NaN|     2.465|     0.666|     0.763|     0.807|Forecast         |
|9  |open-nodes      | 3012|   NaN|     2.716|     0.613|     0.714|     0.760|Forecast         |
|4  |leaf-frequency  | 3012| 0.151|     2.755|     0.600|     0.714|     0.754|SearchCompletion |
|10 |ssg             | 3012|   NaN|     3.016|     0.603|     0.698|     0.739|Forecast         |
|7  |gap             | 3012|   NaN|     3.308|     0.539|     0.658|     0.708|Forecast         |
|5  |tree-profile    | 3012|   NaN|     4.140|     0.463|     0.580|     0.638|Custom           |
|2  |ssg             | 3012| 0.068|     6.301|     0.631|     0.701|     0.758|SearchCompletion |
|3  |gap             | 3012| 0.259|     9.378|     0.497|     0.594|     0.652|SearchCompletion |

For four methods, there are two ways to compute an estimation of tree size: either as an approximation
of search completion, or by computing a time series forecast that takes into account the most recent values and trend development.
Therefore, these methods appear twice in the table, and the column "Group" shows whether a forecast or
a search completion approximation was used.
The column n contains the number of records that were aggregated for this table.

The methods are sorted by their geometric mean approximation ratio of the true tree size.
We normalize each ratio such that it is bounded from below by 1, which would correspond to a perfect estimate.
The method at the top, Random.Forest, is the method that yields
the most accurate estimate of search tree size at termination.
For methods that approximate search completion, the Mean Squared Error of the approximation is shown in column "MSE".
The columns 2Accurate etc. give the fraction of records that are within a factor of 2 (3,4) of the actual tree size at termination.
Better methods reach higher values in those columns.

Disclaimer: This table has been produced as a showcase on a subset of 91 publicly available instances almost all of which are solved by SCIP in less than 100 seconds.
The figures therein are not representative beyond the data set, and the ranking of the methods may substantially change on other data sources.


Using the Trained Regression in SCIP
------------------------------------

In the above table, the two learned methods "Random.Forest" and "linear.monotone" outperform the other methods (out of which they are constructed).
In the output directory, the file "monotone.set" contains the linear regression coefficients and can be input into SCIP like a normal settings file.
For example,

```
scip -s output/monotone.set
```

launches a SCIP interactive shell with the coefficients preloaded. Since the coefficients are user parameters,
SCIP comes with reasonable default values that have been calibrated on a large data set.
We call this a monotone regression because it combines the two values "tree weight" and "SSG", which are monotone.

The trained regression forest can also be loaded into SCIP. The file "rf_model.rfcsv" in the output directory
contains all learned splits on the input features for all trees in its ensemble as one
large data file in an extended CSV format with additional storage size information.
The location of this file must be explicitly specified for SCIP by setting the string parameter "estimation/regforestfilename",
e.g., via

```
estimation/regforestfilename = "/path/to/rf_model.rfcsv"
```
as part of settings file for SCIP.
In the interactive shell, type `set estimation regforestfilename output/rf_model.rfcsv` to tell SCIP
that it should load the regression forest.
A loaded regression forest automatically takes precedence over all other methods and is used as search completion approximation in the new
display column "compl.".
In contrast to the monotone regression, the regression forest is not guaranteed to be monotone.

Note that most of the other methods in the table are also readily available. If the tree size should be estimated by an SSG forecasting
method, use the "method" parameter

```
estimation/method = s
```

If the training suggests that tree weight search completion should be used instead, combine the two parameters

```
estimation/completiontype = w
estimation/method = c
```
