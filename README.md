# FRHPJRSSC
All of the scripts (and data) used in "A Fractional Model for Earthquakes" submitted to JRSSC Series C

The fractional Hawkes process is now referred to as FHP.
Notes:
These are the data files (JT.csv, LM.csv HM.csv and KaikoraQuakes.csv) that were used. To replicate results in the manuscript, one needs to remove events in LM.csv that have magnitude <2.75 and then only use events 1 to 1348 of those left and only use events 1 to 1023 in HM.csv (parameters for this are in the input files and are done automatically when running the scripts).

Scripts within the "Consistency" file are what was used for the simulation study, along with the "ml.m" function and its C++ adaption as mentioned. 
"ml.m" and "MLapp.m" are needed for all of the MATLAB scripts to run, and is available in the "Consistency" file.
Info for some of this code below.
#

## Overview of Consistency

## Compiling the C++ code

Some parts of the code have been re-implemented in C++ and you may need to run
```
mex -setup
```
You will need to compile the C++ code prior to running the code. If you wish to enable OpenMP threading for fine grain parallelism,
then you will need to pass the appropriate OpenMP option
(-fopenmp below for the GNU compiler) as well as the OpenMP library to link against (-lgomp in the case below).
The full compilation step then becomes:
```
mex CXXFLAGS="$CXXFLAGS -fopenmp -fPIC" -lgomp LTInversionArray.cpp
```
Note: [You will need to load the appropriate environment on Mahuika](#environment-on-mahuika)

If you cannot compile the C++ code, in line 100 of `ml.m' set cpp_code=0;

The function "MLapp.m" incorporates a Poincare asymptotic expansion for large time, and a Power series approximation for small time.

## Paramaters

### Model Paramaters

| Name | Type | Default | Desc. |
| --- | --- | --- | --- |
| `NSIM` | int[] | `[]` | Number of simulations to run |
| `alpha` | double | `exp(-2)` | set value of x(1) in simulation study |
| `beta` | double | `0.7` | set value of x(4) in simulation study |
| `gamma` | double | `1.5` | set value of x(3) in simulation study|
| `theta` | int | `1` | Parameter for magnitude distribution, do not change in simulation study. | 
| `lambda` | int | `3` | Parameter for magnitude distribution, do not change in simulation study. |
| `mu` | double | `1` | x(2), i.e 2nd parameter to be estimated via maximisation, corresponds to background rate (\lambda_0 in the manuscript) in simulation study. |
| `c` | double | `1` | set value of x(5) in simulation study|
| `M0` | double | `2.5` | Cutoff magnitude in data set |
| `MMAX` | double | `9.5` |Maximum permissible magnitude |
| `Nrand` | int | `30` | Number of initial guess for minimisation |
| `RepNum` | int | `500` | Number of repeats of simulating data and estimating it |
| `Batchsize`| int | `10` | Size of parallel batches in simulation study |
| `A`| int | `1` | First event to be considered when fitting actual data. |
| `B`| int | `1` | Last event to be considered when fitting actual data. |
| `PAR`| int | `0` | Swicth for minimising log-likelihood in parallel or not. Set to 0 if not using a parallel computer. |

### Execution parameters
| Name | Type | Default | Desc. |
| --- | --- | --- | --- |
| `parallel` | struct | | Struct containing info on how to run this job in paralell  |
|  | `outer` | `maxNumCompThreads > 2` | Type of parallel to use on outer loop (NSIMs)|
|  |  |  | `0` = Serial. |
|  |  |  | `1` = Parpool. |
|  | `inner` | `0` | Type of parallel to use on inner loop (consRep) |
|  |  |  | `0` = Serial. |
|  |  |  | `1` = Parpool. |
|  |  |  | `2` = Slurm. |
|  | `batchSize` | `1` | Number of times to run `consRep` on each worker. Output of `cons` should be `[4 x batchsize, 1]`  |
|  | `slurmParameters` | `{}` | Keyword + arguments to pass onto the slurm header. e.g. `{'--time', '01:00:00', '--mem', '1G'}`|
| `saveOutput` | bool / string(path) | `False` | 'False' for no output, or path to output. |

# FHP Simulation
To run the simulation load the file "input_files/simIN.m" then run "JRSSCSim4GITHUB/ExpConsCheck.m". It is recommended to use a small value of RepNum, Nrand, and NSIM if on a desktop computer. You will also need the parallel inner and outer structure above set to 0. 

# FHP data estimation
"DataEstimationLogit.m" and "DataEstimation0Lamb0Hessian.m" were used to find parameter and standard error estimates of the FHP and of the restricted FHP.
Ensure you have the files "JRSSCSim4GITHUB/LTInversionArray.cpp", "JRSSCSim4GITHUB/LTInversionArray.mexa64.m", "JRSSCSim4GITHUB/ml.m" and "JRSSCSim4GITHUB/MLapp.m" in the working directory of MATLAB when estimating data.

"CFHPEstimation.m" estimates MLE for the version of the fractional Hawkes process presented in Chen et al. (2021).


# ETAS estimation
This is done using "FitTimeETAS.R" and follows closely what was done by D. Harte in his cited guide. To run it enter in the information in the input_files for "JT", "LM" and "HM" verbatim. For "KK" follow the instructions in the script to uncomment certain sections of code.


# Predictive Capability
Use "PredCapETAS.m" to compute $bar p$ as in Section 3.4 and "PredictionCapability.m" to compute bar $p$ using "input_files/PredINLM.m". Then use "Compute information gains" of "miscScriptJRSSCFraqEq.m" to compute IGPE. To compute other IGPe use MLE as presented in the paper with input files for corresponding data set.


# Misc
Intensity of the FHP was computed using "Compute FHP intensity vector" section in the "miscScriptJRSSCFraqEq.m" script.
Transformed time residuals for the FHP was computed in the same script in the "SFHP and restricted SFHP residuals" section.
Transformed time residuals and intensity for the ETAS model is done using the "residuals" and "etas_gif" both in the PtProcess CRAN package by Harte (2010).

