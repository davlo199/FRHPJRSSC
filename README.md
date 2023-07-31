# FRHPJRSSC
All of the scripts (and data) used in "A Fractional Model for Earthquakes" submitted to JRSSC Series C

The fractional Hawkes process is now referred to as FHP.
Notes:
These are the data files (JT.csv, LM.csv HM.csv and KaikoraQuakes.csv) that were used. To replicate results in the manuscript, one needs to remove events in LM.csv that have magnitude <2.75 and then only use events 1 to 1348 of those left and only use events 1 to 1023 in HM.csv (parameters for this are in the input files).

Scripts within the "Consistency" file are what was used for the simulation study, along with the "ml.m" function and its C++ adaption as mentioned. 
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



## Paramaters

### Model Paramaters

| Name | Type | Default | Desc. |
| --- | --- | --- | --- |
| `NSIM` | int[] | `[]` | Number of simulations to run |
| `alpha` | double | `exp(-2)` | set value of x(1) |
| `beta` | double | `0.7` | set value of x(4)  |
| `gamma` | double | `1.5` | set value of x(3) |
| `theta` | int | `1` | Parameter for magnitude distribution, do not change. | 
| `lambda` | int | `3` | Parameter for magnitude distribution, do not change. |
| `mu` | double | `1` | x(2), i.e 2nd parameter to be estimated via maximisation, corresponds to background rate (\lambda_0 in the manuscript). |
| `c` | double | `1` | set value of x(5) |
| `M0` | double | `2.5` | Cutoff magnitude in data set |
| `MMAX` | double | `9.5` |Maximum permissible magnitude |
| `Nrand` | int | `30` | Number of initial guess for minimisation |
| `RepNum` | int | `500` | Number of repeats of simulating data and estimating it |
| `Batchsize`| int | `10` | Size of parallel batches |


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



# FHP data estimation
"DataEstimationLogit.m" and "DataEstimation0Lamb0Hessian.m" were used to find parameter and standard error estimates of the FHP and of the restricted FHP.

# ETAS estimation
This is done using "FitTimeETAS.R" and follows closely what was done by D. Harte in his cited guide.
# Misc
Intensity of the FHP was computed using "Compute FHP intensity vector" section in the "Creatboxplots.m" script.
Transformed time for the FHP was computed in the same script in the "Compute transformed time residuals" section.
Transformed time for the ETAS model was computed in "ResidualonSamePlotWCl.R".
Intensity of the ETAS model is computed in "IntPlot.R".
