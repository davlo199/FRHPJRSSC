#!/bin/bash

#SBATCH --mem               2G
#SBATCH --cpus-per-task     2           # This should be equal to 'length(NSIM) * 2'
#SBATCH --time              00:40:00    # Should be about 120% of the time you estimate needed
#SBATCH --output            PredKKSFHP.txt  # What would usually be printed to terminal will go here.

# Above is the 'Slurm header', it defines the metaparameters relating to your job.

# Any line starting with a '#' is a comment and will not be read executed. Same as matlab '%'

module load MATLAB/2020b
# Makes matlab available to use.

# Best practise to choose as close to the version on your local machine.
# e.g.
# module load MATLAB/2021b

matlab -nodisplay -r "PredictionCapability input_files/PredInKK"

# To run this script.
# 'cd <path/to/code>'
# 'sbatch submit.sl'
# Then type `sacct` to see status of job.
