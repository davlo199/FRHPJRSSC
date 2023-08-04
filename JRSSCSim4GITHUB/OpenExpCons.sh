#!/bin/bash

module load MATLAB

matlab -nosplash -nodisplay -r "ExpConsCheck input_files/testInputSerial; quit"
