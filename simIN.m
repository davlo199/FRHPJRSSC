% Number of simulations to run
NSIM = 30;
RepNum=500;
Nrand=5;
saveOutput = 'JRSSCSim30.mat'; %String. Path to '.mat'
parallel.batchSize=5;
parallel.inner=0;
parallel.outer=0;
parallel.slurmParameters={'--mem-per-cpu', '1250m', '--time', '16:00:00','--cpus-per-task','12'};
