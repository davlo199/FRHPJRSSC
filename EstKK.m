M0=2.75;
dataname='KaikoraQuakes';
Nrand=30;
CSV=1;
PAR=0;
A=1;
B=2398;
SEED=1;
Nparam=7;
oldMLE=[];
parallel.slurmParameters={'--mem-per-cpu', '500m', '--time', '5:00:00','--cpus-per-task','12'};
