M0=2.5;
batchsize=1;
dataname='JT';
Nrand=30;
PAR=0;
A=1;
CSV=1;
MaxFval=3000;
B=761;

Nparam=5;

oldMLE=[];

parallel.slurmParameters={'--mem-per-cpu', '300m', '--time', '00:45:00','--cpus-per-task','12'};
