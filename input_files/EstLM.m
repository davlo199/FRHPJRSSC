M0=2.75;
batchsize=1;
dataname='LM';
CSV=1;
Nrand=50;
PAR=0;
M=13;
Tprime=10;
maxFval=3000;
A=1;
B=1348;
oldMLE=[];
parallel.slurmParameters={'--mem-per-cpu', '200m', '--time', '00:40:00','--cpus-per-task','12'};