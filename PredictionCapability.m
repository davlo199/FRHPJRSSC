function PredictionCapability(input_file)
run(input_file)
%%
%Reading in datafile
Data=strcat(dataname,'.csv');
Points=readtable(Data);
Time=Points(:,1);
Mag=Points(:,5);
MAG=table2array(Mag);
MAG=MAG(1:end);
TimeCell=table2array(Time);
TimeNum=cell2mat(TimeCell);
datestr=TimeNum(:,1:10);
timestr=TimeNum(:,12:(end-1));
for i=1:height(Time)
datestr(i,1:11)=pad(datestr(i,:),11);
end
dtstr=[datestr,timestr];
numdatetime=datenum(dtstr); %Absolute time of events in serial date which is number of days
Events=numdatetime; %Removes event at time 0, relative number of days since first event
idx=MAG>=M0;
Events=Events(idx).';
MAG=MAG(idx).';
Events=Events(A:B);
MAG=MAG(A:B);
Events=Events-Events(1);

%Setting up simulation
%%
PredInt=[1e-6:0.25:Events(end), Events(end)]; %End points of time bins, each is 6 hours 
%PredInt=linspace(1e-6,Events(end),2000) for 2000 evenly sized intervals

NInt=length(PredInt)-1; %Number of intervals for prediction

batchsize=10;
NRep=1000; %Number of simulations for each interal
if PAR==0
PSFHP=zeros(batchsize,ceil(NInt/batchsize));
    for K=1:ceil(NInt/batchsize)
    PSFHP(:,K)=PredictionLoop(alp, bet, gam, m, M0, c, Events, MAG, NRep, K, PredInt,batchsize);
    if floor(K/10)==K/10
        K
    end
    end
else
    functionhandle=@(K)PredictionLoop(alp, bet, gam, m, M0, c, Events, MAG, NRep, K, PredInt,batchsize);
    PSFP=slurmParallel(functionhandle,1:ceil(NInt/batchsize),parallel.slurmParameters);
end
basefilename=strcat(dataname,'2sfhpPrediction','A',sprintf('%d',A),'B',sprintf('%d',B));

save(basefilename)

end

function OUT=PredictionLoop(alpha, beta, gamma, mu, M0, c, Events, MAG, NRep, kInt, PredInt,batchsize)
OUT=zeros(1,batchsize);
if kInt~=ceil((length(PredInt)-1)/batchsize)
for i=1:batchsize
    rng((kInt-1)*batchsize+i) %For reproducibility
    T0=PredInt((kInt-1)*batchsize+i);
    TF=PredInt((kInt-1)*batchsize+i+1);
    %History fot his interval
    TmpEvent=Events(Events<=T0);
    TmpMAG=MAG(Events<=T0);
    Prop=zeros(1,NRep);
        for L=1:NRep
        Prop(L)=SimPrediction(alpha, beta, gamma, mu, M0, c,TmpEvent,TmpMAG,T0,TF);
        end
    OUT(i)=sum(Prop)/NRep;
end
else
    for i=1:(1-kInt)*batchsize+length(PredInt)-1
        rng((kInt-1)*batchsize+i) %For reproducibility
    T0=PredInt((kInt-1)*batchsize+i);
    TF=PredInt((kInt-1)*batchsize+i+1);
    %History fot his interval
    TmpEvent=Events(Events<=T0);
    TmpMAG=MAG(Events<=T0);
    Prop=zeros(1,NRep);
        for L=1:NRep
        Prop(L)=SimPrediction(alpha, beta, gamma, mu, M0, c,TmpEvent,TmpMAG,T0,TF);
        end
    OUT(i)=sum(Prop)/NRep;
    end
end
end

function N = SimPrediction(alpha, beta, gamma, mu, M0, c,Events,MAG,T0,TF)
%Events and MAG constitute history, (T0,TF] is time interval we want to
%predict in
%Rest are parameters in the simulation
%Only need to determine whether there is one point in (T0,TF] since
%considering unmarked case.
    t = T0;
    epsilon = 1e-10;
    SimPoints = Events;
    N=0; %Number of simulated points
    % Simulate Points
   
    while t<=TF
        M = mu + alpha * ((exp(gamma * (MAG - M0))) ...
            * (MLapp(t+epsilon-SimPoints,c,beta,beta).*(t+epsilon-SimPoints).^(beta-1)).');
        E = exprnd(1 / M, 1, 1);
        t = t + E;
        if t>TF
            break
        end
        U = rand;

        if (U < (mu + alpha * ((exp(gamma * (MAG - M0))) * (MLapp(t-SimPoints,c,beta,beta) ...
                .* (t - SimPoints).^(beta - 1)).')) / M)
            N=N+1;
            
            break %Stops after we have simulate one event, and has N=1
        end

    end
    
end
