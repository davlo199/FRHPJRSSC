function PredCapETAS(input_file)
run(input_file)
%%
%Reading in datafile
Data=strcat(dataname,'.csv');
Points=readtable(Data);
% Time=Points(:,1); %For USGS data i.e. JT,LM,HM
% Mag=Points(:,5);
Time=Points(:,3); %For GEONET data i.e. KK
Mag=Points(:,7);
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
Events=Events(AA:B);
MAG=MAG(AA:B);
Events=Events-Events(1);

%Setting up simulation
%%
PredInt=[1e-6:0.25:Events(end), Events(end)]; %End points of time bins, each is 6 hours 


NInt=length(PredInt)-1; %Number of intervals fr prediction


NRep=1000; %Number of simulations for each interal

P=zeros(1,NInt);
    for K=1:NInt
    P(K)=PredictionLoop(mu0, A, del, p, M0, cE, Events, MAG, NRep, K, PredInt);
    if floor(K/50)==K/50
        K
    end
    end

basefilename=strcat(dataname,'ETASPrediction','A',sprintf('%d',AA),'B',sprintf('%d',B));

save(basefilename)



end



function OUT=PredictionLoop(mu0, A, del, p, M0, cE, Events, MAG, NRep, kInt, PredInt)

    rng(kInt) %For reproducibility
    T0=PredInt(kInt);
    TF=PredInt(kInt+1);
    %History for this interval
    TmpEvent=Events(Events<=T0);
    TmpMAG=MAG(Events<=T0);
    Prop=zeros(1,NRep);
        for L=1:NRep
        Prop(L)=SimPrediction(mu0, A, del, p, M0, cE,TmpEvent,TmpMAG,T0,TF);
        end
    OUT=sum(Prop)/NRep;

end


function N = SimPrediction(mu0, A, del, p, M0, cE,Events,MAG,T0,TF)
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
        M = mu0 + A * exp(del * (MAG - M0))*((1+(t+epsilon-SimPoints)/cE).^(-p)).' ;
        E = exprnd(1 / M, 1, 1);
        t = t + E;
        if t>TF
            break
        end
        U = rand;

        if (U < (mu0 + A * exp(del * (MAG - M0)) * ((1+(t-SimPoints)/cE).^(-p)).' / M))
            N=N+1;
            
            break %Stops after we have simulate one event, and has N=1
        end

    end
    
end
