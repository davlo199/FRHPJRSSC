function CFHPEstimation(input_file)
%MLE for Chen et al. FHP model.
run(input_file)
%%
if CSV==1
Data=strcat(dataname,'.csv');
Points=readtable(Data);
Time=Points(:,1);
Mag=Points(:,5);
%Comment 8/9 and uncomment 11/12 for KK
% Time=Points(:,3);
% Mag=Points(:,7);
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
else
    Data=strcat(dataname,'.txt');
    Points=readtable(Data);
    Time=Points(:,1);
    Mag=Points(:,2);
    MAG=table2array(Mag);
    MAG=MAG(1:end);
    TimeCell=table2array(Time);
%     TimeNum=cell2mat(TimeCell);
%     datestr=TimeNum(:,1:10);
%     timestr=TimeNum(:,12:(end-1));
%     for i=1:height(Time)
%         datestr(i,1:11)=pad(datestr(i,:),11);
%     end
%     dtstr=[datestr,timestr];
%     numdatetime=datenum(dtstr); %Absolute time of events in serial date which is number of days
    Events=TimeCell; %Removes event at time 0, relative number of days since first event
    idx=MAG>=M0;
    Events=Events(idx).';
    MAG=MAG(idx).';
    Events=Events(A:B);
    MAG=MAG(A:B);
    Events=Events-Events(1);
end

%%
start=clock;


D=Events; 

%Switch for parallelisation options on NeSI HPC. Use PAR=0 on a standard PC
switch PAR
    case 1
    MLLik=zeros(Nrand,4);
    parpool(nesiSafeParpool(Nrand));
    parfor K=1:Nrand
    MLLik(K,:)=Minimisation(D,MAG,M0);
    end
    case 2        
        functionhandle=@(Nrand)MinimisationLoop(D,MAG,M0,Nrand,oldMLE);
        MLLikTemp=slurmParallel(functionhandle,1:Nrand,parallel.slurmParameters);
        MLLik=MLLikTemp.';
    otherwise
    MLLik=zeros(Nrand,4);
    for K=1:Nrand
    MLLik(K,:)=MinimisationLoop(D,MAG,M0,K,oldMLE);
    end
end



%% %%Check MLE
Likelihood=MLLik(:,end);

MINlik=min(Likelihood);
idx = Likelihood==MINlik;
if sum(idx)>1
AA=MLLik(idx,1:3); %Selecting only one set in case minimiser finds the same minimum multiple times
BestMLE=[AA(1,:) MINlik]; %Vector of MLE and the likelihood
else
BestMLE=[MLLik(idx,1:3) MINlik];
end

BestMinlik=min(BestMLE(:,end));
FinalParameters=BestMLE((BestMLE(:,end)==BestMinlik),:);
%Back transforming the parameters
FinalParameters(2)=exp(FinalParameters(2));
FinalParameters(1)=exp(FinalParameters(1))/(1+exp(FinalParameters(1)));
FinalParameters(3)=exp(FinalParameters(3))/(1+exp(FinalParameters(3)));


basefilename=strcat(dataname,'ChenFHP','A',sprintf('%d',A),'B',sprintf('%d',B));

endtime=clock;
runtime=endtime-start; %Returns parameters alpha lambda gamma beta negative loglikelihood
save(basefilename)


end


%%
function MLLik=MinimisationLoop(D,MAG,M0,nInitial,oldMLE)
    options=optimoptions('fminunc','MaxFunctionEvaluations',3000);
    
    rng(2*nInitial)
    if ~isempty(oldMLE)
        x0=oldMLE;
    else
    x0=rand(1,3); 
    end
    [MLE,Llhood]=fminunc(@(x)MLEFracEXP(x,D,MAG,M0),x0,options); 
    MLLik=[MLE,Llhood];
    
end

function NegLogLik=MLEFracEXP(x,D,MAG,M0) %Using predetermined matrix
%D events, MAG is magnitude of the events, and M0 minimum magnitude
%%Log likelihood expression
Events=D;
%Transforming parameters for unconstrained minimisation

if isinf(exp(x(1)))
    m=1;
else
    m=exp(x(1))/(1+exp(x(1)));
end

lamb0=exp(x(2));

gam=0; %Restricted case of Chen et al. model
if isinf(exp(x(3)))
    beta=1;
else
beta=exp(x(3))/(1+exp(x(3)));
end
c=1; %Restricted case of Chen et al. model


MAT=zeros(length(Events));
MAT2=MAT;

for j=2:length(Events) 
    MAT(1:j-1,j)=Events(j)-Events(1:j-1);
    
end

IDX2=MAT~=0;

ZMat2=MAT(IDX2); 

TVec1=(c^beta).*(MLapp(ZMat2.',c,beta,beta).*(ZMat2.^(beta-1)).').';

% TempVEC=(c^beta).*ml(-(c.*MAT(IDX2)).^beta,beta,beta,1).*(MAT(IDX2).^(beta-1)).';

MAT2(IDX2)=TVec1.';

%VEC2=ml(-(c*(Events(end)-Events)).^beta,beta);
VEC2=MLapp(Events(end)-Events,c,beta,1);


VECEXP=exp(gam*(MAG-M0));

SUM=sum(log(lamb0+m.*(VECEXP*MAT2)));

COMP=lamb0*D(end)-m.*(VECEXP*((VEC2-1).')); %Vec2 is just sum of ML function at -(T-t_i)^beta

NegLogLik=-(SUM-COMP); %Negative log likelihood
end