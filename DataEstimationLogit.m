function DataEstimationLogit(input_file) %parameter estimation for SFHP
run(input_file)
%% Reading in data
if CSV==1
Data=strcat(dataname,'.csv');
Points=readtable(Data);
Time=Points(:,1);
Mag=Points(:,5);
%Comment lines 7/8 for KK and uncomment 10/11 for KK
%Time=Points(:,3);
%Mag=Points(:,7);

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
%Parallelisation use PAR=0 for desktop PC. PAR=2 for slurm batch parallelisation
switch PAR
    case 1
    MLLik=zeros(Nrand,11);
    parpool(nesiSafeParpool(Nrand));
    parfor K=1:Nrand
    MLLik(K,:)=Minimisation(D,MAG,M0);
    end
    case 2        
        functionhandle=@(Nrand)MinimisationLoop(D,MAG,M0,Nrand,oldMLE);
        MLLikTemp=slurmParallel(functionhandle,1:Nrand,parallel.slurmParameters);
        MLLik=MLLikTemp.';
    otherwise
    MLLik=zeros(Nrand,11);
    for K=1:Nrand
    MLLik(K,:)=MinimisationLoop(D,MAG,M0,K,oldMLE);
    end
end



%% %%Check MLE
Likelihood=MLLik(:,end);
%% Find best mle from all initial values
MINlik=min(Likelihood);
idx = Likelihood==MINlik;
if sum(idx)>1
AA=MLLik(idx,1:10); %Selecting only one set in case minimiser finds the same minimum multiple times
BestMLE=[AA(1,:) MINlik]; %Vector of MLE and the likelihood
else
BestMLE=[MLLik(idx,1:10) MINlik];
end

BestMinlik=min(BestMLE(:,end));
FinalParameters=BestMLE((BestMLE(:,end)==BestMinlik),1:10);
%Back transforming the parameters
FinalParameters(2)=exp(FinalParameters(2));
FinalParameters(3)=exp(FinalParameters(3));
FinalParameters(5)=exp(FinalParameters(5));
FinalParameters(1)=exp(FinalParameters(1))/(1+exp(FinalParameters(1)));
FinalParameters(4)=exp(FinalParameters(4))/(1+exp(FinalParameters(4)));


basefilename=strcat(dataname,'FHP','A',sprintf('%d',A),'B',sprintf('%d',B));

endtime=clock;
runtime=endtime-start; %Returns parameters alpha lambda gamma beta negative loglikelihood
save(basefilename)


end


%%
function MLLik=MinimisationLoop(D,MAG,M0,nInitial,oldMLE)
    options=optimoptions('fminunc','MaxFunctionEvaluations',3000);
    
    rng(2*nInitial)
    if ~isempty(oldMLE)
        x0=oldMLE; %oldMLE is for choosing the initial values if one has found the approximate minima and now just needs to compute it to a lower level of accuracy
    else
    x0=rand(1,5); %Random initial values, could change to randn(1,5) for normally distributed
    end
    [MLE,Llhood,~,~,~,Hessian]=fminunc(@(x)MLEFracEXP(x,D,MAG,M0),x0,options); 
     SE=sqrt(diag(inv(Hessian))).';
    MLLik=[MLE,SE,Llhood];
    
end

function NegLogLik=MLEFracEXP(x,D,MAG,M0) %Using predetermined matrix
%D events, MAG is magnitude of the events, and M0 minimum magnitude
%%Log likelihood expression
Events=D;
%Transforming parameters for unconstrained minimisation
m=exp(x(1))/(1+exp(x(1)));
lamb0=exp(x(2));
gam=exp(x(3)); 
beta=exp(x(4))/(1+exp(x(4)));
c=exp(x(5));


MAT=zeros(length(Events));
MAT2=MAT;

for j=2:length(Events)  %Setting up event matrix where entry i,j is t_j-t_i if this is non-negative (for use in sum of log-intensities)
    MAT(1:j-1,j)=Events(j)-Events(1:j-1);
    
end

IDX2=MAT~=0;

ZMat2=MAT(IDX2); 

TVec1=(c^beta).*(MLapp(ZMat2.',c,beta,beta).*(ZMat2.^(beta-1)).').'; %Calculating sum of log-intensities at event times

% TempVEC=(c^beta).*ml(-(c.*MAT(IDX2)).^beta,beta,beta,1).*(MAT(IDX2).^(beta-1)).';

MAT2(IDX2)=TVec1.';

%VEC2=ml(-(c*(Events(end)-Events)).^beta,beta);
VEC2=MLapp(Events(end)-Events,c,beta,1); 


VECEXP=exp(gam*(MAG-M0));

SUM=sum(log(lamb0+m.*(VECEXP*MAT2))); %The sum of the log intensities at event times

COMP=lamb0*D(end)-m.*(VECEXP*((VEC2-1).')); %Vec2 is just sum of ML function at -(c(T-t_i))^beta which is whats needed for the integral

NegLogLik=-(SUM-COMP); %Negative log likelihood
end
