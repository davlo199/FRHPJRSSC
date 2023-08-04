function DataEstimation0Lamb0Hessian(input_file)
run(input_file)
%%
Data=strcat(dataname,'.csv');
Points=readtable(Data);
        Time=Points(:,1);
        Mag=Points(:,5);
        % Comment lines 6/7 and uncomment 9/10 to run KK
%         Time=Points(:,3);
%         Mag=Points(:,7);
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
        Events=numdatetime(1:end)-numdatetime(1);
        
idx=MAG>=M0;
Events=Events(idx).';
MAG=MAG(idx).';
Events=Events(A:B);
MAG=MAG(A:B);

%%
start=clock;


D=Events; 





switch PAR
    case 1
    MLLik=zeros(Nrand,5);
    parpool(nesiSafeParpool(Nrand));
    parfor K=1:Nrand
    MLLik(K,:)=MinimisationLoop(D,MAG,M0);
    end
    case 2        
        functionhandle=@(Nrand)MinimisationLoop(D,MAG,M0,Nrand);
        MLLikTemp=slurmParallel(functionhandle,1:Nrand,parallel.slurmParameters);
        MLLik=MLLikTemp.';
    otherwise
    MLLik=zeros(Nrand,9);
    for K=1:Nrand
    [MLLik(K,:)]=MinimisationLoop(D,MAG,M0,K);
    end
end


%% %%Check MLE
Likelihood=MLLik(:,end);
% Find best mle for given beta
MINlik=min(Likelihood);
idx = Likelihood==MINlik;
if sum(idx)>1
AA=MLLik(idx,1:4);
BestMLE=[AA(1,:) MINlik]; %Vector of MLE and the likelihood
else
BestMLE=[MLLik(idx,1:4) MINlik];
end

BestMinlik=min(BestMLE(:,end));
FinalParameters=BestMLE((BestMLE(:,end)==BestMinlik),1:4);
FinalParameters(1)=exp(FinalParameters(2))/(1+exp(FinalParameters(2)));
FinalParameters(2)=exp(FinalParameters(2));
FinalParameters(3)=exp(FinalParameters(3))/(1+exp(FinalParameters(3)));
FinalParameters(4)=exp(FinalParameters(4));

%%
basefilename=strcat(dataname,'FHP0Lamb0','A',sprintf('%d',A),'B',sprintf('%d',B));

endtime=clock;
runtime=endtime-start; 
save(basefilename)


end

%% Auxilliary functions


function MLLik=MinimisationLoop(D,MAG,M0,nInitial)
    
    options=optimoptions('fminunc','MaxFunctionEvaluations',3000);
    rng(2*nInitial)
    x0=rand(1,4); 
    [MLE,Llhood,~,~,~,Hessian]=fminunc(@(x)MLEFracEXP(x,D,MAG,M0),x0,options); 
    SE=sqrt(diag(inv(Hessian))).';
    MLLik=[MLE,SE,Llhood];
    
end

function NegLogLik=MLEFracEXP(x,D,MAG,M0) %Using predetermined matrix
%D events, MAG is magnitude of the events, and M0 minimum magnitude
%%Log likelihood expression
Events=D;
m=exp(x(1))/(1+exp(x(1)));
gam=exp(x(2)); 
beta=exp(x(3))/(1+exp(x(3)));
c=exp(x(4));

VEC=zeros(length(Events));



MAT=zeros(length(Events));

MAT3=zeros(length(Events));

for j=2:length(Events)  

    MAT(1:j-1,j)=Events(j)-Events(1:j-1);
    
end

IDX2=MAT~=0;

ZMat2=MAT(IDX2); 

TVec1=(c^beta).*(MLapp(ZMat2.',c,beta,beta).*(ZMat2.^(beta-1)).').';


MAT3(IDX2)=TVec1.';

VEC2=MLapp(Events(end)-Events,c,beta,1);

%integral over time domain of intensity is unchanged since the lower limit of integration is still the first event 
%%Finding the first sum


VECEXP=exp(gam*(MAG-M0));
A=VECEXP*MAT3;
SUM=sum(log(m.*A(2:end)));

%Calculating intensity at points after first since, left limit of intensity
%at t_1 is 0 and so log is infinite

COMP=-m.*(VECEXP*((VEC2-1).'));


NegLogLik=-(SUM-COMP); %Negative log likelihood
end