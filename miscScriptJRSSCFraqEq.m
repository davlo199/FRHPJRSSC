%% Compute information gains
%Using the output from "PredCapETAS.m" and "PredicitionCapability.m" below
%computes IGPe
Tmp=size(PSFHP);
PTest=reshape(PSFHP,1,Tmp(1)*Tmp(2)); %For output of PredictionCapability.m
% PTest=P; %For output of PredCapETAS.m
X=zeros(1,NInt);
for I=1:NInt
T0=PredInt(I);
TF=PredInt(I+1);
if any(Events>T0&(Events<=TF))
    X(I)=1;
end
end
Xlogic=logical(X);
delI=diff(PredInt);

PPois=expcdf(delI,Events(end)/length(Events)); %MATLAB uses mean as parameter
PTest=PTest(1:length(PPois)); %Removes zeros on end from batchsizing
IGPe=(1/length(Events))*(sum(log(PTest(Xlogic)./PPois(Xlogic)))+sum(log((1-PTest(~Xlogic))./(1-PPois(~Xlogic)))));

%% SFHP and restricted SFHP residuals
Comp=[];
if length(FinalParameters)==4
m=FinalParameters(1);
gam=FinalParameters(2);
beta=FinalParameters(3);
c=FinalParameters(4);

N=length(Events);
for i=1:N
VECEXP=exp(gam*(MAG(1:i)-M0));
VEC2=MLapp(Events(i)-Events(1:i),c,beta,1);
Comp(i)=-m.*(VECEXP*((VEC2-1).'));
end

else
m=FinalParameters(1);
lamb0=FinalParameters(2);
gam=FinalParameters(3);
beta=FinalParameters(4);
c=FinalParameters(5);

N=length(Events);
for i=1:N
VECEXP=exp(gam*(MAG(1:i)-M0));

VEC2=MLapp(Events(i)-Events(1:i),c,beta,1);
Comp(i)=-m.*(VECEXP*((VEC2-1).'));
end
Comp=lamb0.*Events+Comp;

end

%% Create FHP intensity vector
%Set up time vector 
eps=1e-3;
A=1;
B=1348;
N=500;
Events2=Events(A:B);
dt=(Events2(end)-Events2(1))/(N-1);
MAG2=MAG(A:B);
t1=Events2+eps;
t2=Events2(2:end)-eps;
t3=Events2(1):dt:Events2(end);
tvec=sort([t1,t2,t3,Events2]);
%tvec=Events2;
%Create intensity vector
m=0.402;
lamb0=0.414;
gam=1.273;
beta=0.665;
c=0.849;


VEC=zeros((B-A)+1,length(tvec));
VecTime=VEC;

VECEXP=exp(gam*(MAG2-M0));

for j=1:length(tvec) 
    JDX=Events2<tvec(j);
    VecTime(JDX,j)=tvec(j)-Events2(JDX);  
end

TDX=VecTime==0;
VEC=((c^beta).*reshape(ml(-(c*VecTime).^beta,beta,beta,1),size(VecTime)).*(VecTime.^(beta-1)));
VEC(TDX)=0;

Intensity=lamb0+m.*(VECEXP*VEC);

