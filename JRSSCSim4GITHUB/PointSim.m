%clear
%filename='Ogata06MLESIM';
%%%% To simulate points from Ogata2006MLE
NSIM=1348;
Tfinal=180;
lambda=7.908;
theta=16115.796;
%a=1;
M0=2.75; %Cutoff magnitude in data set
MMAX=10; %Maximum permissible magnitude 
M1=7.3; %Magnitude of first event
N=NSIM;
mu=theta;
K=1/lambda;
sigma=M0/lambda;
A=sigma/K;
magexp=exprnd(mu,1,N)+M0; %%W-a random variable
magPar=gprnd(K,sigma,A,1,N);
mag=min(MMAX,min(magexp,magPar));


bHat=log10(exp(1))./(mean(mag)-M0);

%%
% for II=1:length(NSIM)



n=0;
alpha=0.097; 
mu=0.4763;
epsilon=1e-10;
beta=0.5282;
c=0.8223;
gamma=2.1013;



SimPoints=0;  %%Sets initial event large to check explosion
t=SimPoints;
MAG2=M1;
while t<=Tfinal
t;
M=mu+alpha*((exp(gamma*(MAG2-M0)))*c^beta*(ml(-(c.*(t+epsilon-SimPoints)).^beta,beta,beta,1).*(t+epsilon-SimPoints).^(beta-1)).');
E=exprnd(1/M,1,1);
t=t+E;
U=rand;
if (U<(mu+c^beta.*alpha*((exp(gamma*(MAG2-M0)))*(ml(-(c.*(t-SimPoints)).^beta,beta,beta,1)...
        .*(t-SimPoints).^(beta-1)).'))/M)
n=n+1;
SimPoints = [SimPoints, t];
magexp=exprnd(mu,1,1)+M0; %%W-a random variable
magPar=gprnd(K,sigma,A,1,1);
mag=min(MMAX,min(magexp,magPar));
MAG2=[MAG2,mag];
end
end
%%
IDX=Events<=180;
subplot(2,1,1)
plot(SimPoints(1:end-1),MAG2(1:end-1),'x')
title("Simulated")
subplot(2,1,2)
plot(Events(IDX),MAG(IDX),'x')
title("Actual")

%%
N=[100,500,750,1023,2500,3750,5000];
lambda=2.7179;
theta=0.78053;
M0=2.5; %Cutoff magnitude in data set
MMAX=10; %Maximum permissible magnitude 
M1=7.3; %Magnitude of first event
mu=theta;
K=1/lambda;
sigma=M0/lambda;
A=sigma/K;
Mean=[];
Std=[];

for i=1:length(N)
    
    magexp=exprnd(mu,1e3,N(i))+M0; %%W-a random variable
    magPar=gprnd(K,sigma,A,1e3,N(i));
    mag=min(MMAX,min(magexp,magPar));
    Mean(i,:)=max(mag.');
                
end
%%
boxplot(Mean.',N)
xlabel("Number of simulated events")
ylabel("Maximum magntiude of trial")


