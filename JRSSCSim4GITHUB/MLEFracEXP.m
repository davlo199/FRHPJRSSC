function NegLogLik=MLEFracEXP(x,D,MAG,M0,Tprime,M) 
Events=D;

m=exp(x(1))/(1+exp(x(1)));

lamb0=exp(x(2));

gam=exp(x(3)); 

beta=exp(x(4))/(1+exp(x(4)));

c=exp(x(5));

VECEXP=exp(gam*(MAG-M0));

% Int=zeros(size(Events));
% 
% cTime=(c*(Events(end)-Events)).^beta;
% 
% Intdx1=cTime>=Tprime; %Events to use Poin Asym App.
% 
% Int(Intdx1)=AppML(M,beta,-(c*(Events(end)-Events(Intdx1))).^beta,1);

Int=MLapp(Events(end)-Events,c,beta,1); %E_{\beta,1} used in the calculation of the integral

MAT=zeros(length(Events));

MAT3=zeros(length(Events));

for j=2:length(Events)  %Creating upper diagonal matrix for event times where entry (i,j) is t_j-t_i if this quantity is non-negative.

    MAT(1:j-1,j)=Events(j)-Events(1:j-1);
    
end

IDX2=MAT~=0;

% ZMat=(c.*MAT(IDX2)).^beta;

ZMat2=MAT(IDX2); 

% SumDx1=ZMat>=Tprime;

% TVec1=zeros(size(ZMat2));

% TVec1(SumDx1)=(c^beta).*(AppML(M,beta,-ZMat(SumDx1).',beta).*(ZMat2(SumDx1).^(beta-1)).').';

TVec1=(c^beta).*(MLapp(ZMat2.',c,beta,beta).*(ZMat2.^(beta-1)).').';

% Determining when to use power series approximation or the `ml.m' function
% if beta>=0.5
% 
%     Intdx2=cTime<=4;
% 
%     SumDx2=ZMat<=4;
% 
%     Int(Intdx2)=PSML(120,beta,-(c*(Events(end)-Events(Intdx2))).^beta,1);
% 
%     Intdx3=~(Intdx1|Intdx2);
% 
%     Int(Intdx3)=ml(-(c*(Events(end)-Events(Intdx3))).^beta,beta);
% 
%     TVec1(SumDx2)=(c^beta).*(PSML(120,beta,-ZMat(SumDx2).',beta).*(ZMat2(SumDx2).^(beta-1)).').';
% 
%     SumDx3=~(SumDx1|SumDx2);
% 
%     TVec1(SumDx3)=(c^beta).*ml(-ZMat(SumDx3),beta,beta).*(ZMat2(SumDx3).^(beta-1)).';
% 
% elseif beta<0.5 && beta>=0.25
% 
%     Intdx2=cTime<=2;
% 
%     SumDx2=ZMat<=2;
% 
%     Int(Intdx2)=PSML(350,beta,-(c*(Events(end)-Events(Intdx2))).^beta,1);
% 
%     Intdx3=~(Intdx1|Intdx2);
% 
%     Int(Intdx3)=ml(-(c*(Events(end)-Events(Intdx3))).^beta,beta);
% 
%     TVec1(SumDx2)=(c^beta).*(PSML(350,beta,-ZMat(SumDx2).',beta).*(ZMat2(SumDx2).^(beta-1)).').';
% 
%     SumDx3=~(SumDx1|SumDx2);
% 
%     TVec1(SumDx3)=(c^beta).*ml(-ZMat(SumDx3),beta,beta).*(ZMat2(SumDx3).^(beta-1)).';
% 
% else
% 
%     Intdx2=cTime<=1;
% 
%     SumDx2=ZMat<=1;
% 
%     Int(Intdx2)=PSML(100,beta,-(c*(Events(end)-Events(Intdx2))).^beta,1);
% 
%     Intdx3=~(Intdx1|Intdx2);
% 
%     Int(Intdx3)=ml(-(c*(Events(end)-Events(Intdx3))).^beta,beta);
% 
%     TVec1(SumDx2)=(c^beta).*(PSML(100,beta,-ZMat(SumDx2).',beta).*(ZMat2(SumDx2).^(beta-1)).').';
% 
%     SumDx3=~(SumDx1|SumDx2);
% 
%     TVec1(SumDx3)=(c^beta).*ml(-ZMat(SumDx3),beta,beta).*(ZMat2(SumDx3).^(beta-1)).';
% 
% end
% 
% 
% 
% %Computed after the cases are done
MAT3(IDX2)=TVec1.';

SUM=sum(log(lamb0+m.*(VECEXP*MAT3))); %Sum term of the log-likelihood
COMP=lamb0*D(end)-m.*(VECEXP*((Int-1).')); %Integral term of the log-likelihood
NegLogLik=-(SUM-COMP); 
end

function OUT=AppML(M,beta1,Z,beta2)
% OUT=0;
% for II=1:M
% OUT=OUT-(Z.^(-II))./gamma(beta-II.*beta);
% end
% Below is vectorised version of above and is the Poincare asymptotic
% series of the Mittag Leffler function
OUT=-(gamma(beta2-beta1.*(1:M)).^(-1))*(Z.^(-(1:M).'));

end
function OUT2=PSML(M,beta1,Z,beta2)
% OUT=0;
% for k=0:M
% OUT=OUT+(Z.^k)./gamma(beta2+k.*beta1);
% end
%Below is vectorised version of above and is the power series approximation
% to the Mittag-Leffler function. Note that since Z<0 log(Z)=|Z|+i\pi hence
% the use of the real function to remove any small imaginary part due to
% finite approximations.
OUT2=real(sum(exp((0:M).*log(Z).'-gammaln(beta2+beta1.*(0:M))),2)).';
OUT2(abs(Z)<=1e-15)=1./gamma(beta2);
end
