function  [n,nchi,nchi3,Table3,Table5,RiskPremia,CIbRiskPremia,CIuRiskPremia,nuT,CIbnuT,CIunuT,TimeVarTest]=EstTimeVariantMod(K,p,q,n,T,ExcR,MC,I,F,Z,Zi,dataset)
% This function estimates time-varying risk premia. 
% OUTPUT: 
% n=cross-sectional dimension; 
% nchi=trimmed cross-sectional dimension used in the estimation approach;
% nchi3=trimmed cross-sectional dimension used to compute statistics;
% Table3 is a (Kp)x6 matrix that contains estimated annualized components
% and their CI of vec(F') (col 1,2,3) and nu (col 4,5,6);
% Table5 is a 2x2 matrix that contains statistics (col 1) based on Q_{e}
% (1st row) and on Q_{a} (2nd row) and their p-values (col 2);
% RiskPremia is a Kx(T-1) matrix of estimated annualized time-varying
% risk premia;
% CIbRiskPremia and CIuRiskPremia contain CI of risk premia;
% TimeVarTest contains stat and pv of time variation test.

% In the conditional model we loose one observation
ExcR=ExcR(2:T,:); I=I(2:T,:); Ti=sum(I)';

% Definition of dimension d,d1,d2 - Section 3.1
% \b_{i,t}=B_{i}*Z_{t-1}+C_{i}*Z_{i,t-1}; \lambda_{t}=\Lamdba*Z_{t-1}
d1=p*(p+1)/2+p*q;
d2=K*(p+q);
d=d1+d2;

%First pass (see Section 3.2)
%Estimation of betas 
[beta1,beta2,res,CN,V,X,~,~,~,~]=TSRegress_CN_V_cond(ExcR,F,Z,Zi,I,Ti,n,T,d1,d2,d);        

%Defition of tau
tau=(T-1)./Ti;      %nx1
% Trimming Approach
chi1=15;
chi2=(T-1)/60;
ind=find(CN<=chi1);
beta1=beta1(:,ind); beta2=beta2(:,ind); 
V=V(ind,:); X=X(:,:,ind); I=I(:,ind); res=res(:,ind); tau=tau(ind,:); Ti=Ti(ind,:); 
ind=find(tau<=chi2);
beta1=beta1(:,ind); beta2=beta2(:,ind); 
V=V(ind,:); X=X(:,:,ind); I=I(:,ind); res=res(:,ind); tau=tau(ind,:); Ti=Ti(ind,:); 
nchi=size(Ti,1);   %trimmed cross-sectional dimension

%Second Pass            (see eq. (13))
%Cross-sectional regression
[nu,vecbeta3]=CR_reg_cond(beta1,beta2,V,tau,K,p,q,nchi,d1,d2,d,[],0);
%Definition of \hat{F} 
Fhat=F(2:T,:)'*Z(1:T-1,:)*(inv(Z(1:T-1,:)'*Z(1:T-1,:)));        %Kxp

if dataset~=3;      %Portfolios: n small, T large         
    %Time-varying risk premia estimator 
    LAMBDAv=nu+vec(Fhat');           %Kpx1
    LAMBDA=(invVec(LAMBDAv,p,K))';   %Kxp;
    lambda=LAMBDA*(Z(1:T-1,:)');     %Kx(T-1)
    %Time-varying nu estimator
    NU=(invVec(nu,p,K)');   %Kxp;
    nuT=NU*(Z(1:T-1,:)');   %Kx(T-1)
    %Confidence Intervals of estimators vec(F'), nu and lambda (see Prop.5)
    [CIvecF,CInu,CIb,CIu,CIbnuT,CIunuT,Tij,SIGMAF,SIGMAnu]=CIlambda_nu_F_nfix(Ti,I,res,V,X,F,Z,[],vecbeta3,nu,nuT,lambda,Fhat,nchi,T,K,p,q,d1,d2,d,0);
    %Statistics: Section 3.5
    nchi3=nchi;        %no trimming for panel of portfolios
    [stat,pv,stat_a,pv_a]=GOS_statistics(nu,beta1,beta2,vecbeta3,V,X,I,res,Tij,T,d1,d2,d,nchi3,p,q,K);
else                %Individual stocks: n large, T large
    %Unbiased estimator for nu (see Prop. 4, 5)
    [nuB,~]=UnbiasedEst_Cond(nu,vecbeta3,V,K,p,q,nchi,T,tau,d1,d2,d,[],0);
    %Unbiased estimator of risk premia 
    LAMBDAv=nuB+vec(Fhat');          %Kpx1
    LAMBDA=(invVec(LAMBDAv,p,K))';   %Kxp;
    lambda=LAMBDA*(Z(1:T-1,:)');     %Kx(T-1)
    %Time-varying nu estimator
    NU=(invVec(nuB,p,K)');   %Kxp;
    nuT=NU*(Z(1:T-1,:)');   %Kx(T-1)
    %Confidence Intervals of estimators vec(F'), nu and lambda (see Prop.5)
    [CIvecF,CInu,CIb,CIu,CIbnuT,CIunuT,tau,SIGMAF,SIGMAnu]=CIlambda_nu_F_nTinfty(Ti,I,res,V,X,F,Z,[],vecbeta3,nuB,nuT,lambda,Fhat,nchi,T,K,p,q,d1,d2,d,chi2,0);    
    %Statistics: Section 3.5
    % Trimming (see Section 4.3 and results on Monte Carlo simulations in supplementary materials)
    chi3=(T-1)/240;
    ind=find(diag(tau)<=chi3);
    beta1=beta1(:,ind); beta2=beta2(:,ind);
    V=V(ind,:); X=X(:,:,ind); I=I(:,ind); res=res(:,ind);Ti=Ti(ind,:);
    nchi3=size(Ti,1);   %trimmed cross-sectional dimension to compute the statistics
    %Re-definition of matrices tau and Tij
    [Tij,tau]=Tijmat(I,Ti,nchi3,T);
    %Re-definiton of cross-sectional estimator
    [nu2,~]=CR_reg_cond(beta1,beta2,V,diag(tau),K,p,q,nchi3,d1,d2,d,[],0);
    %Definition of the weights for the stat
    [~,w]=UnbiasedEst_Cond(nu2,vecbeta3,V,K,p,q,nchi3,T,diag(tau),d1,d2,d,[],0);
    [stat,pv,stat_a,pv_a,~]=CondStat(beta1,beta2,nu2,V,X,w,I,res,tau,Tij,p,q,K,d1,d2,d,nchi3,T,chi3);
end;

% Time variation test statistic (Section 4.3)
e=[zeros(p-1,1),eye(p-1,p-1)];
A=kron(eye(K),e);
% H_{0}: Avec(F')=0
gammaF=A*vec(Fhat');
XiF=T*gammaF'*(inv(A*SIGMAF*A'))*gammaF;    
pvF=1-chi2cdf(XiF,K*(p-1));
% H_{0}: Anu=0  %nu is constant over time    
% H_{0}: nu=0   %nu is zero over time    
if dataset~=3;
    gammaNU=A*nu;
    XiNU=T*gammaNU'*(inv(A*SIGMAnu*A'))*gammaNU;
    pvNU=1-chi2cdf(XiNU,K*(p-1));    
    
    XiNU2=T*nu'*(inv(SIGMAnu))*nu;
    pvNU2=1-chi2cdf(XiNU2,K*p); 
else
    gammaNU=A*nuB;
    XiNU=nchi*T*gammaNU'*(inv(A*SIGMAnu*A'))*gammaNU;
    pvNU=1-chi2cdf(XiNU,K*(p-1));
    
    XiNU2=nchi*T*nuB'*(inv(SIGMAnu))*nuB;
    pvNU2=1-chi2cdf(XiNU2,K*p); 
end;
   
%RESULTS
if dataset~=3;
    Table3=[vec(Fhat'),CIvecF,nu,CInu].*12*100/10;      %annualized vec(F'), nu
else
    Table3=[vec(Fhat'),CIvecF,nuB,CInu].*12*100/10;     %annualized vec(F'), nu
end;
Table5=[stat,pv;stat_a,pv_a];                      %statisticsoRiskPremia=lambda.*12*100/10;                      
RiskPremia=lambda.*12*100/10;                      %annualized risk premia
CIbRiskPremia=CIb.*12*100/10;                      %annualized CI below of risk premia
CIuRiskPremia=CIu.*12*100/10;                      %annualized CI up of risk premia
nuT=nuT.*12*100/10;                                %annualized time-varying nu  Kx(T_1)
CIbnuT=CIbnuT.*12.*100/10;                         %annualized CI below of nu_{t}
CIunuT=CIunuT.*12.*100/10;                         %annualized CI up of nu_{t}
TimeVarTest=[XiF,pvF;XiNU,pvNU;XiNU2,pvNU2];       %time variation test stat

