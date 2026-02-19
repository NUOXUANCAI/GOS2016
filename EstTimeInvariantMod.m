function [n,nchi,nchi3,Table1,Table2,Table4,IC]=EstTimeInvariantMod(n,T,K,ExcR,MC,F,I,Ti,dataset)                                                
% This function estimates time-invariant risk premia. 
% OUTPUT: 
% n=cross-sectional dimension; 
% nchi=trimmed cross-sectional dimension used in the estimation approach;
% nchi3=trimmed cross-sectional dimension used to compute statistics;
% Table1 is a Kx3 matrix that contains estimated (bias corrected) annualized risk
% (col 1) premia and their confidence intervals (col 2,3);
% Table2 is a Kx3 matrix that contains estimated (bias corrected) annualized nu
% (col 1) and their confidence intervals (col 2,3);
% Table4 is a 2x2 matrix that contains statistics (col 1) based on Q_{e}
% (1st row) and on Q_{a} (2nd row) and their p-values (col 2)

%First pass             (see Section 3.2)
%Time-series regression
[a,b,res,CN,V,Vxi,tau,~,~,~,~]=TSRegress_CN_V(ExcR,F,I,Ti,n,T,K);

%Trimming Approach
chi1=15;
chi2=T/12;
ind=find(CN<chi1);
a=a(ind,:); b=b(ind,:); V=V(ind,:); Vxi=Vxi(ind,:); I=I(:,ind); res=res(:,ind); tau=tau(ind,:); Ti=Ti(ind,:); 
ind=(tau<=chi2);    
a=a(ind,:); b=b(ind,:); V=V(ind,:); Vxi=Vxi(ind,:); I=I(:,ind); res=res(:,ind); tau=tau(ind,:); Ti=Ti(ind,:); 

nchi=size(Ti,1);    %trimmed cross-sectional dimension

%Second pass            (see eq. (14))
%Cross-sectional regression
nu=CR_reg(a,b,V,tau,nchi,K,[],0);

if dataset~=3;      %Portfolios: n small, T large         
    %Risk premia estimator eq. (14)
    lambda=nu+mean(F)';
    %Confidence intervals of estimators nu and lambda
    [CInu,CIlambda,tau,Tij]=ciLambdaNu_Uncondnfix(res,V,b,nu,lambda,F,[],I,Ti,K,nchi,T,0);
    %Statistics: Gibbons, Ross and Shanken (1989) stat.
    nchi3=nchi;        %no trimming for panel of portfolios
    [stat,pv,stat_a,pv_a]=GRS_Uncstatistics(tau,Tij,F,V,res,a,b,nchi3,T,K);
else                %Individual stocks: n large, T large
    %Unbiased estimator of nu (see Prop. 4, 5)
    [nu,w,v,~]=UnbiasedEst(nu,b,V,Vxi,tau,nchi,T,K,[],0);
    %Unbiased estimator of risk premia 
    lambda=nu+mean(F)';   
    %Confidence intervals of estimators nu and lambda (see Section 4.4)
    [CInu,CIlambda,tau,Tij]=ciLambdaNu_UncondnTinfty(res,b,w,v,nu,lambda,F,I,Ti,K,nchi,T,chi2);
    %Statistics: Section 3.5
    % Trimming (see Section 4.3 and results on Monte Carlo simulations in supplementary materials)
    chi3=T/240;
    ind=find(diag(tau)<=chi3);
    a=a(ind,:); b=b(ind,:); V=V(ind,:); I=I(:,ind); res=res(:,ind); tau=tau(ind,ind); Tij=Tij(ind,ind);
    nchi3=size(tau,1);  %trimmed cross-sectional dimensinìon to compute the statistics
    %Cross-sectional regression on the trimmed sample 
    nu2=CR_reg(a,b,V,diag(tau),nchi3,K,[],0);    
    %Definition of the weights on the trimmed sample 
    [v,w,V]=weightsSTAT(res,I,diag(Tij),nu2,F,diag(tau),nchi3,T,K);    
    [stat,pv,stat_a,pv_a,~]=UncStat(a,b,F,K,nu2,V,w,v,I,res,nchi3,T,tau,Tij,chi3);
end;


%RESULTS
Table1=[lambda,CIlambda].*12*100/10;   %annualized risk premia
Table2=[nu,CInu].*12*100/10;           %annualized nu
Table4=[stat,pv;stat_a,pv_a];          %statistics by GRS


