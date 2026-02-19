function [CInuB,CIlambda,tau,Tij]=ciLambdaNu_UncondnTinfty(res,b,w,v,nu,lambda,F,I,Ti,K,n,T,chi2)
%This function computes the CI of unbiased nu and lambda estimators
%We apply the asymptotic results for n and T\rightarrow\infty for the 
%time-invariant model (See Sections 3.3, 3.4)

% Matrices tau and Tij
[Tij,tau]=Tijmat(I,Ti,n,T); 

%% Confidence Intervals for unbiased nu
% Thresholding coefficients (we define M using a cross-validation method as
% proposed by Bickel and Levina (2008), see Supplementary Mat.)
if K==1;     
    M=0.0780;
else
    M=0.0570;
end;
CInuB=CInuUnbiased(res,I,tau,Tij,F,K,b,w,v,nu,n,T,M,chi2);

%% Confidence Intervals for lambda 
SigmaF=varcovNW(T,3,F);         %Var/Cov matrix of factors
varL=diag(SigmaF)/T;            %Asymptotic variance of lambda estimator
CIlambda=[lambda-1.96*sqrt(varL),lambda+1.96*sqrt(varL)];  


