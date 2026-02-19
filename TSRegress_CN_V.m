function [a,b,res,CN,V,Vxi,tau,R2,R2adj,SYSRisk,IDIRisk]=TSRegress_CN_V(Y,F,I,Ti,n,T,K)
%This function estimates the time-invariant betas by TS regression
%see Section 3.2

X=[ones(T,1),F];                    %regressor matrix, Tx(K+1)
Qx=X'*X/T;                          %Qx
tau=T./Ti;                          %tau
a=zeros(n,1);                       %intercept, nx1
b=zeros(n,K);                       %factor loadings, nxK
CN=zeros(n,1);                      %condition number for each asset i 
V=zeros(n,(K+1)*(K+1+1)/2);         %V stores (vech((\hatQ^-1)Sii(\hatQ^-1)))' for each i -
Vxi=zeros(n,(K+1)*(K+1+1)/2);       %Vxi stores (vech((\hatQ_{xi}^-1)Sii(\hatQ_{xi}^-1)))' for each i - 
res=zeros(T,n);                     %matrix of residuals, Txn

R2=zeros(n,1);                      %coeffients of determination - Definitions in Appendix 11 SM
R2adj=zeros(n,1);                   %adjusted coefficients of determination
SYSRisk=zeros(n,1);                 %systematic risk
IDIRisk=zeros(n,1);                 %idiosyncratic risk

for i=1:n;   
    Ii=repmat(I(:,i),1,K+1);        %indicator matrix Ii, Tx(K+1)
    Xi=Ii.*X;                       %regressor matrix Xi
    Qxi=(Xi'*Xi)/Ti(i);             %Qxi, (K+1)x(K+1)
    %Beta - OLS estimators
    beta=(inv(Qxi))*(Xi'*Y(:,i)/Ti(i));  %(K+1)x1
    a(i,:)=beta(1,1);               %1x1
    b(i,:)=beta(2:K+1,1);           %Kx1
    res(:,i)=Y(:,i)-Xi*beta;        %Fitted residual    
    
    %Empirical measures of idiosyncratic risk (Appendix 11 SM)
    YFit=Xi*beta;                   
    mYFit=sum(YFit)/Ti(i);
    mY=sum(Y(:,i))/Ti(i);
    TSS=sum((Y(:,i)-repmat(mY,T,1)).^2);
    ESS=sum((YFit-repmat(mYFit,T,1)).^2);
    RSS=sum(res(:,i).^2);
    R2(i,1)=ESS/TSS;
    R2adj(i,1)=1-(Ti(i)-1)/(Ti(i)-(K+1))*(1-R2(i,1));
    SYSRisk(i,1)=sqrt(ESS/Ti(i));
    IDIRisk(i,1)=sqrt(RSS/Ti(i));
        
    %Condition number - Greene (2008)
    %X is scaled to have unit column length
    Xtildei=zeros(T,(K+1));     
    for s=1:(K+1);
        Xtildei(:,s)=Xi(:,s)/norm(Xi(:,s)); 
    end;
    Qxtildei=(Xtildei'*Xtildei)/Ti(i);  
    %We compute CN_{i} only if Qxtildei is a p.d. matrix 
    l=eig(Qxtildei);
    if l>0;                     
        CN(i,1)=sqrt(max(l)/min(l));
    else                          %If the matrix is not p.d.
        CN(i,1)=100;
    end;
    %Defintiion of Sii matrix, V and Vxi
    ri=repmat(res(:,i),1,K+1).*Xi;  
    Sii=(ri'*ri)/Ti(i);            %(K+1)x(K+1)
    Vmat=(inv(Qx))*Sii*(inv(Qx));
    V(i,:)=(vech(Vmat))';    
    Vmat=(inv(Qxi))*Sii*(inv(Qxi));
    Vxi(i,:)=(vech(Vmat))';
end;
