function [nu,vecbeta3]=CR_reg_cond(beta1,beta2,V,tau,K,p,q,n,d1,d2,d,MC,weights)
%This function computes the nu WLS estimator - Equation (13) for the
%time-varying model                                
D=duplication(p);       %(p*p)x(p*(p+1)/2)
Dplus=(inv(D'*D))*D';   %Moore-Penrose inverse of D, (p*(p+1)/2)x(p*p)
Wpq=commutation(p,q);   %commutation matrix, (p*q)x(p*q)
%Definition of nu1 (OLS estimator) and vecbeta3
sum1=zeros(K*p,K*p);
sum2=zeros(K*p,1);
vecbeta3=zeros(d1*K*p,n);
for i=1:n;
    Bi=(invVec(beta2(1:K*p,i),p,K))';     %Kxp
    Ci=(invVec(beta2(K*p+1:end,i),q,K))'; %Kxq
    beta3=[(Dplus*kron(Bi',eye(p)));(Wpq*kron(Ci',eye(p)))]; %d1xKp
    sum1=sum1+beta3'*beta3;               %KpxKp
    sum2=sum2+beta3'*beta1(:,i);          %Kpx1
    vecbeta3(:,i)=vec(beta3');            %Kpd1x1 vec(beta3 TRANSPOSE)
end;

if weights==0;                  %WLS weights
    nu1=sum1\sum2;                  %Kpx1, OLS estimator
    [C,~]=Cmatrix(nu1,p,q,K,d1,d2); %C is dxd1 matrix - it depends on nu1
    %Definition of WLS estimator, eq.(13)
    sum1=zeros(K*p,K*p);
    sum2=zeros(K*p,1);
    D=duplication(d);               %d*dx(d*(d+1)/2)
    for i=1:n
        vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
        W=invVec(vecV,d,d);     %dxd
        v=tau(i,1)*C'*W*C;      %d1xd1 symmetric matrix
        w=eye(d1,d1)/(diag(diag(v)));           %d1xd1 weighting matrix  inv(diag(diag(v))
        beta3=(invVec(vecbeta3(:,i),K*p,d1))';  %d1xKp
        sum1=sum1+beta3'*w*beta3;
        sum2=sum2+beta3'*w*beta1(:,i);
    end;
    nu=sum1\sum2;

elseif weights==1;      %Value-weighted estimator
    MCb=(mean(MC).*tau')';
    %Definition of WLS estimator, eq.(13)
    sum1=zeros(K*p,K*p);
    sum2=zeros(K*p,1);
    for i=1:n        
        w=MCb(i,1)./sum(MCb);                  %scalar, weight
        beta3=(invVec(vecbeta3(:,i),K*p,d1))';  %d1xKp
        sum1=sum1+beta3'*w*beta3;
        sum2=sum2+beta3'*w*beta1(:,i);
    end;
    nu=sum1\sum2;

elseif weights==2;      %OLS estimator
    nu=sum1\sum2;                  %Kpx1, OLS estimator
    
end;
