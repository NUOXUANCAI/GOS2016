function [nuB,w,v,Bnu]=UnbiasedEst(nu,b,V,Vxi,tau,n,T,K,MC,weights)
% This function computes the unbiased estimator of nu for the
% time-invariant model - see Prop. 4

c=[1;-nu];              %(K+1)x1, nu is the WLS estimator
Bs1=zeros(K,1);         %sum in the bias_nu
v=zeros(n,1);           %vector of vii
w=zeros(n,1);           %vector of weigths w
E2=[zeros(1,K); eye(K)];%selection matrix E2
D=duplication(K+1);     %duplication matrix

if weights==1;
    MCb=(mean(MC).*tau')';
end;
%Definition of internal sum of Bias term
for i=1:n
    vecV=D*V(i,:)';             %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
    W=invVec(vecV,K+1,K+1);     %W=inv(Qx)*Sii*inv(Qx)- matrix(K+1)x(K+1)
    v(i,1)=tau(i,1)*c'*W*c;     %scalar
    if weights==0;              %WLS weights
        w(i,1)=1/v(i,1);        %scalar, weight
    elseif weights==1;          %VW weights
        w(i,1)=MCb(i,1)./sum(MCb);
    elseif weights==2;          %OLS weights
        w(i,1)=1;
    end;
    vecVxi=D*Vxi(i,:)';
    Wxi=invVec(vecVxi,K+1,K+1); %Wxi=inv(Qxi)*Sii*inv(Qxi)- matrix(K+1)x(K+1)
    Bs1=Bs1+w(i,1)*tau(i,1)*E2'*Wxi*c; %sum in the bias term of nu - Prop.4
end;

Wn=diag(w);            %nxn diagonal matrix
Qb=(b'*Wn*b)/n;        %weighted second moment of estimated b
Bnu=Qb\(Bs1/n);        %Bias term
nuB=nu-Bnu/T;          %unbiased estimator of nu

