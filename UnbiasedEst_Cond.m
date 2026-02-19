function [nuB,w]=UnbiasedEst_Cond(nu,vecbeta3,V,K,p,q,n,T,tau,d1,d2,d,MC,weights)
% This function computes the unbiased estimator of nu for the
% time-varying model - see Prop. 4

[C,Ja]=Cmatrix(nu,p,q,K,d1,d2);     %C is dxd1 and Ja is pKd1xd2
Jb=kron(vec(eye(d1))',eye(K*p))*kron(eye(d1),Ja); %Jb is kpxd1d2
E2=[zeros(d1,d2); eye(d2)];         %selction matrix, dxd2
if weights==0;
    w=zeros(d1,n);              %diagonal elements of weighting matrix
elseif weights==1;
    w=zeros(1,n);
elseif weights==2;
    w=ones(1,n);
end;
D=duplication(d);           %d*dxd*(d+1)/2
Bs1=zeros(d1*d2,1);         %term of sum of bias_nu
Qb3W=zeros(K*p,K*p);        %weighted second moment of beta3

if weights==0;              %WLS weights
    for i=1:n
        vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
        W=invVec(vecV,d,d);     %dxd
        v=tau(i,1)*C'*W*C;      %d1xd1 symmetric matrix
        w(:,i)=diag(v);         %d1x1
        wi=eye(d1,d1)/(diag(w(:,i)));   %d1xd1 weighting matrix
        Bs1=Bs1+tau(i,1)*vec(E2'*W*C*wi);
        beta3=(invVec(vecbeta3(:,i),K*p,d1))';  %d1xKp
        Qb3W=Qb3W+beta3'*wi*beta3;              %KpxKp
    end;

elseif weights==1;          %VW weights
    MCb=(mean(MC).*tau')';
    for i=1:n
        vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
        W=invVec(vecV,d,d);     %dxd
        w(:,i)=MCb(i,1)./sum(MCb);            %scalar, weight        
        Bs1=Bs1+tau(i,1)*vec(E2'*W*C*w(:,i));
        beta3=(invVec(vecbeta3(:,i),K*p,d1))';  %d1xKp
        Qb3W=Qb3W+beta3'*w(:,i)*beta3;              %KpxKp
    end;
    
elseif weights==2;          %OLS weights
    for i=1:n
        vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
        W=invVec(vecV,d,d);     %dxd        
        Bs1=Bs1+tau(i,1)*vec(E2'*W*C);
        beta3=(invVec(vecbeta3(:,i),K*p,d1))';  %d1xKp
        Qb3W=Qb3W+beta3'*beta3;              %KpxKp
    end;
    
end;

Qb3W=Qb3W/n;                %Qb3W
Bnu=Qb3W\Jb*(Bs1/n);        %(inv(Qb3W))*Jb*(Bs1/n)
nuB=nu-Bnu/T;    
