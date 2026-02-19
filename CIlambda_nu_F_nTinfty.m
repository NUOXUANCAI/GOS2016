function [CIvecF,CInu,CIb,CIu,CIbnuT,CIunuT,tau,SIGMAF,SIGMAnu]=CIlambda_nu_F_nTinfty(Ti,I,res,V,X,F,Z,MC,vecbeta3,nu,nuT,lambda,Fhat,n,T,K,p,q,d1,d2,d,chi2,weights)
                                                                                      
%This function computes the CI of nu and lambda. 
%We apply the asymptotic results for n and T\rightarrow\infty for the 
%time-varying model (See Sections 3.3, 3.4)

%% Var-cov of vec(F')
u=zeros(T,K);                           %TxK - matrix of fitted residuals
u(2:T,:)=F(2:T,:)-Z(1:T-1,:)*Fhat';     %Fitted residuals of SUR
su=zeros(K*p,K*p);
for t=2:T;
    su=su+kron(u(t,:)'*u(t,:),Z(t-1,:)'*Z(t-1,:));
end;
SIGMAu=su/(T-1);                        %KpxKp
Qz=(Z(1:T-1,:)'*Z(1:T-1,:))/(T-1);      %pxp

SIGMAF=(kron(eye(K),(inv(Qz))))*SIGMAu*(kron(eye(K),(inv(Qz))));
stF=sqrt(diag(SIGMAF)/T);
CIvecF(:,1)=vec(Fhat')-1.96*stF;
CIvecF(:,2)=vec(Fhat')+1.96*stF;

%% Var-cov of nu
% Tij and tau matrix
[Tij,tau]=Tijmat(I,Ti,n,T-1);
% Definition of C
[C,~]=Cmatrix(nu,p,q,K,d1,d2);    %dxd1 matrix
% Definition of Qb3
D=duplication(d);               %d*dxd*(d+1)/2
Qb3=zeros(K*p,K*p);             %Qb3 matrix
Qx=zeros(d,d,n);             %Qx_{i}
beta3=zeros(d1,K*p,n);       %beta3_{i}
if weights==0;               %WLS weights
    w=zeros(d1,d1,n);       
elseif weights==1;           %VW weights
    w=zeros(1,1,n);
    MCb=(mean(MC).*(diag(tau))')';
elseif weights==2;           %OLS weights
    w=ones(1,1,n);
end;

for i=1:n
    beta3(:,:,i)=(invVec(vecbeta3(:,i),K*p,d1))';  %d1xKp    
    if weights==0;
        vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
        W=invVec(vecV,d,d);     %dxd
        v=tau(i,i)*C'*W*C;      %d1xd1 symmetric matrix
        w(:,:,i)=eye(d1,d1)/(diag(diag(v)));  %d1xd1 weighting matrix    
    elseif weights==1;        
        w(:,:,i)=MCb(i,1)./sum(MCb);                  %scalar, weight
    end;
    Qx(:,:,i)=(X(:,:,i)*X(:,:,i)')/Tij(i,i);
    Qb3=Qb3+beta3(:,:,i)'*w(:,:,i)*beta3(:,:,i);%KpxKp matrix
end;
Qb3=Qb3/n;

Ssum=zeros(K*p,K*p);
% Thresholding coefficients (we define M using a cross-validation method as
% proposed by Bickel and Levina (2008), see Supplementary Mat.)
if K==1;      
    M=0.0750;  
elseif K==3;
    M=0.0580;
else
    M=0.0670;
end;
eta=1;
%Definition of \tilde{Sigma_{\nu}} 
for i=1:n;
    %display(i);
    for j=i:n;       %symmetric property
        if i~=j && tau(i,j)~=0 && tau(i,j)<=chi2;
            Iij=I(:,i).*I(:,j);             %Indicator Iij
            ri=repmat(Iij.*res(:,i),1,d).*X(:,:,i)';    %(T-1)xd
            rj=repmat(Iij.*res(:,j),1,d).*X(:,:,j)';    %(T-1)xd
            Sij=(ri'*rj)/(Tij(i,j));       %dxd matrix
            kappa=M*sqrt(log(n)./(Tij(i,j).^eta));   %kappa depends on Tij
            THR=norm(Sij)>=kappa;          %it takes value 1 or 0 - comparison on the correlation            
            if THR==1;
                Sijtilde=THR.*Sij;          %thresholded covariance values
                Sjitilde=Sijtilde';
                
                CSQ_ijC=C'*(inv(Qx(:,:,i)))*Sijtilde*(inv(Qx(:,:,j)))*C;
                CSQ_jiC=C'*(inv(Qx(:,:,j)))*Sjitilde*(inv(Qx(:,:,i)))*C;
                
                Ssum=Ssum+(tau(i,i)*tau(j,j)/tau(i,j))*beta3(:,:,i)'*w(:,:,i)*CSQ_ijC*w(:,:,j)*beta3(:,:,j)+...
                    (tau(i,i)*tau(j,j)/tau(i,j))*beta3(:,:,j)'*w(:,:,j)*CSQ_jiC*w(:,:,i)*beta3(:,:,i);
            end;
            
        elseif i==j
            vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
            CSQ_iiC=C'*invVec(vecV,d,d)*C;     %dxd
            Ssum=Ssum+tau(i,i)*beta3(:,:,i)'*w(:,:,i)*CSQ_iiC*w(:,:,i)*beta3(:,:,i);
        end;
    end;
end;
SIGMAnu=(inv(Qb3))*Ssum/(n)*(inv(Qb3));      
stnu=sqrt(diag(SIGMAnu)/(n*(T-1)));
CInu(:,1)=nu-1.96*stnu;
CInu(:,2)=nu+1.96*stnu;

%% Var-cov of time-varying lambda and nu
SIGMALAMBDA=SIGMAF;                     %See Prop.4
stlambda=zeros(K,T);                    %Standard error for any lambda_t - KxT
CIb=zeros(K,T-1);                       %below
CIu=zeros(K,T-1);                       %up
stnuT=zeros(K,T);                       %Standard error for any nu_t - KxT
CIbnuT=zeros(K,T-1);                    %CI below
CIunuT=zeros(K,T-1);                    %CI up
Wpk=commutation(p,K);
for t=2:T;
    SIGMAlambda=(kron(Z(t-1,:),eye(K)))*Wpk*SIGMALAMBDA*Wpk'*(kron(Z(t-1,:)',eye(K)));
    stlambda(:,t-1)=sqrt(diag(SIGMAlambda)/(T-1));
    CIb(:,t-1)=(lambda(:,t-1)-1.96*stlambda(:,t-1));
    CIu(:,t-1)=(lambda(:,t-1)+1.96*stlambda(:,t-1));
    SIGMAnuT=(kron(Z(t-1,:),eye(K)))*Wpk*SIGMAnu*Wpk'*(kron(Z(t-1,:)',eye(K)));
    stnuT(:,t-1)=sqrt(diag(SIGMAnuT)/(n*(T-1)));
    CIbnuT(:,t-1)=(nuT(:,t-1)-1.96*stnuT(:,t-1));
    CIunuT(:,t-1)=(nuT(:,t-1)+1.96*stnuT(:,t-1));
end;
