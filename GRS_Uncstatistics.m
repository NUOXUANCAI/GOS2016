function [GRSe,pv1,GRSa,pv2]=GRS_Uncstatistics(tau,Tij,F,V,res,a,b,n,T,K)
%This function computes the Gibbons, Ross and Shanken (1989) statistics. 
%We consider invariant model and small cross-sectional dimension.

%Definition of Var-Cov matrix: SIGMA
SIGMA=zeros(n*(K+1),n*(K+1));
D=duplication(K+1);         %duplication matrix
X=[ones(T,1),F];            %regressor matrix
Qx=(X'*X)/T;                %Qx,(K+1)x(K+1)
p2=(K+1:K+1:n*(K+1))';      %p1,p2: indices to ref to blocks of matrix SIGMA
p1=(1:K+1:n*(K+1))';
for i=1:n;
    vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
    W=invVec(vecV,K+1,K+1); %(K+1)x(K+1)
    SIGMA(p1(i):p2(i),p1(i):p2(i))=tau(i,i)*W;  %diagonal block
    for j=i+1:n;            %symmetric matrix, no-diagonal block
        if tau(i,j)~=0;
            sum1=zeros(K+1,K+1);
            for t=1:T;
                sum1=sum1+res(t,i)*res(t,j)*X(t,:)'*X(t,:);
            end;
            Sij=sum1/(Tij(i,j));    
            v=(tau(i,i)*tau(j,j)/tau(i,j))*(inv(Qx))*Sij*(inv(Qx));
            SIGMA(p1(i):p2(i),p1(j):p2(j))=v;
            SIGMA(p1(j):p2(j),p1(i):p2(i))=v;
        end;
    end;
end;

%Definition of the OMEGA var-cov depending on nu_OLS estimator 
%Feasible GLS estimator
nu1=regress(a,b);   %nu OLS estimator
c=[1;-nu1];         %vector c
OMEGA=kron(eye(n),c')*SIGMA*kron(eye(n),c);     %OMEGA is n(K+1)xn(K+1)
nu=(inv(b'*(inv(OMEGA))*b))*(b'*(inv(OMEGA))*a);%FGLS estimator of nu
e=a-b*nu;           %fitted cross-sect errors of FGLS regression

%GRS stat H_{0}:e=0
GRSe=T*e'*(inv(OMEGA))*e;
pv1=1-chi2cdf(GRSe,n-K);

%GRS stat H_{0}:a=0
E1=[1;zeros(K,1)];
OMEGAa=kron(eye(n),E1')*SIGMA*kron(eye(n),E1);     %OMEGA is n(K+1)xn(K+1)
GRSa=T*a'*(inv(OMEGAa))*a;
pv2=1-chi2cdf(GRSa,n);