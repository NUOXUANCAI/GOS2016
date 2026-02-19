function [CInu,CIlambda,tau,Tij]=ciLambdaNu_Uncondnfix(res,V,b,nu,lambda,F,MC,I,Ti,K,n,T,weights)
%This function computes the CI of nu and lambda. 
%We apply the asymptotic results for n fixed and T\rightarrow\infty for the 
%time-invariant model (See Sections 4.3, 4.4)

%Tij and tau matrix (nxn matrices)
[Tij,tau]=Tijmat(I,Ti,n,T);

%% Variance-Cov of factors
SigmaF=varcovNW(T,3,F);

%% Variance-Cov of estimators nu - Under heteroskedasticy assumption
% see Section 3.4
X=[ones(T,1),F];    %regressor matrix
Qx=(X'*X)/T;        %(K+1)x(K+1)
Bn=b;               %nxK matrix of factor loadings
c=[1;-nu];          %(K+1)x1
D=duplication(K+1); %duplication matrix
w=zeros(n,1);       %weights 
Vn=zeros(n,n);      %symmetric matrix with element vij
for i=1:n;
        vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
        W=invVec(vecV,K+1,K+1); %(K+1)x(K+1)
        v=tau(i,i)*c'*W*c;      %scalar
    if weights==0;              %WLS weights
        w(i,1)=1/v;             %scalar
    elseif weights==1;          %VW weights
        mcB=tau(i,1)*mean(MC(:,i));      %scalar
        w=mcB;
    elseif weights==2;          %OLS weights
        w=1;
    end;        
    Vn(i,i)=v;              %diagonal element
    %Definition of Sij, and elemnts of Vn
    for j=i+1:n;            %symmetric matrix
        if tau(i,j)~=0;     
            sum1=zeros(K+1,K+1);
            for t=1:T;
                sum1=sum1+res(t,i)*res(t,j)*X(t,:)'*X(t,:);
            end;
            Sij=sum1/(Tij(i,j)); 
            v=(tau(i,i)*tau(j,j)/tau(i,j))*c'*(inv(Qx))*Sij*(inv(Qx))*c;
            Vn(i,j)=v;      
            Vn(j,i)=v;
        end;
    end;
end;
Wn=diag(w);                 %diagonal matrix of weights
SigmaNu_N=(inv((Bn'*Wn*Bn)/n))*((Bn'*Wn*Vn*Wn*Bn)/(n^2))*(inv((Bn'*Wn*Bn)/n)); 
var_nu_N=diag(SigmaNu_N)/T;
st_nu_N=sqrt(var_nu_N);                             %stand.dev of nu
CInu=[nu-1.96*st_nu_N, nu+1.96*st_nu_N];            %conf.int. for nu

%% Variance-Cov of risk premia 
var_L_N=diag(SigmaF+SigmaNu_N)/T;                   
st_L_N=sqrt(var_L_N);                               %stand.dev of lambda
CIlambda=[lambda-1.96*st_L_N,lambda+1.96*st_L_N];   %conf.int. for lambda
