function nu=CR_reg(a,b,V,tau,n,K,MC,weights)
%This function computes the nu WLS estimator 

if weights==0;          %WLS weights - Equation (14)    
    nu1=regress(a,b);   %OLS estimator, Kx1
    c=[1;-nu1];         %(K+1)x1
    sum1=zeros(K,K);    %KxK
    sum2=zeros(K,1);    %Kx1
    D=duplication(K+1);
    for i=1:n
        vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
        W=invVec(vecV,K+1,K+1); %(K+1)x(K+1), inv(Qx)*Sii*inv(Qx)
        v=tau(i,1)*c'*W*c;      %scalar
        w=1/v;                  %scalar, weight
        sum1=sum1+w*b(i,:)'*b(i,:);    %KxK
        sum2=sum2+w*b(i,:)'*a(i,:);    %Kx1
    end;
    nu=sum1\sum2;            %WLS of nu

elseif weights==1;      %Value-weighted estimator
    MCb=(mean(MC).*tau')';
    sum1=zeros(K,K);    %KxK
    sum2=zeros(K,1);    %Kx1
    for i=1:n
        w=MCb(i,1)./sum(MCb);          %scalar, weight
        sum1=sum1+w*b(i,:)'*b(i,:);    %KxK
        sum2=sum2+w*b(i,:)'*a(i,:);    %Kx1
    end;
    nu=sum1\sum2;            %WLS of nu

elseif weights==2;      %OLS estimator
    nu=regress(a,b);   
end;


