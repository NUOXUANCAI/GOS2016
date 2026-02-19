function [Tij, tau]=Tijmat(I,Ti,n,T)
%This function computes matrices Tij and tau 
Tij=zeros(n,n);     
for i=1:n-1;
    for j=i+1:n;    %symmetric matrix
        Tij(i,j)=sum(I(:,i).*I(:,j));  
        Tij(j,i)=Tij(i,j);
    end;   
end;
Tij=Tij+diag(Ti);
tau=T.*ones(n,n)./Tij;
ind=isinf(tau);     %tau is inf if Tij=0
tau(ind)=0;



