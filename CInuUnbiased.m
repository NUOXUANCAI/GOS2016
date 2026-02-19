function CInuB=CInuUnbiased(res,I,tau,Tij,F,K,b,w,v,nuB,n,T,M,chi2)
%This function computes the CI for the unbiased estimator nu
%We apply the asymptotic results obtained for n and T \righarrow\infty.
%see Prop. 4 and 5, and eq.(15)

Bn=b;               %nxK matrix of factor loadings
c=[1;-nuB];         %(K+1)x1
X=[ones(T,1),F];    %regressor matrix
Qx=(X'*X)/T;        %(K+1)x(K+1)
Wn=diag(w);         %diagonal matrix of weights
Vn=zeros(n,n);      %variance matrix with elements vij
eta=1;              
%\tilde{S_{ij}} and elements of Vn
for i=1:n;
    Vn(i,i)=v(i,1); %diagonal elements of Vn 
    for j=i+1:n;    %symmetric property
        if tau(i,j)~=0 && tau(i,j)<=chi2; 
            Iij=I(:,i).*I(:,j);             %Indicator Iij
            ri=repmat(Iij.*res(:,i),1,K+1).*X;
            rj=repmat(Iij.*res(:,j),1,K+1).*X;
            Sij=(ri'*rj)/(Tij(i,j));        %\hat{Sij}
            kappa=M*sqrt(log(n)./(Tij(i,j).^eta));  %thresholding parameter
            THR=norm(Sij,'fro')>=kappa;     %THR is a matrix with 1 and 0 computed on the correlation matrix
            Sij=(THR.*Sij);                 %\tilde{Sij}
            vij=(tau(i,i)*tau(j,j)/tau(i,j))*c'*(inv(Qx))*Sij*(inv(Qx))*c;  %see def. vij 
            Vn(i,j)=vij;
            Vn(j,i)=Vn(i,j);                
        end;
    end;        
end;

SIGMAnu=(inv((Bn'*Wn*Bn)/n))*((Bn'*Wn*Vn*Wn*Bn)/n)*(inv((Bn'*Wn*Bn)/n));    %eq. (15)
varN=diag(SIGMAnu)/(T*n);
CInuB(:,1)=nuB-1.96*sqrt(varN);
CInuB(:,2)=nuB+1.96*sqrt(varN);
