function [stat,pv,stat_a,pv_a,SIGMA_xi]=UncStat(a,b,F,K,nu,V,w,v,I,res,n,T,tau,Tij,chi3)
%This function computes the statistic and its p-values for the two null
%hypothesis in Section 3.5.
%We compute statistics based on Q_{e} and on Q_{a} for the
%time-invariant model  when n and T \rightarrow \infty
beta=[a';b'];       %(K+1)xn, coefficients
c=[1;-nu];          %(K+1)x1
X=[ones(T,1),F];    %regressor matrix
Qx=(X'*X)/T;        %(K+1)x(K+1)
D=duplication(K+1); %duplication matrix
E1=[1;zeros(K,1)];  %selection matrix

%vector of weights - for stat based on Qa
va=zeros(n,1);
wa=zeros(n,1);      
for i=1:n;
    vecV=D*V(i,:)';             %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
    W=invVec(vecV,K+1,K+1);     %W=inv(Qx)*Sii*inv(Qx)- matrix(K+1)x(K+1)
    va(i,1)=tau(i,i)*E1'*W*E1;  %scalar
    wa(i,1)=1/va(i,1);          %scalar
end;

Wn=diag(w);          %diagonal matrix with elements w
Wn_a=diag(wa);       %diagonal matrix with elements wa
Vn=zeros(n,n);       %variance matrix with elements vij
Vn_a=zeros(n,n);     %variance matrix with elements vaij

% Thresholding coefficients (we define M using a cross-validation method as proposed by Bickel and Levina (2008) 
if K==1;     
    M=0.0780;
else
    M=0.0570;
end;
eta=1;
%Definiton of \tilde{S_{ij}}, and elements of Vn and Vn_a
for i=1:n;
    Vn(i,i)=v(i,1);     %diagonal elements of Vn    
    Vn_a(i,i)=va(i,1);
    for j=i+1:n;        %Sij tilde is a symmetric matrix
        if tau(i,j)<=chi3 && tau(i,j)~=0; 
            Iij=I(:,i).*I(:,j);             %Indicator Iij
            ri=repmat(Iij.*res(:,i),1,K+1).*X;%ri=Iij.*(res(:,i); 
            rj=repmat(Iij.*res(:,j),1,K+1).*X;%rj=Iij.*res(:,j);
            Sij=(ri'*rj)/(Tij(i,j));        %\hat{Sij}
            kappa=M*sqrt(log(n)./(Tij(i,j).^eta));  %thresholding parameter
            THR=norm(Sij,'fro')>=kappa;     %it is a scalar with 1 and 0 computed on the correlation matrix
            Sij=THR.*Sij;                   %\tilde{Sij}

            vij=(tau(i,i)*tau(j,j)/tau(i,j))*c'*(inv(Qx))*Sij*(inv(Qx))*c; %see def. vij 
            Vn(i,j)=vij;
            Vn(j,i)=Vn(i,j);                
            
            vaij=(tau(i,i)*tau(j,j)/tau(i,j))*E1'*(inv(Qx))*Sij*(inv(Qx))*E1; %see def. vij 
            Vn_a(i,j)=vaij;
            Vn_a(j,i)=Vn_a(i,j);            
        end;
    end;        
end;

%Statistic based on Qe
e=(c'*beta)';                   %nx1
Qe=(e'*Wn*e)/n;                 %Qe
Bxi=1;                          %bias term of statistic xi 
KXi=T*sqrt(n)*(Qe-Bxi/T);       %chi
SIGMA_xi=2*(w'*(Vn.^2)*w)/n;    %variance of stat 
stat=(SIGMA_xi)^(-1/2)*KXi;     %statistic based on Qe
pv=(1-normcdf(stat,0,1));       %p-value

%Statistic based on Qa
Qa=(a'*Wn_a*a)/n;               %Qa
Bxi_a=1;                        %bias term of statistic xi 
KXi_a=T*sqrt(n)*(Qa-Bxi_a/T);   %chi_a
SIGMA_xia=2*(wa'*(Vn_a.^2)*wa)/n;%variance of stat
stat_a=(SIGMA_xia)^(-1/2)*KXi_a;%statistic based on Qa    
pv_a=(1-normcdf(stat_a,0,1));   %p-value



