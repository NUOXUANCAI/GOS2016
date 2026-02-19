function [stat,pv,stat_a,pv_a,SIGMA_xi]=CondStat(beta1,beta2,nu,V,X,w,I,res,tau,Tij,p,q,K,d1,d2,d,n,T,chi3)
%This function computes the statistic and its p-values for the two null
%hypothesis in Section 3.5.
%We compute statistics based on Q_{e} and on Q_{a} for the
%time-varying model  when n and T \rightarrow \infty

[C,~]=Cmatrix(nu,p,q,K,d1,d2);    %dxd1 matrix
%Definition of Qe and Qa
E1=[eye(d1);zeros(d2,d1)];        %dxd1
D=duplication(d);           %d*dxd*(d+1)/2
wa=zeros(d1,n);             %diagonal elements of weighting matrix
ewe=zeros(n,1);             %it contains e'*wi*e at each row (Qe)
awa=zeros(n,1);             %it contains a'*wa*a at each row (Qa)
for i=1:n
    beta=[beta1(:,i);beta2(:,i)];   %dx1
    e=C'*beta;              %d1x1
    wi=eye(d1,d1)/(diag(w(:,i))); %d1xd1 weighting matrix where w(:,i) is diagonal elements of weighting matrix
    ewe(i,1)=e'*wi*e;       %scalar
    
    vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
    W=invVec(vecV,d,d);     %dxd
    va=tau(i,i)*E1'*W*E1;   %d1xd1 symmetric matrix
    wa(:,i)=diag(va);       %d1x1
    wai=eye(d1,d1)/(diag(wa(:,i)));%d1xd1 weighting matrix
    awa(i,1)=beta1(:,i)'*wai*beta1(:,i);    %scalar
end;
Qe=sum(ewe)/n;
Qa=sum(awa)/n;

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
% See Definition of \tilde{\Sigma}_{\xi} 
Vdiag=0;                     %diagonal term of variance of xi
Vcross=0;                    %cross term of variance of xi
Vdiag_a=0;                   %diagonal term of variance of xi (dep on Qa)
Vcross_a=0;                  %cross term of variance of xi (dep on Qa)
for i=1:n    
    %display(i);
    vecV=D*V(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
    W=invVec(vecV,d,d);     %dxd
    vi=C'*W*C;
    via=E1'*W*E1;
    wi=eye(d1,d1)/(diag(w(:,i)));
    wai=eye(d1,d1)/(diag(wa(:,i)));
    Vdiag=Vdiag+(tau(i,i)^2)*trace(wi*vi*wi*vi);
    Vdiag_a=Vdiag_a+(tau(i,i)^2)*trace(wai*via*wai*via);
    for j=i+1:n               %Sij tilde is different to Sji        
        if tau(i,j)~=0 && tau(i,j)<=chi3;
            Iij=I(:,i).*I(:,j);             %Indicator Iij            
            ri=repmat(Iij.*res(:,i),1,d).*X(:,:,i)';    %(T-1)xd
            rj=repmat(Iij.*res(:,j),1,d).*X(:,:,j)';    %(T-1)xd            
            Sij=(ri'*rj)/(Tij(i,j));       %dxd matrix
            Sji=(rj'*ri)/(Tij(j,i));       %dxd matrix
            kappa=M*sqrt(log(n)./(Tij(i,j).^eta));
            THR1=norm(Sij,'fro')>=kappa;    %THR is a matrix of 1 and 0 computed on the correlation matrix
            Sij=THR1.*Sij;                  %thresholded estimator
            THR2=norm(Sji,'fro')>=kappa;
            Sji=THR2.*Sji;
            
            if THR1~=0 && THR2~=0;                                       
                Qxi=(X(:,:,i)*X(:,:,i)')/Tij(i,i);
                Qxj=(X(:,:,j)*X(:,:,j)')/Tij(j,j);             
                
                vij=C'*(inv(Qxi))*Sij*(inv(Qxj))*C;
                vji=C'*(inv(Qxj))*Sji*(inv(Qxi))*C;
                
                vija=E1'*(inv(Qxi))*Sij*(inv(Qxj))*E1;
                vjia=E1'*(inv(Qxj))*Sji*(inv(Qxi))*E1;
                
                wi=eye(d1,d1)/(diag(w(:,i)));
                wj=eye(d1,d1)/(diag(w(:,j)));
                wai=eye(d1,d1)/(diag(wa(:,i)));
                waj=eye(d1,d1)/(diag(wa(:,j)));
                
                Vcross=Vcross+2*(tau(i,i)^2*tau(j,j)^2/tau(i,j)^2)*trace(wi*vij*wj*vji);
                Vcross_a=Vcross_a+2*(tau(i,i)^2*tau(j,j)^2/tau(i,j)^2)*trace(wai*vija*waj*vjia);
            end;
        end;
    end;
end

%Definition of Statistic based on Qe - Prop.6
Bxi=d1;
KXi=(T-1)*sqrt(n)*(Qe-(Bxi/(T-1)));
SIGMA_xi=2*((Vdiag+Vcross)/n);
stat=(sqrt(SIGMA_xi))\KXi;
pv=(1-normcdf(stat,0,1)); 

%Definition of Statistic based on Qa
Bxi_a=d1; 
KXi_a=(T-1)*sqrt(n)*(Qa-(Bxi_a/(T-1)));
SIGMA_xia=2*((Vdiag_a+Vcross_a)/n);
stat_a=(sqrt(SIGMA_xia))\KXi_a;
pv_a=(1-normcdf(stat_a,0,1)); 



