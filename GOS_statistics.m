function [Stat,pv1,Stat_a,pv2]=GOS_statistics(nu,beta1,beta2,vecbeta3,Z,X,I,res,Tij,T,d1,d2,d,n,p,q,K)
%This function computes the statistic and its p-values for the two null
%hypothesis in Section 3.5.
%We compute statistics based on Q_{e} and on Q_{a} for the
%time-varying model  when n is fixed and T \rightarrow \infty.

[C,~]=Cmatrix(nu,p,q,K,d1,d2);    %dxd1 matrix
tau=(T-1).*ones(n,n)./Tij;        %tau matrix
ind=isinf(tau);                   %tau is inf if Tij=0
tau(ind)=0;

%Definition of Qe and Qa
E1=[eye(d1);zeros(d2,d1)];  %dxd1
D=duplication(d);           %d*dxd*(d+1)/2
w=zeros(d1,n);              %matrix of weight
wa=zeros(d1,n);             %diagonal elements of weighting matrix
ewe=zeros(n,1);             %it contains e'*wi*e at each row (Qe)
awa=zeros(n,1);             %it contains a'*wa*a at each row (Qa)
for i=1:n
    beta=[beta1(:,i);beta2(:,i)];   %dx1
    e=C'*beta;              %d1x1
    vecV=D*Z(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
    W=invVec(vecV,d,d);     %dxd
    v=tau(i,i)*C'*W*C;      %d1xd1 symmetric matrix
    w(:,i)=diag(v);
    wi=eye(d1,d1)/(diag(w(:,i)));  %d1xd1 weighting matrix  inv(diag(diag(v))      
    ewe(i,1)=e'*wi*e;       %scalar
    
    vecV=D*Z(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
    W=invVec(vecV,d,d);     %dxd
    va=tau(i,i)*E1'*W*E1;   %d1xd1 symmetric matrix
    wa(:,i)=diag(va);       %d1x1
    wai=eye(d1,d1)/(diag(wa(:,i)));%d1xd1 weighting matrix
    awa(i,1)=beta1(:,i)'*wai*beta1(:,i);    %scalar
end;
Qe=sum(ewe)/n;
Qa=sum(awa)/n;

%Definition of statistics
Stat=(T-1)*Qe;
Stat_a=(T-1)*Qa;

%Var-Cov matrix (no thresholding) See Section 3.3 for definitions of
%matrices V,W. We extend these defition to a time-varying setting.
V=zeros(n*d1,n*d1);
V_a=zeros(n*d1,n*d1);
W=zeros(n*d1,n*d1);
W_a=zeros(n*d1,n*d1);
p2=(d1:d1:n*d1)';      %p1,p2: indices to ref to blocks of matrix SIGMA
p1=(1:d1:n*d1)';
B3=zeros(n*d1,K*p);    %matrix B_{n} for time-varying model
for i=1:n;
    beta3=(invVec(vecbeta3(:,i),K*p,d1))';  %d1xKp    
    B3(p1(i):p2(i),:)=beta3;
    
    W(p1(i):p2(i),p1(i):p2(i))=eye(d1,d1)/(diag(w(:,i)));    
    W_a(p1(i):p2(i),p1(i):p2(i))=eye(d1,d1)/(diag(wa(:,i)));
    
    for j=i:n;
        if i~=j && tau(i,j)~=0; 
            Iij=I(:,i).*I(:,j);             
            ri=repmat(Iij.*res(:,i),1,d).*X(:,:,i)';    %(T-1)xd
            rj=repmat(Iij.*res(:,j),1,d).*X(:,:,j)';    %(T-1)xd
            Sij=(ri'*rj)/(Tij(i,j));       %dxd matrix
            
            Qxi=(X(:,:,i)*X(:,:,i)')/Tij(i,i);
            Qxj=(X(:,:,j)*X(:,:,j)')/Tij(j,j);
                        
            SQ_ij=(inv(Qxi))*Sij*(inv(Qxj));
            V(p1(i):p2(i),p1(j):p2(j))=tau(i,i)*tau(j,j)/tau(i,j)*C'*SQ_ij*C;  %d1xd1
            V(p1(j):p2(j),p1(i):p2(i))=V(p1(i):p2(i),p1(j):p2(j));
            
            V_a(p1(i):p2(i),p1(j):p2(j))=tau(i,i)*tau(j,j)/tau(i,j)*E1'*SQ_ij*E1;  %d1xd1
            V_a(p1(j):p2(j),p1(i):p2(i))=V_a(p1(i):p2(i),p1(j):p2(j));
            
        elseif i==j
            vecV=D*Z(i,:)';         %vec(V)=D*vech(V) where D is the duplication matrix - pag.57 MN
            SQ_ii=invVec(vecV,d,d);     
            
            V(p1(i):p2(i),p1(i):p2(i))=tau(i,i)*C'*SQ_ii*C;  %d1xd1            
            V_a(p1(i):p2(i),p1(i):p2(i))=tau(i,i)*E1'*SQ_ii*E1;  %d1xd1
        end;
    end;    
end;

% Definition of p.values
D=sqrtm(V)*(W-W*B3*(inv(B3'*W*B3))*B3'*W)*sqrtm(V); %Extend expression 
egD1=sort(eig(D),'descend'); 
egD2=sort(eig(V_a),'descend'); 
M=100000; %simulation number
CHI1=chi2rnd(1,n*d1-p*K,M);  
CHI2=chi2rnd(1,n*d1,M);
MM1=zeros(M,1);
MM2=zeros(M,1);
for m=1:M;
    for j=1:n*p-d*K;
        MM1(m,1)=MM1(m,1)+egD1(j)*CHI1(j,m);
    end;
    for j=1:n*p;
        MM2(m,1)=MM2(m,1)+egD2(j)*CHI2(j,m);
    end;
    MM1(m,1)=MM1(m,1)/n;
    MM2(m,1)=MM2(m,1)/n;
end;
pv1=sum(MM1>Stat)/M;
pv2=sum(MM2>Stat_a)/M;