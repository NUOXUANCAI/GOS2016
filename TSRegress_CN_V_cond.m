function [beta1,beta2,res,CN,V,X,R2,R2adj,SYSRisk,IDIRisk]=TSRegress_CN_V_cond(Y,F,Z,Zi,I,Ti,n,T,d1,d2,d)         
%This function estimates the time-varyig betas by TS regression
%see Section 3.2

%regressor matrix - it is defined for any asset i
X=CondReg(n,T,d1,d2,Z,Zi,F);     %X is dx(T-1)xn
    
beta1=zeros(d1,n);      %d1xn
beta2=zeros(d2,n);      %d2xn
CN=zeros(n,1);          %Condition number, nx1
V=zeros(n,d*(d+1)/2);   %V stores vech((inv(Qxi))*Sii*(inv(Qxi)))
res=zeros(T-1,n);       %matrix of residuals,(T-1)xn

R2=zeros(n,1);          %coeffients of determination - Definitions in SM
R2adj=zeros(n,1);       %adjusted coefficients of determination
SYSRisk=zeros(n,1);     %systematic risk
IDIRisk=zeros(n,1);     %idiosyncratic risk

for i=1:n;   
        Ii=repmat(I(:,i),1,d); %indicator matrix Ii, Tx(K+1)
        Xi=Ii.*X(:,:,i)';      %regressor matrix Xi, (T-1)xd
        Qxi=(Xi'*Xi)/Ti(i);    %Qxi, dxd
        %Beta - OLS estimators
        beta=(inv(Qxi))*((Xi'*Y(:,i))/Ti(i));  %dx1
        beta1(:,i)=beta(1:d1,1);        %d1x1
        beta2(:,i)=beta(d1+1:end,1);    %d2x1
        res(:,i)=Y(:,i)-Xi*beta;        %Fitted residual
        
        %Empirical measures of idiosyncratic risk (SM)
        YFit=Xi*beta;
        mYFit=sum(YFit)/Ti(i);
        mY=sum(Y(:,i))/Ti(i);
        TSS=sum((Y(:,i)-repmat(mY,T-1,1)).^2);
        ESS=sum((YFit-repmat(mYFit,T-1,1)).^2);
        RSS=sum(res(:,i).^2);
        R2(i,1)=ESS/TSS;
        R2adj(i,1)=1-(Ti(i)-1)/(Ti(i)-d)*(1-R2(i,1));
        SYSRisk(i,1)=sqrt(ESS/Ti(i));
        IDIRisk(i,1)=sqrt(RSS/Ti(i));    
        
        %Condition number - Greene (2008)
        %X is scaled to have unit column length
        Xtildei=zeros(T-1,d);          
        for s=1:d;
            Xtildei(:,s)=Xi(:,s)/(norm(Xi(:,s)));
        end;        
        Qxtildei=(Xtildei'*Xtildei)/Ti(i);
        %We compute CN_{i} only if Qxtildei is a p.d. matrix
        l=eig(Qxtildei);
        if l>0;
            CN(i,1)=sqrt(max(l)/min(l));                       
        else                     %If the matrix is not p.d.
            CN(i,1)=100;
        end;
        %Definition of Sii matrix and V
        ri=repmat(res(:,i),1,d).*Xi;  
        Sii=(ri'*ri)/Ti(i);             
        Vmat=(inv(Qxi))*Sii*(inv(Qxi)); %dxd 
        V(i,:)=(vech(Vmat))';           %it is a vector (d*(d+1)/2)x1 for any i
end;

 