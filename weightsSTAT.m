function [v,w,V]=weightsSTAT(res,I,Ti,nu,F,tau,n,T,K)
%This function computes the WLS weights used to compute the statistic in
%Section 3.5. 

c=[1;-nu];              %(K+1)x1, nu is the WLS estimator
v=zeros(n,1);           %vector of vii
w=zeros(n,1);           %vector of weigths w
V=zeros(n,(K+1)*(K+1+1)/2);%V stores (vech((\hatQ^-1)Sii(\hatQ^-1)))' for each i -
Tinew=zeros(n,1);       %vector of Ti
X=[ones(T,1),F];        %regressor matrix
Qx=(X'*X)/T;            %Qx, (K+1)x(K+1)   
resORD=sort(abs(res),'descend'); %sort residuals
for i=1:n;
    %trimming on the residuals
    perc=ceil(0.025*Ti(i,1));
    delRes=resORD(1:perc,i);       %residuals that should be deleted, 0.25% of the big res
    for m=1:perc;
        x=find(repmat(delRes(m,1),T,1)==abs(res(:,i)));
        I(x,i)=0;
        res(x,i)=0;        
    end;
    Tinew(i,1)=sum(I(:,i));     %Ti    
    Ii=repmat(I(:,i),1,K+1);    %indicator matrix, Tx(K+1)             
    Xi=Ii.*X;                   %regressor matrix Xi, Tx(K+1)             
    %Definition of Sii, V, v and w
    ri=repmat(res(:,i),1,K+1).*Xi;  
    Sii=(ri'*ri)/Tinew(i,1);             
    Vmat=(inv(Qx))*Sii*(inv(Qx));
    V(i,:)=(vech(Vmat))';
    v(i,1)=tau(i,1)*c'*Vmat*c;
    w(i,1)=1/(v(i,1));    
end;

