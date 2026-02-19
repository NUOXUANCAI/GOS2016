function X=CondReg(n,T,d1,d2,Z,Zi,F)
% The function computes the regressor for conditional model    

X=zeros(d1+d2,T-1,n);
for i=1:n;
    x1=zeros(d1,T);
    x2=zeros(d2,T);
    for t=2:T;
        %Xt is a matrix s.t. if i==j, Xt(i,j)=Zi,t-1^2, otherwise
        %Xt(i,j)=2*Zi,t-1*Zj,t-1
        Xt=(Z(t-1,:)'*Z(t-1,:))+tril((Z(t-1,:)'*Z(t-1,:)),-1)+triu((Z(t-1,:)'*Z(t-1,:)),1); %Xt - pxp at any t
        Xit=Zi(t-1,:,i)'*Z(t-1,:);   %qxp

        x1(:,t)=[vech(Xt);vec(Xit)]';                                     %d1x1 for any i,t
        x2(:,t)=[kron(F(t,:)',Z(t-1,:)');kron(F(t,:)',Zi(t-1,:,i)')]';    %d2x1 for any i,t
    end;
    x1(:,1)=[]; %d1x(T-1)
    x2(:,1)=[]; %d2x(T-1)
    X(:,:,i)=[x1; x2];
end;
        
