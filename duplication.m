function D=duplication(n)
%The function computes the duplication matrix D - for square matrices nxn -
%ref Magnus and Neudecker (2007)

a=tril(ones(n));    %triangular matrix
i=find(a);          %i index of cell not equal zero

a(i)=1:length(i);
a=a+tril(a,-1)';
j=vec(a);

m=n*(n+1)/2;
D=zeros(n*n,m);

for r=1:n*n;
    D(r,j(r))=1;
end;
